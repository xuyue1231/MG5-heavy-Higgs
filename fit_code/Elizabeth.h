#include "TROOT.h"
#include "TMath.h"
//#define M_PI 3.14159265359

// individual terms of the Poisson distribution
Double_t SUSYStat_PoissonTerm(Double_t lambda, Int_t i) {
  if(lambda<0) {
    printf("Error: input lambda can not be negative !\n");
    return 0;
  }
  if(i<0) return 0;
  else if(i==0) return exp(-lambda);
  else {
    Double_t t = exp(-lambda);
    for(Int_t k=i;k>=1;k--) t *= lambda/k;
    return t;
  }
}

// forward declaration
Double_t SUSYStat_PoissonCCDF(Double_t lambda, Double_t n, Double_t precision=1e-6);

// the left tail integral (cdf) of the Poisson distribution
Double_t SUSYStat_PoissonCDF(Double_t lambda, Double_t n, Double_t precision=1e-6) {
  if(lambda<0) {
    printf("Error: input lambda can not be negative !\n");
    return 0;
  }
  if(n<0) return 0;
  else if(n>lambda) return 1-SUSYStat_PoissonCCDF(lambda,n);
  Int_t i = Int_t(n);
  Double_t t = SUSYStat_PoissonTerm(lambda,i);
  Double_t p = t;
  while(i>0 && t/p>precision) {
    t *= i/lambda;
    p += t;
    i--;
  }
//   Int_t lo = Int_t(n);
//   Double_t x = n-lo;
//   if(x>0.5) p += SUSYStat_PoissonTerm(lambda,lo+1)*(x-0.5);
//   else p -= SUSYStat_PoissonTerm(lambda,lo)*(0.5-x);
  return p;
}

// the right tail integral (ccdf) of the Poisson distribution
Double_t SUSYStat_PoissonCCDF(Double_t lambda, Double_t n, Double_t precision) {
  if(lambda<0) {
    printf("Error: input lambda can not be negative !\n");
    return 0;
  }
  if(n<0) return 1;
  if(n<=lambda) return 1-SUSYStat_PoissonCDF(lambda,n);
  Int_t i = Int_t(n)+1;
  Double_t t = SUSYStat_PoissonTerm(lambda,i);
  Double_t p = t;
  while(t/p>precision) {
    i++;
    t *= lambda/i;
    p += t;
  }
//   Int_t up = Int_t(n)+1;
//   Double_t x = up-n;
//   if(x>0.5) p += SUSYStat_PoissonTerm(lambda,up-1)*(x-0.5);
//   else p -= SUSYStat_PoissonTerm(lambda,up)*(0.5-x);
  return p;
}

// convert p-value to standard deviation (significance)
Double_t SUSYStat_SigFromP(Double_t p) {
  if(p<=0 || p>=1) {
    printf("Error: input p-value is invalid !\n");
    return 0;
  }
  Double_t pval = (p<0.5?p:1-p);
  Double_t sig = 0;
  if(pval>6.221*10e-16) sig = sqrt(2.)*TMath::ErfcInverse(2*pval);
  else { // asymptotic formula for very small pval
    Double_t u = -2*log(sqrt(2*M_PI)*pval);
    sig = sqrt(u-log(u));
  }
  return (p<0.5?sig:-sig);
}

// convert standard deviation to p-value
Double_t SUSYStat_Pval(Double_t Sig) {
  // asymptotic formula for very large sig
  if(Sig>20) return exp(-Sig*Sig/2)/sqrt(2*M_PI)/Sig;
  else if(Sig<-20) return 1+exp(-Sig*Sig/2)/sqrt(2*M_PI)/Sig;
  else return 0.5*TMath::Erfc(Sig/sqrt(2));
}

// LLR-based discovery significance with sys. error db on b0
// b0 is the expected background, n is the total observed
Double_t SUSYStat_DiscSig_LLR(Double_t b0, Double_t n, Double_t db) {
  if(b0<0) {
    printf("Error: input b0 can not be negative !\n");
    return 0;
  }
  if(db<0) {
    printf("Error: input db can not be negative !\n");
    return 0;
  }
  if(n<=b0) return 0;
  if(db==0) {
    Double_t s = n-b0;
    return sqrt(2*(s+b0)*log(1+s/b0)-2*s);
  }
  Double_t tmp = b0-db*db;
  Double_t b = 0.5*(tmp+sqrt(pow(tmp,2)+4*db*db*n));
  return sqrt(2*n*log(n/b)-n+b-(b-b0)*b0/db/db);
}

// integration-based discovery significance with sys. error db on b0
// b0 is the expected background, n is the total observed
Double_t SUSYStat_DiscSig_Int(Double_t b0, Double_t n, Double_t db, Double_t precision=1e-6) {
  if(b0<0) {
    printf("Error: input b0 can not be negative !\n");
    return 0;
  }
  if(db<0) {
    printf("Error: input db can not be negative !\n");
    return 0;
  }
  if(n<=b0) return 0;
  if(db==0) return SUSYStat_SigFromP(SUSYStat_PoissonCCDF(b0,n));
  Double_t nsig = SUSYStat_SigFromP(precision);
  Double_t p = 0;
  Double_t C = 1/sqrt(2*M_PI)/db;
  Double_t dx = 0.01;
  Double_t x = -nsig;
  while(x<nsig) {
    Double_t b = b0+db*x;
    p += C*exp(-0.5*pow((b-b0)/db,2))*db*dx*SUSYStat_PoissonCCDF(b<0?0:b,n);
    x += dx;
  }
  return SUSYStat_SigFromP(p);
}

// LLR-based exclusion significance with sys. errors db and dmu (relative)
// on the expected background b0 and signal s0, rho is the correlation
// coefficient between db and dmu
Double_t SUSYStat_ExclSig_LLR(Double_t b0, Double_t s0, Double_t n, Double_t db, Double_t dmu, Double_t rho=0) {
  if(b0<0 || s0<0) {
    printf("Error: input b0 or s0 can not be negative !\n");
    return 0;
  }
  if(db<0) {
    printf("Error: input db can not be negative !\n");
    return 0;
  }
  if(dmu<0) {
    printf("Error: input dmu can not be negative !\n");
    return 0;
  }
  if(fabs(rho)>1) {
    printf("Error: input rho is invalid !\n");
    return 0;
  }
  Double_t sign = (n<=s0+b0?1:-1);
  if(db==0 && dmu==0) {
    Double_t term = 0;
    if(n!=0) term = 2*n*log(n/(s0+b0));
    return sign*sqrt(term+2*(s0+b0)-2*n);
  }
  else {
    Double_t d2 = db*db+s0*s0*dmu*dmu+2*rho*db*dmu*s0;
    Double_t lambda = -0.5*(d2-s0-b0)+0.5*sqrt(pow(d2-s0-b0,2)+4*d2*n);
    Double_t term = 0;
    if(n!=0) term = 2*n*log(n/lambda);
    return sign*sqrt(term+lambda-n+(s0+b0-lambda)/d2*(s0+b0));
  }
}

// integration-based exclusion significance with sys. errors db and dmu (relative)
// on the expected background b0 and signal s0, rho is the correlation
// coefficient between db and dmu
Double_t SUSYStat_ExclSig_Int(Double_t b0, Double_t s0, Double_t n, Double_t db, Double_t dmu, Double_t rho=0, Double_t precision=1e-6) {
  if(b0<0 || s0<0) {
    printf("Error: input b0 or s0 can not be negative !\n");
    return 0;
  }
  if(db<0) {
    printf("Error: input db can not be negative !\n");
    return 0;
  }
  if(dmu<0) {
    printf("Error: input dmu can not be negative !\n");
    return 0;
  }
  if(fabs(rho)>1) {
    printf("Error: input rho is invalid !\n");
    return 0;
  }
  if(n>=s0+b0) return 0;
  Double_t nsig = SUSYStat_SigFromP(precision);
  Double_t p = 0;
  if(db==0 && dmu==0) {
    return SUSYStat_SigFromP(SUSYStat_PoissonCDF(s0+b0,n));
  }
  else if(dmu==0) {
    Double_t C = 1/sqrt(2*M_PI)/db;
    Double_t dx = 0.01;
    Double_t x = -nsig;
    while(x<nsig) {
      Double_t b = b0+db*x;
      p += C*exp(-0.5*pow((b-b0)/db,2))*db*dx*SUSYStat_PoissonCDF(s0+(b<0?0:b),n);
      x += dx;
    }
  }
  else if(db==0) {
    Double_t C = 1/sqrt(2*M_PI)/dmu;
    Double_t dx = 0.01;
    Double_t x = -nsig;
    while(x<nsig) {
      Double_t mu = 1+dmu*x;
      p += C*exp(-0.5*pow((mu-1)/dmu,2))*dmu*dx*SUSYStat_PoissonCDF(s0*(mu<0?0:mu)+b0,n);
      x += dx;
    }
  }
  else if(fabs(rho)==1) {
    Double_t C = 1/sqrt(2*M_PI)/dmu;
    Double_t dx = 0.01;
    Double_t x = -nsig;
    while(x<nsig) {
      Double_t mu = dmu*x;
      p += C*exp(-0.5*pow(mu/dmu,2))*dmu*dx*SUSYStat_PoissonCDF((s0+mu*s0<0?0:s0+mu*s0)+(b0+mu*rho*b0*db/dmu<0?0:b0+mu*rho*b0*db/dmu),n);
      x += dx;
    }
  }
  else {
    Double_t Dmu = dmu*sqrt(1-rho*rho);
    Double_t C1 = 1/sqrt(2*M_PI)/db;
    Double_t C2 = 1/sqrt(2*M_PI)/Dmu;
    Double_t dx1 = 0.01;
    Double_t dx2 = 0.01;
    Double_t x1 = -nsig;
    while(x1<nsig) {
      Double_t b = b0+db*x1;
      Double_t mean = 1+rho*(b-b0)*dmu/db;
      Double_t x2 = -nsig;
      while(x2<nsig) {
	Double_t mu = mean+Dmu*x2;
	p += C1*exp(-0.5*pow((b-b0)/db,2))*C2*exp(-0.5*pow((mu-mean)/Dmu,2))*db*dx1*Dmu*dx2*SUSYStat_PoissonCDF(s0*(mu<0?0:mu)+(b<0?0:b),n);
	x2 += dx2;
      }
      x1 += dx1;
    }
  }
  return SUSYStat_SigFromP(p);
}
