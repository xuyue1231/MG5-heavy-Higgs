#include "TROOT.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "THStack.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "../reader.h"
#define GeV 1000
#define MASSLO 70
#define MASSHI 110
#define lumi 300 //fb-1

void calibrate_jet() {
  Double_t thr[7] = {30,50,80,120,200,400,800};
  Double_t sf[3][7] = {
    {1.00113, 1.0039,  1.00243, 1.00129, 1.00032, 1.00126, 1.00133},
    {1.05701, 1.03567, 1.02281, 1.01503, 1.00959, 1.00691, 1.00511},
    {1.04688, 1.02907, 1.021  , 1.01858, 1.02032, 1.     , 1.     }
  };

  for(Int_t i=0; i<jet_n; i++) {
    Double_t SF = 1.0;
    if(fabs(jet_eta[i])<1.7) {
      for(Int_t j=6; j>=0; j--) {
	if(jet_pt[i]/GeV>=thr[j]) {
	  SF = 1/sf[0][j]; break;
	}
      }
    }
    else if(fabs(jet_eta[i])<3.2) {
      for(Int_t j=6; j>=0; j--) {
	if(jet_pt[i]/GeV>=thr[j]) {
	  SF = 1/sf[1][j]; break;
	}
      }
    }
    else {
      for(Int_t j=6; j>=0; j--) {
	if(jet_pt[i]/GeV>=thr[j]) {
	  SF = 1/sf[2][j]; break;
	}
      }
    } 
    jet_E[i] *= SF;
    jet_m[i] *= SF;
    jet_pt[i] *= SF;
  }  
}

void calibrate_fatjet() {
  Double_t thr[7] = {40,50,80,120,200,400,800};
  Double_t sf[3][7] = {
    {0.835528, 0.876777, 0.914269, 0.940094, 0.962477, 0.981316, 0.990554},
    {0.959362, 0.968022, 0.976831, 0.983729, 0.989864, 0.995393, 0.997024},
    {0.979653, 0.984607, 0.990887, 0.99717,  1.00537 , 1.      , 1.      }
  };

  for(Int_t i=0; i<fatjet_n; i++) {
    Double_t SF = 1.0;
    if(fabs(fatjet_eta[i])<1.7) {
      for(Int_t j=6; j>=0; j--) {
	if(fatjet_pt[i]/GeV>=thr[j]) {
	  SF = 1/sf[0][j]; break;
	}
      }
    }
    else if(fabs(fatjet_eta[i])<3.2) {
      for(Int_t j=6; j>=0; j--) {
	if(fatjet_pt[i]/GeV>=thr[j]) {
	  SF = 1/sf[1][j]; break;
	}
      }
    }
    else {
      for(Int_t j=6; j>=0; j--) {
	if(fatjet_pt[i]/GeV>=thr[j]) {
	  SF = 1/sf[2][j]; break;
	}
      }
    } 
    fatjet_E[i] *= SF;
    fatjet_m[i] *= SF;
    fatjet_pt[i] *= SF;
    fatjet_m_sj1[i] *= SF;
    fatjet_pt_sj1[i] *= SF;
    fatjet_m_sj2[i] *= SF;
    fatjet_pt_sj2[i] *= SF;
  }
}

void SetMax(TH1* h1, TH1* h2, Double_t scale=1.0) {
  h1->SetMaximum(scale*TMath::Max(h1->GetMaximum(),h2->GetMaximum()));
  h2->SetMaximum(scale*TMath::Max(h1->GetMaximum(),h2->GetMaximum()));
}

Float_t DR(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2) {
  Float_t dphi = fabs(phi1-phi2);
  if(dphi>M_PI) dphi = 2*M_PI-dphi;
  Float_t deta = fabs(eta1-eta2);
  return sqrt(dphi*dphi+deta*deta);
}

Double_t NormPhi(Double_t phi) {
  if(phi>=M_PI) return phi - 2*M_PI;
  else if(phi<-M_PI) return phi + 2*M_PI;
  else return phi;
}
Double_t significance(Double_t b0, Double_t s0, Double_t db) {
  if(db==0) return sqrt(2*(s0+b0)*log(1+s0/b0)-2*s0);
  else {
    Double_t tmp = b0-db*db;
    Double_t b = 0.5*(tmp+sqrt(pow(tmp,2)+4*db*db*(b0+s0)));
    return sqrt(2*(b0+s0)*log((b0+s0)/b)-(b0+s0)+b-(b-b0)*b0/db/db);
  }
}

class binning {
public:
  Int_t nbin[50];
  Float_t xlo[50];
  Float_t xhi[50];
  std::string titleX[50];
  Float_t* var1[50];
  Int_t* var2[50];
  Bool_t MeVtoGeV[50];
  Int_t nvar;
  
  binning();
  virtual ~binning();

  void add(Int_t nbin_, Float_t xlo_, Float_t xhi_, const char* titleX_, Int_t* var_);
  void add(Int_t nbin_, Float_t xlo_, Float_t xhi_, const char* titleX_, Float_t* var_, Bool_t MeVtoGeV_);
  Float_t getVal(Int_t i);
};

binning::binning() {
  nvar = 0;
  for(Int_t i=0; i<50; i++) {
    nbin[i] = 1; xlo[i] = 0; xhi[i] = 1; var1[i] = 0; var2[i] = 0; MeVtoGeV[i] = 0;
  }
}

binning::~binning() {}

void binning::add(Int_t nbin_, Float_t xlo_, Float_t xhi_, const char* titleX_, Int_t* var_) {
  if(nvar>=0 && nvar<50) {
    nbin[nvar] = nbin_; xlo[nvar] = xlo_; xhi[nvar] = xhi_;
    titleX[nvar] = titleX_; var2[nvar] = var_;
    nvar++;
  }
}

void binning::add(Int_t nbin_, Float_t xlo_, Float_t xhi_, const char* titleX_, Float_t* var_, Bool_t MeVtoGeV_) {
  if(nvar>=0 && nvar<50) {
    nbin[nvar] = nbin_; xlo[nvar] = xlo_; xhi[nvar] = xhi_;
    titleX[nvar] = titleX_; var1[nvar] = var_; MeVtoGeV[nvar] = MeVtoGeV_;
    nvar++;
  }
}

Float_t binning::getVal(Int_t i) {
  Float_t tmp = -999999;
  if(i>=0  && i<nvar) {
    if(var1[i]) tmp = MeVtoGeV[i] ? *var1[i]/GeV : *var1[i];
    else if(var2[i]) tmp = *var2[i];
  }
  if(tmp >= xhi[i]) tmp = xhi[i]*0.999999;
  if(tmp < xlo[i]) tmp = xlo[i];
  return tmp;
}

binning* bn(0);

// new vars:
TLorentzVector v1, v2, v3;
TLorentzVector vec_nu;
Bool_t isFatJet_v3;
Float_t tau21_v3;
Int_t reg;

Float_t mLL;
Float_t ptLL;
Float_t ptL3;
TLorentzVector V0, V1, V2;
Bool_t isLL[3];
Float_t mV[3];
Float_t pV[3];
Float_t ptV[3];
Float_t mHH;
Float_t pHH;
Float_t ptHH;
Float_t mT;

void initialize() {
  bn = new binning();

  //bn->add(40,0,2000,"m_{WW}^{truth} [GeV]",&mvv_truth,true);
  bn->add(100,50,150,"m_{LL} [GeV]",&mLL,true);
  bn->add(50,0,1500,"p_{T,V0} [GeV]",&ptV[0],true);
  bn->add(50,0,500,"m_{V2} [GeV]",&mV[2],true);
  bn->add(25,0,1000,"p_{T,V2} [GeV]",&ptV[2],true);
  bn->add(36,200,2000,"m_{VV} [GeV]",&mHH,true);
  bn->add(50,0,2000,"p_{T,VV} [GeV]",&ptHH,true);
  bn->add(50,0,1500,"MET [GeV]",&MET_et,true);
  bn->add(25,0,250,"m_{T} [GeV]",&mT,true);
  bn->add(30,0,1500,"p_{T,l3} [GeV]",&ptL3,true);
}

Float_t figSigSF = 1; //sig
Float_t figMax = 1.2;

Bool_t Pass(Int_t ich) {
  //Bool_t cut = mLL>0;
  //Bool_t cut = mLL/GeV>80 && mLL/GeV<100 && mHH>0 && ptL3/GeV>50 && (reg==3 || reg==4 || reg==7);
  //Bool_t cut = mLL/GeV>80 && mLL/GeV<100 && mHH>0 && ptL3/GeV>50 && (reg==1 || reg==2 || reg==6);
  Bool_t cut = mLL/GeV>80 && mLL/GeV<100 && mHH>0 && ptL3/GeV>50;
  return cut;
}

void GetVars(Int_t ich) {
  TLorentzVector vec_l1;
  TLorentzVector vec_l2;
  TLorentzVector vec_l3;
  if(el_n==3 && mu_n==0) {
    vec_l1.SetPtEtaPhiE(el_pt[0],el_eta[0],el_phi[0],el_E[0]);
    vec_l2.SetPtEtaPhiE(el_pt[1],el_eta[1],el_phi[1],el_E[1]);
    vec_l3.SetPtEtaPhiE(el_pt[2],el_eta[2],el_phi[2],el_E[2]);
    if(el_charge[0]*el_charge[1]<0 && el_charge[0]*el_charge[2]<0) {
      if(vec_l1.DeltaR(vec_l3)<vec_l1.DeltaR(vec_l2)) {
	TLorentzVector tmp = vec_l2; vec_l2 = vec_l3; vec_l3 = tmp;
      }
      mLL = (vec_l1+vec_l2).M();
      ptLL = (vec_l1+vec_l2).Pt();
    }
    else if(el_charge[0]*el_charge[1]<0 && el_charge[1]*el_charge[2]<0) {
      if(vec_l2.DeltaR(vec_l3)<vec_l1.DeltaR(vec_l2)) {
	TLorentzVector tmp = vec_l1; vec_l1 = vec_l3; vec_l3 = tmp;
      }
      mLL = (vec_l1+vec_l2).M();
      ptLL = (vec_l1+vec_l2).Pt();
    }
    else if(el_charge[0]*el_charge[2]<0 && el_charge[1]*el_charge[2]<0) {
      if(vec_l1.DeltaR(vec_l3)<vec_l2.DeltaR(vec_l3)) {
	TLorentzVector tmp = vec_l2; vec_l2 = vec_l3; vec_l3 = tmp;
      }
      else {
	TLorentzVector tmp = vec_l1; vec_l1 = vec_l3; vec_l3 = tmp;
      }
      mLL = (vec_l1+vec_l2).M();
      ptLL = (vec_l1+vec_l2).Pt();
    }
    else {
      mLL = 0;
      ptLL = 0;
    }
  }
  else if(mu_n==3 && el_n==0) {
    vec_l1.SetPtEtaPhiE(mu_pt[0],mu_eta[0],mu_phi[0],mu_E[0]);
    vec_l2.SetPtEtaPhiE(mu_pt[1],mu_eta[1],mu_phi[1],mu_E[1]);
    vec_l3.SetPtEtaPhiE(mu_pt[2],mu_eta[2],mu_phi[2],mu_E[2]);
    if(mu_charge[0]*mu_charge[1]<0 && mu_charge[0]*mu_charge[2]<0) {
      if(vec_l1.DeltaR(vec_l3)<vec_l1.DeltaR(vec_l2)) {
	TLorentzVector tmp = vec_l2; vec_l2 = vec_l3; vec_l3 = tmp;
      }
      mLL = (vec_l1+vec_l2).M();
      ptLL = (vec_l1+vec_l2).Pt();
    }
    else if(mu_charge[0]*mu_charge[1]<0 && mu_charge[1]*mu_charge[2]<0) {
      if(vec_l2.DeltaR(vec_l3)<vec_l1.DeltaR(vec_l2)) {
	TLorentzVector tmp = vec_l1; vec_l1 = vec_l3; vec_l3 = tmp;
      }
      mLL = (vec_l1+vec_l2).M();
      ptLL = (vec_l1+vec_l2).Pt();
    }
    else if(mu_charge[0]*mu_charge[2]<0 && mu_charge[1]*mu_charge[2]<0) {
      if(vec_l1.DeltaR(vec_l3)<vec_l2.DeltaR(vec_l3)) {
	TLorentzVector tmp = vec_l2; vec_l2 = vec_l3; vec_l3 = tmp;
      }
      else {
	TLorentzVector tmp = vec_l1; vec_l1 = vec_l3; vec_l3 = tmp;
      }
      mLL = (vec_l1+vec_l2).M();
      ptLL = (vec_l1+vec_l2).Pt();
    }
    else {
      mLL = 0;
      ptLL = 0;
    }
  }
  else if(el_n==2 && el_charge[0]*el_charge[1]<0 && mu_n==1) {
    vec_l1.SetPtEtaPhiE(el_pt[0],el_eta[0],el_phi[0],el_E[0]);
    vec_l2.SetPtEtaPhiE(el_pt[1],el_eta[1],el_phi[1],el_E[1]);
    vec_l3.SetPtEtaPhiE(mu_pt[0],mu_eta[0],mu_phi[0],mu_E[0]);
    mLL = (vec_l1+vec_l2).M();
    ptLL = (vec_l1+vec_l2).Pt();
  }
  else if(mu_n==2 && mu_charge[0]*mu_charge[1]<0 && el_n==1) {
    vec_l1.SetPtEtaPhiE(mu_pt[0],mu_eta[0],mu_phi[0],mu_E[0]);
    vec_l2.SetPtEtaPhiE(mu_pt[1],mu_eta[1],mu_phi[1],mu_E[1]);
    vec_l3.SetPtEtaPhiE(el_pt[0],el_eta[0],el_phi[0],el_E[0]);
    mLL = (vec_l1+vec_l2).M();
    ptLL = (vec_l1+vec_l2).Pt();
  }
  else {
    mLL = 0;
    ptLL = 0;
  }
  v1 = vec_l1+vec_l2;
  ptL3 = vec_l3.Pt();
  
  vec_nu.SetXYZT(0,0,0,0);
  if(vec_l3.Pt()>0) {
    Float_t phiX = atan2(MET_ey,MET_ex);
    Float_t dphi = fabs(vec_l3.Phi() - phiX);
    if(dphi>M_PI) dphi = 2*M_PI-dphi;
    mT = sqrt(2*vec_l3.Pt()*MET_et*(1-cos(dphi)));
    Float_t pzMiss = 0;
    if(mT/GeV<80.4) {
      Float_t ptX = vec_l3.Px()*MET_ex+vec_l3.Py()*MET_ey;
      Float_t tmp = pow(80.4*GeV,2)+2*ptX;
      Float_t pzMiss1 = (vec_l3.Pz()*tmp-vec_l3.E()*sqrt(pow(tmp,2)-pow(2*vec_l3.Pt()*MET_et,2)))/2/vec_l3.Pt()/vec_l3.Pt();
      Float_t pzMiss2 = (vec_l3.Pz()*tmp+vec_l3.E()*sqrt(pow(tmp,2)-pow(2*vec_l3.Pt()*MET_et,2)))/2/vec_l3.Pt()/vec_l3.Pt();
      // TLorentzVector tmp1;
      // TLorentzVector tmp2;
      // tmp1.SetXYZM(MET_ex,MET_ey,pzMiss1,0);
      // tmp2.SetXYZM(MET_ex,MET_ey,pzMiss2,0);
      // if(vec_l3.DeltaR(tmp1)<vec_l3.DeltaR(tmp2)) pzMiss = pzMiss2;
      // else pzMiss = pzMiss1;
      
      // if(ich==0) {
      // 	Float_t pzMiss_truth(0);
      // 	for(Int_t i=0; i<truthV_n; i++) {
      // 	  if(abs(truthV_pdgId[i])==24) {
      // 	    if(abs(truthV_dau1_pdgId[i])==12 || abs(truthV_dau1_pdgId[i])==14 || abs(truthV_dau1_pdgId[i])==16) {
      // 	      pzMiss_truth = truthV_dau1_pt[i]*sinh(truthV_dau1_eta[i]);
      // 	      break;
      // 	    }
      // 	    if(abs(truthV_dau2_pdgId[i])==12 || abs(truthV_dau2_pdgId[i])==14 || abs(truthV_dau2_pdgId[i])==16) {
      // 	      pzMiss_truth = truthV_dau2_pt[i]*sinh(truthV_dau2_eta[i]);
      // 	      break;
      // 	    }
      // 	  }
      // 	}
      // 	pzMiss = fabs(pzMiss1-pzMiss_truth)<fabs(pzMiss2-pzMiss_truth) ? pzMiss1 : pzMiss2;
      // }
      // else pzMiss = fabs(pzMiss1)<fabs(pzMiss2) ? pzMiss1 : pzMiss2;

      pzMiss = fabs(pzMiss1)<fabs(pzMiss2) ? pzMiss1 : pzMiss2;
    }
    else {
      pzMiss = vec_l3.Pz()*MET_et/vec_l3.Pt();
    }
    vec_nu.SetXYZM(MET_ex,MET_ey,pzMiss,0);
  }
  v2 = vec_l3+vec_nu;

  Int_t ifj1(-1);
  Float_t maxPt = -999;
  for(Int_t i=0; i<fatjet_n; i++) {
    if(fatjet_pt[i]/GeV<50) continue;
    if(fatjet_pt[i]>maxPt) {
      maxPt = fatjet_pt[i];
      ifj1 = i;
    }
  }

  Int_t ij[2] = {-1,-1};
  maxPt = -999;
  for(Int_t i=0; i<jet_n; i++) {
    if(jet_pt[i]/GeV<30) continue;
    if(jet_pt[i]>maxPt) {
      maxPt = jet_pt[i];
      ij[0] = i;
    }
  }
  maxPt = -999;
  for(Int_t i=0; i<jet_n; i++) {
    if(jet_pt[i]/GeV<30) continue;
    if(i==ij[0]) continue;
    if(jet_pt[i]>maxPt) {
      maxPt = jet_pt[i];
      ij[1] = i;
    }
  }
  
  isFatJet_v3 = false;
  tau21_v3 = -999;
  v3.SetXYZT(0,0,0,0);

  reg = 0;
  if(v1.Pt()>0 && v2.Pt()>0) {
     if(v3.Pt()==0) {
       if(v2.Pt()/GeV>600 && ij[0]>=0 && jet_m[ij[0]]/GeV>60 && jet_m[ij[0]]/GeV<160 && jet_tau2[ij[0]]/jet_tau1[ij[0]]<0.6) {
	 TLorentzVector tmp;
	 tmp.SetPtEtaPhiM(jet_pt[ij[0]],jet_eta[ij[0]],jet_phi[ij[0]],jet_m[ij[0]]);
	 if(v1.DeltaR(tmp)<v1.DeltaR(v2) && v1.DeltaR(tmp)<tmp.DeltaR(v2) && (v1+tmp).Pt()/GeV>600) {
	   v3 = tmp;
	   isFatJet_v3 = true;
	   tau21_v3 = jet_tau2[ij[0]]/jet_tau1[ij[0]];
	   reg = 3;
	 }
       }
     }
    if(v3.Pt()==0) {
      if(v2.Pt()/GeV>600 && ifj1>=0 && fatjet_m[ifj1]/GeV>70 && fatjet_m[ifj1]/GeV<140 && fatjet_tau2[ifj1]/fatjet_tau1[ifj1]<0.5) {
        TLorentzVector tmp;
      	tmp.SetPtEtaPhiM(fatjet_pt[ifj1],fatjet_eta[ifj1],fatjet_phi[ifj1],fatjet_m[ifj1]);
	if(v1.DeltaR(tmp)<v1.DeltaR(v2) && v1.DeltaR(tmp)<tmp.DeltaR(v2) && (v1+tmp).Pt()/GeV>600) {
	  v3 = tmp;
	  isFatJet_v3 = true;
	  tau21_v3 = fatjet_tau2[ifj1]/fatjet_tau1[ifj1];
	  reg = 4;
        }
      }
    }
    if(v3.Pt()==0) {
      if(v2.Pt()/GeV>600 && ij[0]>=0 && ij[1]>=0) {
        TLorentzVector tmp1; TLorentzVector tmp2;
	tmp1.SetPtEtaPhiM(jet_pt[ij[0]],jet_eta[ij[0]],jet_phi[ij[0]],jet_m[ij[0]]);
	tmp2.SetPtEtaPhiM(jet_pt[ij[1]],jet_eta[ij[1]],jet_phi[ij[1]],jet_m[ij[1]]);
	TLorentzVector tmp = tmp1+tmp2;
	if(tmp.M()/GeV>60 && tmp.M()/GeV<120 && v1.DeltaR(tmp)<v1.DeltaR(v2) && v1.DeltaR(tmp)<tmp.DeltaR(v2) && (v1+tmp).Pt()/GeV>600) {
	  v3 = tmp;
	  isFatJet_v3 = false;
	  reg = 7;
        }
      }
    }

    if(v3.Pt()==0) {
      if(v1.Pt()/GeV>600 && ij[0]>=0 && jet_m[ij[0]]/GeV>60 && jet_m[ij[0]]/GeV<160 && jet_tau2[ij[0]]/jet_tau1[ij[0]]<0.6) {
	TLorentzVector tmp;
	tmp.SetPtEtaPhiM(jet_pt[ij[0]],jet_eta[ij[0]],jet_phi[ij[0]],jet_m[ij[0]]);
	if(v2.DeltaR(tmp)<v1.DeltaR(v2) && v2.DeltaR(tmp)<tmp.DeltaR(v1) && (v2+tmp).Pt()/GeV>600) {
	  v3 = tmp;
	  isFatJet_v3 = true;
	  tau21_v3 = jet_tau2[ij[0]]/jet_tau1[ij[0]];
	  reg = 1;
	}
      }
    }
    if(v3.Pt()==0) {
      if(v1.Pt()/GeV>600 && ifj1>=0 && fatjet_m[ifj1]/GeV>70 && fatjet_m[ifj1]/GeV<140 && fatjet_tau2[ifj1]/fatjet_tau1[ifj1]<0.5) {
        TLorentzVector tmp;
      	tmp.SetPtEtaPhiM(fatjet_pt[ifj1],fatjet_eta[ifj1],fatjet_phi[ifj1],fatjet_m[ifj1]);
	if(v2.DeltaR(tmp)<v1.DeltaR(v2) && v2.DeltaR(tmp)<tmp.DeltaR(v1) && (v2+tmp).Pt()/GeV>600) {
	  v3 = tmp;
	  isFatJet_v3 = true;
	  tau21_v3 = fatjet_tau2[ifj1]/fatjet_tau1[ifj1];
	  reg = 2;
        }
      }
    }
    if(v3.Pt()==0) {
      if(v1.Pt()/GeV>600 && ij[0]>=0 && ij[1]>=0) {
        TLorentzVector tmp1; TLorentzVector tmp2;
	tmp1.SetPtEtaPhiM(jet_pt[ij[0]],jet_eta[ij[0]],jet_phi[ij[0]],jet_m[ij[0]]);
	tmp2.SetPtEtaPhiM(jet_pt[ij[1]],jet_eta[ij[1]],jet_phi[ij[1]],jet_m[ij[1]]);
	TLorentzVector tmp = tmp1+tmp2;
	if(tmp.M()/GeV>60 && tmp.M()/GeV<120 && v2.DeltaR(tmp)<v1.DeltaR(v2) && v2.DeltaR(tmp)<tmp.DeltaR(v1) && (v2+tmp).Pt()/GeV>600) {
	  v3 = tmp;
	  isFatJet_v3 = false;
	  reg = 6;
        }
      }
    }
  }

  V0.SetPxPyPzE(0,0,0,0); V1.SetPxPyPzE(0,0,0,0); V2.SetPxPyPzE(0,0,0,0);
  isLL[0] = isLL[1] = isLL[2] = false;
  mV[0] = mV[1] = mV[2] = 0;
  pV[0] = pV[1] = pV[2] = 0;
  ptV[0] = ptV[1] = ptV[2] = 0;
  mHH = pHH = ptHH = 0;
  if(reg==3 || reg==4 || reg==7) {
    V0 = v2; V1 = v1; V2 = v3;
    isLL[1] = true;
    mV[0] = V0.M(); mV[1] = V1.M(); mV[2] = V2.M();
    pV[0] = V0.P(); pV[1] = V1.P(); pV[2] = V2.P();
    ptV[0] = V0.Pt(); ptV[1] = V1.Pt(); ptV[2] = V2.Pt();
    mHH = (V1+V2).M();
    pHH = (V1+V2).P();
    ptHH = (V1+V2).Pt();
  }
  else if(reg==1 || reg==2 || reg==6) {
    V0 = v1; V1 = v2; V2 = v3;
    isLL[0] = true;
    mV[0] = V0.M(); mV[1] = V1.M(); mV[2] = V2.M();
    pV[0] = V0.P(); pV[1] = V1.P(); pV[2] = V2.P();
    ptV[0] = V0.Pt(); ptV[1] = V1.Pt(); ptV[2] = V2.Pt();
    mHH = (V1+V2).M();
    pHH = (V1+V2).P();
    ptHH = (V1+V2).Pt();
  }
  else if(v1.Pt()>0 && v2.Pt()>0 && v3.Pt()>0) {
    if(v1.DeltaR(v3)<v2.DeltaR(v3)) {
      V0 = v2; V1 = v1; V2 = v3;
      isLL[1] = true;
    }
    else {
      V0 = v1; V1 = v2; V2 = v3;
      isLL[0] = true;
    }
    mV[0] = V0.M(); mV[1] = V1.M(); mV[2] = V2.M();
    pV[0] = V0.P(); pV[1] = V1.P(); pV[2] = V2.P();
    ptV[0] = V0.Pt(); ptV[1] = V1.Pt(); ptV[2] = V2.Pt();
    mHH = (V1+V2).M();
    pHH = (V1+V2).P();
    ptHH = (V1+V2).Pt();
  }

  if(ich==0)      mc_event_weight *= lumi*0.000234555557e03/0.000234555557e04; // 600, 800, 800
  // if(ich==0)      mc_event_weight *= lumi*4.09409661e-02/4.09409661e-01; // 600, 800, 800
  
  else if(ich==1) mc_event_weight *= lumi*0.701193227/0.34599; //zll_wmlv
  else if(ich==2) mc_event_weight *= lumi*1.41261797/0.67463; //zll_wplv
  else if(ich==3) mc_event_weight *= lumi*0.004726e03/467.61446;//zll_wlv_vjj
  else if(ich==4) mc_event_weight *= lumi*0.014e03/24624.98383; //ttV_3lep
}

void plot() {
  char str[200];

  initialize();
  TString title_sig = "hsig_fw_30_fww_m30";
  TChain* ch[6] = {0,0,0,0,0,0};

  ch[0] = new TChain("tau");
  // ch[0]->Add("../../../rootfiles/test600/ntuple_3lep/S2_red_380.root");
  // ch[0]->Add("../../../rootfiles/rhoH_005/ntuple_3lep/S2_red_265.root");
  ch[0]->Add("../../../rootfiles/rhoH_005/ntuple_3lep/S2_red_1.root");
  Init(ch[0]);
  
  ch[1] = new TChain("tau");
  ch[1]->Add("../../../rootfiles/ntuple_merged2/zll_wmlv_red.root");
  Init(ch[1]);

  ch[2] = new TChain("tau");
  ch[2]->Add("../../../rootfiles/ntuple_merged2/zll_wplv_red.root");
  Init(ch[2]);

  ch[3] = new TChain("tau");
  ch[3]->Add("../../../rootfiles/ntuple_merged2/zll_wlv_vjj_red.root");
  Init(ch[3]);

  ch[4] = new TChain("tau");
  ch[4]->Add("../../../rootfiles/ntuple_merged2/ttV_3lep_red.root");
  Init(ch[4]);

  // sig,bkg
  TH1F* h1[50]; TH1F* h2[50]; TH1F* h3[50];
  
  for(Int_t i=0; i<bn->nvar; i++) {
    sprintf(str,"h1_var%d",i+1);
    h1[i]  = new TH1F(str,"",bn->nbin[i],bn->xlo[i],bn->xhi[i]);
    sprintf(str,"h2_var%d",i+1);
    h2[i]  = new TH1F(str,"",bn->nbin[i],bn->xlo[i],bn->xhi[i]);
    sprintf(str,"h3_var%d",i+1);
    h3[i]  = new TH1F(str,"",bn->nbin[i],bn->xlo[i],bn->xhi[i]);
  }

  TH1F* hh1 = new TH1F("hh1","",100,0,1); 
  TH1F* hh2 = new TH1F("hh2","",100,0,1);

  TH1F* hsig = new TH1F(title_sig,"",40,0,2000);
  TH1F* hbkg1 = new TH1F("hbkg1","zll",40,0,2000);
  TH1F* hbkg2 = new TH1F("hbkg2","ttV",40,0,2000);
  
  TH1F* hz = new TH1F("hz","",1000,0,10); 
  TH1F* hmis = new TH1F("hmis","",200,-400,400); 

  Double_t n_eve[6][5] = {
    {0,0,0,0,0},
    {0,0,0,0,0},
    {0,0,0,0,0},
    {0,0,0,0,0},
    {0,0,0,0,0},
    {0,0,0,0,0}
  };
  
  for(Int_t ich=0; ich<6; ich++) {
    if(!ch[ich]) continue;
    Int_t numberOfEntries = (Int_t)ch[ich]->GetEntries();
    // Loop over all events
    printf(" %d entries to be processed\n",numberOfEntries);
    for(Int_t entry = 0; entry < numberOfEntries; entry++) {
    // Load selected branches with data from specified event
      ch[ich]->GetEntry(entry);
      if((entry+1)%10000000==0) printf(" %d entries processed\n",entry+1);

      calibrate_jet();
      calibrate_fatjet();
      GetVars(ich);

      Float_t pxMiss_truth(0);
      Float_t pzMiss_truth(0);
      
      Double_t wt = mc_event_weight;
      n_eve[ich][0] += wt;

      if(Pass(ich)) {

	n_eve[ich][1] += wt;

	if(ich==0) hsig->Fill(mHH/GeV,wt);
	else if(ich==1||ich==2||ich==3) hbkg1->Fill(mHH/GeV,wt);
	else if(ich==4) hbkg2->Fill(mHH/GeV,wt);



	for(Int_t i=0; i<bn->nvar; i++) {
	  if(ich==0) h1[i]->Fill(bn->getVal(i),wt);
	  else if(ich==1 || ich==2) h2[i]->Fill(bn->getVal(i),wt);
	  else if(ich==3||ich==4) h3[i]->Fill(bn->getVal(i),wt);
	}

	if(ich==0) {
	  if(isFatJet_v3) hh1->Fill(tau21_v3,wt);
	  if(isFatJet_v3) n_eve[0][2] += wt;
	}
	else {
	  if(isFatJet_v3) hh2->Fill(tau21_v3,wt);
	  if(isFatJet_v3) n_eve[1][2] += wt;
	  hz->Fill(wt);
	}

	if(ich==0 && pzMiss_truth!=0) {
	  hmis->Fill(vec_nu.Pz()/GeV-pzMiss_truth/GeV);
	  //hmis->Fill(vec_nu.Px()/GeV-pxMiss_truth/GeV);
	}
      }
    }
  }
  printf("Sig: %.5f %.1f\n",n_eve[0][0],n_eve[0][1]);
  printf("Bkg1: %.5f %.1f\n",n_eve[1][0],n_eve[1][1]);
  printf("Bkg2: %.5f %.1f\n",n_eve[2][0],n_eve[2][1]);
  printf("Bkg3: %.5f %.1f\n",n_eve[3][0],n_eve[3][1]);
  printf("Bkg4: %.5f %.1f\n",n_eve[4][0],n_eve[4][1]);
  printf("HasFatJet: %.1f(sig) %.1f(bkg)\n",n_eve[0][2],n_eve[1][2]);
  printf("significance: %0.3f\n", significance(n_eve[1][1]+n_eve[2][1]+n_eve[3][1]+n_eve[4][1]+n_eve[5][1], n_eve[0][1], 0) );

  TLegend* lg = new TLegend(0.62,0.70,0.90,0.89,"");
  if(figSigSF==1) sprintf(str,"Signal");
  else sprintf(str,"Sig.#times%.1f",figSigSF);
  lg->AddEntry(h1[0],"Signal","F");
  lg->AddEntry(h2[0],"WZ+jets","F");
  lg->AddEntry(h3[0],"other","F");
  lg->SetBorderSize(0);
  lg->SetMargin(0.25);
  lg->SetFillColor(kWhite);

  TCanvas* c01 = new TCanvas("c01","c01",600,600);
  hh1->SetLineColor(kRed);
  hh2->SetLineColor(kBlue);
  hh1->Scale(hh2->Integral()/hh1->Integral());
  SetMax(hh1,hh2,1.1);
  hh1->Draw();
  hh2->Draw("same");
  
  TCanvas* c02 = new TCanvas("c02","c02",600,600);
  hz->Draw();
    
  TCanvas* c03 = new TCanvas("c03","c03",600,600);
  hmis->Draw();
  
  THStack* hsk[50];
  TCanvas* cv[50];

  for(Int_t i=0; i<bn->nvar; i++) {
    h1[i]->Scale(figSigSF);
    h1[i]->SetLineColor(kRed);
    h1[i]->SetFillColor(kRed);
    h2[i]->SetLineColor(kBlue);
    h2[i]->SetFillColor(kBlue);
    h3[i]->SetLineColor(kGreen);
    h3[i]->SetFillColor(kGreen);

    sprintf(str,"hsk_var%d",i+1);
    hsk[i] = new THStack(str,"");
    if(bn->titleX[i]=="m_{VV} [GeV]") {
      hsk[i]->Add(h3[i]);
      hsk[i]->Add(h2[i]);
      hsk[i]->Add(h1[i]);
      // hsk[i]->Add(h1[i]);
      // hsk[i]->Add(h3[i]);
      // hsk[i]->Add(h2[i]);
    }
    else {
      hsk[i]->Add(h1[i]);
      hsk[i]->Add(h3[i]);
      hsk[i]->Add(h2[i]);
    }
    
    sprintf(str,"c%d",i+1);
    cv[i] = new TCanvas(str,str,800,600);
    hsk[i]->Draw("hist");
    hsk[i]->GetXaxis()->SetTitle(bn->titleX[i].c_str());
    sprintf(str,"Events/( %4.2f )",h1[i]->GetBinWidth(1));
    hsk[i]->GetYaxis()->SetTitle(str);
    lg->Draw("same");

    // sprintf(str,"c%d.eps",i+1);
    // cv[i]->SaveAs(str);
}
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Float_t num_s = 0.5;
  Double_t num_s = n_eve[0][1];
  Double_t num_bkg = n_eve[1][1]+n_eve[2][1]+n_eve[3][1];
  Double_t num_sig_bkg = num_s+n_eve[1][1]+n_eve[2][1]+n_eve[3][1];
  TH1F* p_sig = new TH1F("p_sig","",50, 0,50);
  TH1F* p_bkg = new TH1F("p_bkg","",50, 0,50);
  TH1F* p_sig_bkg = new TH1F("p_sig_bkg","",60,0,60);
  for(int i=0; i<30000; i++){
  	p_sig->Fill(gRandom->PoissonD(num_s));
  	p_bkg->Fill(gRandom->PoissonD(num_bkg));
  	p_sig_bkg->Fill(gRandom->PoissonD(num_sig_bkg));
  }
Double_t med_sig_bkg, med_bkg, q;
q=0.5;
p_sig_bkg->GetQuantiles(1, &med_sig_bkg, &q);
p_bkg->GetQuantiles(1, &med_bkg, &q);
cout<<"med sig+bkg:"<<med_sig_bkg<<endl;
cout<<"med bkg:"<<med_bkg<<endl;
 TCanvas *c1 = new TCanvas("c1","demo quantiles",600,600);
   // p_sig_bkg->Draw();
   p_bkg->Draw();
   // show the quantiles in the bottom pad
///////////p value and significance: med_sig_bkg in bkg poisson
TAxis* axis_bkg = p_bkg->GetXaxis();
int bin_bkg = axis_bkg->FindBin(med_sig_bkg+1);
Double_t inte_bkg = p_bkg->Integral(bin_bkg, 50);
// inte_bkg -= p_bkg->GetBinContent(bin_bkg)*(med_sig_bkg-axis_bkg->GetBinLowEdge(bin_bkg))/axis_bkg->GetBinWidth(bin_bkg);
// Double_t p_value = inte_bkg/p_bkg->Integral();
cout<<"bin:"<<bin_bkg<<";low edge:"<<axis_bkg->GetBinLowEdge(bin_bkg)<<";content:"<<p_bkg->GetBinContent(bin_bkg)<<endl;
//////// CL:1-alpha      med_bkg in sig+bkg poisson
TAxis* axis_sig_bkg = p_sig_bkg->GetXaxis();
int bin_sig_bkg = axis_sig_bkg->FindBin(med_bkg);
Double_t inte_sig_bkg = p_sig_bkg->Integral(1, bin_sig_bkg);
// inte_sig_bkg -= p_sig_bkg->GetBinContent(bin_sig_bkg)*(axis_sig_bkg->GetBinUpEdge(bin_sig_bkg)-med_bkg)/axis_sig_bkg->GetBinWidth(bin_sig_bkg);
// Double_t alpha = inte_sig_bkg/p_sig_bkg->Integral();
Double_t alpha = inte_sig_bkg/30000.;
Double_t CL = 1-alpha;
cout<<"CL:"<<CL<<endl;

/////////直接用公式算significance和CL////// CL=1-p-value
Double_t pval = ROOT::Math::chisquared_cdf_c(2.*(num_bkg+num_s), 2*(int(med_bkg)+1));
// p-value = ROOT::Math::normal_cdf_c(d_signifi, 1);
Double_t d_signifi = ROOT::Math::normal_quantile_c(pval,1);

Double_t d_CL = ROOT::Math::chisquared_cdf(2.*(num_bkg+num_s), 2*(int(med_bkg)+1));

cout<<"p value:"<<pval<<endl;
cout<<"fomual signi:"<<d_signifi<<endl;
cout<<"fo CL:"<<d_CL<<endl;

cout<<hsig->Integral()<<endl;
cout<<hbkg1->Integral()<<endl;
cout<<hbkg2->Integral()<<endl;

// TFile f1("mHH_600_hist_bkg_3lep_lumi3000.root","recreate");
// // TFile f1("mHH_600_hist_sig.root","recreate");
// // TFile f1("mHH_600_hist_sig.root","update");
// f1.cd();
// // hsig->Write();
// hbkg1->Write();
// hbkg2->Write();
// f1.Close();


 }
