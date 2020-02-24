#include "TLorentzVector.h"
#include <vector>
#include <cmath>

using namespace std;

// Constants: could only be changed in the code itself.
const double TINY = 1e-20;

// Vec4 class.
// This class implements four-vectors, in energy-momentum space.
// (But can equally well be used to hold space-time four-vectors.)

class Vec4 {

public:

  // Constructors.
  Vec4(double xIn = 0., double yIn = 0., double zIn = 0., double tIn = 0.)
    : xx(xIn), yy(yIn), zz(zIn), tt(tIn) { }
  Vec4(const Vec4& v) : xx(v.xx), yy(v.yy), zz(v.zz), tt(v.tt) { }
  Vec4& operator=(const Vec4& v) { if (this != &v) { xx = v.xx; yy = v.yy;
    zz = v.zz; tt = v.tt; } return *this; }
  Vec4& operator=(double value) { xx = value; yy = value; zz = value;
    tt = value; return *this; }

  // Member functions for input.
  void reset() {xx = 0.; yy = 0.; zz = 0.; tt = 0.;}
  void p(double xIn, double yIn, double zIn, double tIn)
    {xx = xIn; yy = yIn; zz = zIn; tt = tIn;}
  void p(Vec4 pIn) {xx = pIn.xx; yy = pIn.yy; zz = pIn.zz; tt = pIn.tt;}
  void px(double xIn) {xx = xIn;}
  void py(double yIn) {yy = yIn;}
  void pz(double zIn) {zz = zIn;}
  void e(double tIn) {tt = tIn;}

  // Member functions for output.
  double px() const {return xx;}
  double py() const {return yy;}
  double pz() const {return zz;}
  double e() const {return tt;}
  double& operator[](int i) {
    if      (i == 1) return xx;
    else if (i == 2) return yy;
    else if (i == 3) return zz;
    else             return tt;
  }
  double mCalc() const {double temp = tt*tt - xx*xx - yy*yy - zz*zz;
    return (temp >= 0.) ? sqrt(temp) : -sqrt(-temp);}
  double m2Calc() const {return tt*tt - xx*xx - yy*yy - zz*zz;}
  double pT() const {return sqrt(xx*xx + yy*yy);}
  double pT2() const {return xx*xx + yy*yy;}
  double pAbs() const {return sqrt(xx*xx + yy*yy + zz*zz);}
  double pAbs2() const {return xx*xx + yy*yy + zz*zz;}
  double eT() const {double temp = xx*xx + yy*yy;
    return tt * sqrt( temp / (temp + zz*zz) );}
  double eT2() const {double temp = xx*xx + yy*yy;
    return tt*tt * temp / (temp + zz*zz);}
  double theta() const {return atan2(sqrt(xx*xx + yy*yy), zz);}
  double phi() const {return atan2(yy,xx);}
  double thetaXZ() const {return atan2(xx,zz);}
  double pPos() const {return tt + zz;}
  double pNeg() const {return tt - zz;}
  double rap() const {return 0.5 * log( (tt + zz) / (tt - zz) );}
  double eta() const {double xyz = sqrt(xx*xx + yy*yy + zz*zz);
    return 0.5 * log( (xyz + zz) / (xyz - zz) );}

  // Member functions that perform operations.
  void rescale3(double fac) {xx *= fac; yy *= fac; zz *= fac;}
  void rescale4(double fac) {xx *= fac; yy *= fac; zz *= fac; tt *= fac;}
  void flip3() {xx = -xx; yy = -yy; zz = -zz;}
  void flip4() {xx = -xx; yy = -yy; zz = -zz; tt = -tt;}
  void rot(double thetaIn, double phiIn);
  void rotaxis(double phiIn, double nx, double ny, double nz);
  void rotaxis(double phiIn, const Vec4& n);
  void bst(double betaX, double betaY, double betaZ);
  void bst(double betaX, double betaY, double betaZ, double gamma);
  void bst(const Vec4& pIn);
  void bst(const Vec4& pIn, double mIn);
  void bstback(const Vec4& pIn);
  void bstback(const Vec4& pIn, double mIn);

  // Operator overloading with member functions
  inline Vec4 operator-() {Vec4 tmp; tmp.xx = -xx; tmp.yy = -yy; tmp.zz = -zz;
    tmp.tt = -tt; return tmp;}
  inline Vec4& operator+=(const Vec4& v) {xx += v.xx; yy += v.yy; zz += v.zz;
    tt += v.tt; return *this;}
  inline Vec4& operator-=(const Vec4& v) {xx -= v.xx; yy -= v.yy; zz -= v.zz;
    tt -= v.tt; return *this;}
  inline Vec4& operator*=(double f) {xx *= f; yy *= f; zz *= f;
    tt *= f; return *this;}
  inline Vec4& operator/=(double f) {xx /= f; yy /= f; zz /= f;
    tt /= f; return *this;}
  inline Vec4 operator+(const Vec4& v) const {
    Vec4 tmp = *this; return tmp += v;}
  inline Vec4 operator-(const Vec4& v) const {
    Vec4 tmp = *this; return tmp -= v;}
  inline Vec4 operator*(double f) const {
    Vec4 tmp = *this; return tmp *= f;}
  inline Vec4 operator/(double f) const {
    Vec4 tmp = *this; return tmp /= f;}
  inline double operator*(const Vec4& v) const {
    return tt*v.tt - xx*v.xx - yy*v.yy - zz*v.zz;}

  // Operator overloading with friends
  friend Vec4 operator*(double f, const Vec4& v1);

  // Invariant mass of a pair and its square.
  friend double m(const Vec4& v1, const Vec4& v2);
  friend double m2(const Vec4& v1, const Vec4& v2);

  // Scalar and cross product of 3-vector parts.
  friend double dot3(const Vec4& v1, const Vec4& v2);
  friend Vec4 cross3(const Vec4& v1, const Vec4& v2);

  // theta is polar angle between v1 and v2.
  friend double theta(const Vec4& v1, const Vec4& v2);
  friend double costheta(const Vec4& v1, const Vec4& v2);

  // phi is azimuthal angle between v1 and v2 around z axis.
  friend double phi(const Vec4& v1, const Vec4& v2);
  friend double cosphi(const Vec4& v1, const Vec4& v2);

  // phi is azimuthal angle between v1 and v2 around n axis.
  friend double phi(const Vec4& v1, const Vec4& v2, const Vec4& n);
  friend double cosphi(const Vec4& v1, const Vec4& v2, const Vec4& n);

  // R is distance in cylindrical (y/eta, phi) coordinates.
  friend double RRapPhi(const Vec4& v1, const Vec4& v2);
  friend double REtaPhi(const Vec4& v1, const Vec4& v2);

private:

  // The four-vector data members.
  double xx, yy, zz, tt;
};

// Implementation of operator overloading with friends.
inline Vec4 operator*(double f, const Vec4& v1)
  {Vec4 v = v1; return v *= f;}


//--------------------------------------------------------------------------

// Rotation (simple).

void Vec4::rot(double thetaIn, double phiIn) {

  double cthe = cos(thetaIn);
  double sthe = sin(thetaIn);
  double cphi = cos(phiIn);
  double sphi = sin(phiIn);
  double tmpx =  cthe * cphi * xx -    sphi * yy + sthe * cphi * zz;
  double tmpy =  cthe * sphi * xx +    cphi * yy + sthe * sphi * zz;
  double tmpz = -sthe *        xx +                cthe *        zz;
  xx          = tmpx;
  yy          = tmpy;
  zz          = tmpz;

}

//--------------------------------------------------------------------------

// Azimuthal rotation phi around an arbitrary axis (nz, ny, nz).

void Vec4::rotaxis(double phiIn, double nx, double ny, double nz) {

  double norm = 1./sqrt(nx*nx + ny*ny + nz*nz);
  nx         *= norm;
  ny         *= norm;
  nz         *= norm;
  double cphi = cos(phiIn);
  double sphi = sin(phiIn);
  double comb = (nx * xx + ny * yy + nz * zz) * (1. - cphi);
  double tmpx = cphi * xx + comb * nx + sphi * (ny * zz - nz * yy);
  double tmpy = cphi * yy + comb * ny + sphi * (nz * xx - nx * zz);
  double tmpz = cphi * zz + comb * nz + sphi * (nx * yy - ny * xx);
  xx          = tmpx;
  yy          = tmpy;
  zz          = tmpz;

}

//--------------------------------------------------------------------------

// Azimuthal rotation phi around an arbitrary (3-vector component of) axis.

void Vec4::rotaxis(double phiIn, const Vec4& n) {

  double nx   = n.xx;
  double ny   = n.yy;
  double nz   = n.zz;
  double norm = 1./sqrt(nx*nx + ny*ny + nz*nz);
  nx         *= norm;
  ny          *=norm;
  nz          *=norm;
  double cphi = cos(phiIn);
  double sphi = sin(phiIn);
  double comb = (nx * xx + ny * yy + nz * zz) * (1. - cphi);
  double tmpx = cphi * xx + comb * nx + sphi * (ny * zz - nz * yy);
  double tmpy = cphi * yy + comb * ny + sphi * (nz * xx - nx * zz);
  double tmpz = cphi * zz + comb * nz + sphi * (nx * yy - ny * xx);
  xx          = tmpx;
  yy          = tmpy;
  zz          = tmpz;

}

//--------------------------------------------------------------------------

// Boost (simple).

void Vec4::bst(double betaX, double betaY, double betaZ) {

  double beta2 = betaX*betaX + betaY*betaY + betaZ*betaZ;
  double gamma = 1. / sqrt(1. - beta2);
  double prod1 = betaX * xx + betaY * yy + betaZ * zz;
  double prod2 = gamma * (gamma * prod1 / (1. + gamma) + tt);
  xx += prod2 * betaX;
  yy += prod2 * betaY;
  zz += prod2 * betaZ;
  tt = gamma * (tt + prod1);

}

//--------------------------------------------------------------------------

// Boost (simple, given gamma).

void Vec4::bst(double betaX, double betaY, double betaZ, double gamma) {

  double prod1 = betaX * xx + betaY * yy + betaZ * zz;
  double prod2 = gamma * (gamma * prod1 / (1. + gamma) + tt);
  xx += prod2 * betaX;
  yy += prod2 * betaY;
  zz += prod2 * betaZ;
  tt = gamma * (tt + prod1);

}

//--------------------------------------------------------------------------

// Boost given by a Vec4 p.

void Vec4::bst(const Vec4& pIn) {

  double betaX = pIn.xx / pIn.tt;
  double betaY = pIn.yy / pIn.tt;
  double betaZ = pIn.zz / pIn.tt;
  double beta2 = betaX*betaX + betaY*betaY + betaZ*betaZ;
  double gamma = 1. / sqrt(1. - beta2);
  double prod1 = betaX * xx + betaY * yy + betaZ * zz;
  double prod2 = gamma * (gamma * prod1 / (1. + gamma) + tt);
  xx          += prod2 * betaX;
  yy          += prod2 * betaY;
  zz          += prod2 * betaZ;
  tt           = gamma * (tt + prod1);

}

//--------------------------------------------------------------------------

// Boost given by a Vec4 p and double m.

void Vec4::bst(const Vec4& pIn, double mIn) {

  double betaX = pIn.xx / pIn.tt;
  double betaY = pIn.yy / pIn.tt;
  double betaZ = pIn.zz / pIn.tt;
  double gamma = pIn.tt / mIn;
  double prod1 = betaX * xx + betaY * yy + betaZ * zz;
  double prod2 = gamma * (gamma * prod1 / (1. + gamma) + tt);
  xx          += prod2 * betaX;
  yy          += prod2 * betaY;
  zz          += prod2 * betaZ;
  tt           = gamma * (tt + prod1);

}

//--------------------------------------------------------------------------

// Boost given by a Vec4 p; boost in opposite direction.

void Vec4::bstback(const Vec4& pIn) {

  double betaX = -pIn.xx / pIn.tt;
  double betaY = -pIn.yy / pIn.tt;
  double betaZ = -pIn.zz / pIn.tt;
  double beta2 = betaX*betaX + betaY*betaY + betaZ*betaZ;
  double gamma = 1. / sqrt(1. - beta2);
  double prod1 = betaX * xx + betaY * yy + betaZ * zz;
  double prod2 = gamma * (gamma * prod1 / (1. + gamma) + tt);
  xx          += prod2 * betaX;
  yy          += prod2 * betaY;
  zz          += prod2 * betaZ;
  tt           = gamma * (tt + prod1);

}

//--------------------------------------------------------------------------

// Boost given by a Vec4 p and double m; boost in opposite direction.

void Vec4::bstback(const Vec4& pIn, double mIn) {

  double betaX = -pIn.xx / pIn.tt;
  double betaY = -pIn.yy / pIn.tt;
  double betaZ = -pIn.zz / pIn.tt;
  double gamma = pIn.tt / mIn;
  double prod1 = betaX * xx + betaY * yy + betaZ * zz;
  double prod2 = gamma * (gamma * prod1 / (1. + gamma) + tt);
  xx          += prod2 * betaX;
  yy          += prod2 * betaY;
  zz          += prod2 * betaZ;
  tt           = gamma * (tt + prod1);

}

double pow2(const double& x) {return x*x;}

//--------------------------------------------------------------------------

// The invariant mass of two four-vectors.

double m(const Vec4& v1, const Vec4& v2) {
  double m2 = pow2(v1.tt + v2.tt) - pow2(v1.xx + v2.xx)
     - pow2(v1.yy + v2.yy) - pow2(v1.zz + v2.zz);
  return (m2 > 0.) ? sqrt(m2) : 0.;
}

//--------------------------------------------------------------------------

// The squared invariant mass of two four-vectors.

double m2(const Vec4& v1, const Vec4& v2) {
  double m2 = pow2(v1.tt + v2.tt) - pow2(v1.xx + v2.xx)
     - pow2(v1.yy + v2.yy) - pow2(v1.zz + v2.zz);
  return m2;
}

//--------------------------------------------------------------------------

// The scalar product of two three-vectors.

double dot3(const Vec4& v1, const Vec4& v2) {
  return v1.xx*v2.xx + v1.yy*v2.yy + v1.zz*v2.zz;
}

//--------------------------------------------------------------------------

// The cross product of two three-vectors.

Vec4 cross3(const Vec4& v1, const Vec4& v2) {
  Vec4 v;
  v.xx = v1.yy * v2.zz - v1.zz * v2.yy;
  v.yy = v1.zz * v2.xx - v1.xx * v2.zz;
  v.zz = v1.xx * v2.yy - v1.yy * v2.xx; return v;
}

//--------------------------------------------------------------------------

// Opening angle between two three-vectors.

double theta(const Vec4& v1, const Vec4& v2) {
  double cthe = (v1.xx * v2.xx + v1.yy * v2.yy + v1.zz * v2.zz)
    / sqrt( (v1.xx*v1.xx + v1.yy*v1.yy + v1.zz*v1.zz)
    * (v2.xx*v2.xx + v2.yy*v2.yy + v2.zz*v2.zz) );
  cthe = max(-1., min(1., cthe));
  return acos(cthe);
}

//--------------------------------------------------------------------------

// Cosine of the opening angle between two three-vectors.

double costheta(const Vec4& v1, const Vec4& v2) {
  double cthe = (v1.xx * v2.xx + v1.yy * v2.yy + v1.zz * v2.zz)
    / sqrt( (v1.xx*v1.xx + v1.yy*v1.yy + v1.zz*v1.zz)
    * (v2.xx*v2.xx + v2.yy*v2.yy + v2.zz*v2.zz) );
  cthe = max(-1., min(1., cthe));
  return cthe;
}

//--------------------------------------------------------------------------

// Azimuthal angle between two three-vectors.

double phi(const Vec4& v1, const Vec4& v2) {
  double cphi = (v1.xx * v2.xx + v1.yy * v2.yy) / sqrt( max( TINY,
    (v1.xx*v1.xx + v1.yy*v1.yy) * (v2.xx*v2.xx + v2.yy*v2.yy) ));
  cphi = max(-1., min(1., cphi));
  return acos(cphi);
}

//--------------------------------------------------------------------------

// Cosine of the azimuthal angle between two three-vectors.

double cosphi(const Vec4& v1, const Vec4& v2) {
  double cphi = (v1.xx * v2.xx + v1.yy * v2.yy) / sqrt( max( TINY,
    (v1.xx*v1.xx + v1.yy*v1.yy) * (v2.xx*v2.xx + v2.yy*v2.yy) ));
  cphi = max(-1., min(1., cphi));
  return cphi;
}

//--------------------------------------------------------------------------

// Azimuthal angle between two three-vectors around a third.

double phi(const Vec4& v1, const Vec4& v2, const Vec4& n) {
  double nx = n.xx; double ny = n.yy; double nz = n.zz;
  double norm = 1. / sqrt(nx*nx + ny*ny + nz*nz);
  nx *= norm; ny *=norm; nz *=norm;
  double v1s = v1.xx * v1.xx + v1.yy * v1.yy + v1.zz * v1.zz;
  double v2s = v2.xx * v2.xx + v2.yy * v2.yy + v2.zz * v2.zz;
  double v1v2 = v1.xx * v2.xx + v1.yy * v2.yy + v1.zz * v2.zz;
  double v1n = v1.xx * nx + v1.yy * ny + v1.zz * nz;
  double v2n = v2.xx * nx + v2.yy * ny + v2.zz * nz;
  double cphi = (v1v2 - v1n * v2n) / sqrt( max( TINY,
    (v1s - v1n*v1n) * (v2s - v2n*v2n) ));
  cphi = max(-1., min(1., cphi));
  return acos(cphi);
}

//--------------------------------------------------------------------------

// Cosine of the azimuthal angle between two three-vectors around a third.

double cosphi(const Vec4& v1, const Vec4& v2, const Vec4& n) {
  double nx = n.xx; double ny = n.yy; double nz = n.zz;
  double norm = 1. / sqrt(nx*nx + ny*ny + nz*nz);
  nx *= norm; ny *=norm; nz *=norm;
  double v1s = v1.xx * v1.xx + v1.yy * v1.yy + v1.zz * v1.zz;
  double v2s = v2.xx * v2.xx + v2.yy * v2.yy + v2.zz * v2.zz;
  double v1v2 = v1.xx * v2.xx + v1.yy * v2.yy + v1.zz * v2.zz;
  double v1n = v1.xx * nx + v1.yy * ny + v1.zz * nz;
  double v2n = v2.xx * nx + v2.yy * ny + v2.zz * nz;
  double cphi = (v1v2 - v1n * v2n) / sqrt( max( TINY,
    (v1s - v1n*v1n) * (v2s - v2n*v2n) ));
  cphi = max(-1., min(1., cphi));
  return cphi;
}

//--------------------------------------------------------------------------

// Distance in cylindrical (y, phi) coordinates.

double RRapPhi(const Vec4& v1, const Vec4& v2) {
  double dRap = abs(v1.rap() - v2.rap());
  double dPhi = abs(v1.phi() - v2.phi());
  if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
  return sqrt(dRap*dRap + dPhi*dPhi);
}

//--------------------------------------------------------------------------

// Distance in cylindrical (eta, phi) coordinates.

double REtaPhi(const Vec4& v1, const Vec4& v2) {
  double dEta = abs(v1.eta() - v2.eta());
  double dPhi = abs(v1.phi() - v2.phi());
  if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
  return sqrt(dEta*dEta + dPhi*dPhi);
}

//==========================================================================

// Thrust class.
// This class finds thrust-related properties of an event.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Minimum number of particles to perform study.
const int    NSTUDYMIN    = 2;

// Major not too low or not possible to find major axis.
const double MAJORMIN     = 1e-10;
  
// http://home.fnal.gov/~mrenna/lutp0613man2/node235.html
// Outcome of analysis - thrust/major/minor
// oblateness: eVal2 - eVal3
double eVal1, eVal2, eVal3;
Vec4   eVec1, eVec2, eVec3;

//--------------------------------------------------------------------------

bool thrust(std::vector<TLorentzVector>& event, double betaX=0, double betaY=0, double betaZ=0) {
  if(event.size()<NSTUDYMIN) {
    printf("thrust: too few particles\n");
    return false;
  }

  // Initial values and counters zero.
  eVal1 = eVal2 = eVal3 = 0.;
  eVec1 = eVec2 = eVec3 = 0.;
  int nStudy = 0;
  vector<Vec4> pOrder;
  Vec4 pSum, nRef, pPart, pFull, pMax;

  // Loop over desired particles in the event.
  for (std::vector<TLorentzVector>::size_type i = 0; i < event.size(); ++i) {
    ++nStudy;
    // Store momenta. Use energy component for absolute momentum.
    Vec4 pNow(event[i].Px(),event[i].Py(),event[i].Pz(),event[i].E());
    pNow.e(pNow.pAbs());
    if(betaX!=0 || betaY!=0 || betaZ!=0) {
      pNow.bst(betaX,betaY,betaZ);
    }
    pSum += pNow;
    pOrder.push_back(pNow);
  }

  // Try all combinations of reference vector orthogonal to two particles.
  for (int i1 = 0; i1 < nStudy - 1; ++i1)
  for (int i2 = i1 + 1; i2 < nStudy; ++i2) {
    nRef = cross3( pOrder[i1], pOrder[i2]);
    nRef /= nRef.pAbs();
    pPart = 0.;

    // Add all momenta with sign; two choices for each reference particle.
    for (int i = 0; i < nStudy; ++i) if (i != i1 && i != i2) {
      if (dot3(pOrder[i], nRef) > 0.) pPart += pOrder[i];
      else                            pPart -= pOrder[i];
    }
    for (int j = 0; j < 4; ++j) {
      if      (j == 0) pFull = pPart + pOrder[i1] + pOrder[i2];
      else if (j == 1) pFull = pPart + pOrder[i1] - pOrder[i2];
      else if (j == 2) pFull = pPart - pOrder[i1] + pOrder[i2];
      else             pFull = pPart - pOrder[i1] - pOrder[i2];
      pFull.e(pFull.pAbs());
      if (pFull.e() > pMax.e()) pMax = pFull;
    }
  }

  // Maximum gives thrust axis and value.
  eVal1 = pMax.e() / pSum.e();
  eVec1 = pMax / pMax.e();
  eVec1.e(0.);

  // Subtract momentum along thrust axis.
  double pAbsSum = 0.;
  for (int i = 0; i < nStudy; ++i) {
    pOrder[i] -= dot3( eVec1, pOrder[i]) * eVec1;
    pOrder[i].e(pOrder[i].pAbs());
    pAbsSum += pOrder[i].e();
  }

  // Simpleminded major and minor axes if too little transverse left.
  if (pAbsSum < MAJORMIN * pSum.e()) {
    if ( abs(eVec1.pz()) > 0.5) eVec2 = Vec4( 1., 0., 0., 0.);
    else                        eVec2 = Vec4( 0., 0., 1., 0.);
    eVec2 -= dot3( eVec1, eVec2) * eVec1;
    eVec2 /= eVec2.pAbs();
    eVec3  = cross3( eVec1, eVec2);
    return true;
  }

  // Try all reference vectors orthogonal to one particles.
  pMax = 0.;
  for (int i1 = 0; i1 < nStudy; ++i1) {
    nRef = cross3( pOrder[i1], eVec1);
    nRef /= nRef.pAbs();
    pPart = 0.;

    // Add all momenta with sign; two choices for each reference particle.
    for (int i = 0; i < nStudy; ++i) if (i != i1) {
      if (dot3(pOrder[i], nRef) > 0.) pPart += pOrder[i];
      else                            pPart -= pOrder[i];
    }
    pFull = pPart + pOrder[i1];
    pFull.e(pFull.pAbs());
    if (pFull.e() > pMax.e()) pMax = pFull;
    pFull = pPart - pOrder[i1];
    pFull.e(pFull.pAbs());
    if (pFull.e() > pMax.e()) pMax = pFull;
  }

  // Maximum gives major axis and value.
  eVal2 = pMax.e() / pSum.e();
  eVec2 = pMax / pMax.e();
  eVec2.e(0.);

  // Orthogonal direction gives minor axis, and from there value.
  eVec3 = cross3( eVec1, eVec2);
  pAbsSum = 0.;
  for (int i = 0; i < nStudy; ++i)
    pAbsSum += abs( dot3(eVec3, pOrder[i]) );
  eVal3 = pAbsSum / pSum.e();

   // Done.
  return true;
}
