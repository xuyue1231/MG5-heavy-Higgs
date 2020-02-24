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
#define MASSLO 75
#define MASSHI 110

#define MASSLO2 70
#define MASSHI2 150

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
Int_t reg;

Float_t mLL;
Float_t ptLL;
TLorentzVector V0, V1, V2;
Bool_t isLL[3];
Float_t tauXX;
Float_t mV[3];
Float_t pV[3];
Float_t ptV[3];
Float_t mHH;
Float_t pHH;
Float_t ptHH;
Float_t EVVV;

void initialize() {
  bn = new binning();

  //bn->add(40,0,2000,"m_{WW}^{truth} [GeV]",&mvv_truth,true);
  //bn->add(150,0,150,"m_{LL} [GeV]",&mLL,true);
  bn->add(50,0,500,"m_{V0} [GeV]",&mV[0],true);
  bn->add(50,0,2000,"p_{T,V0} [GeV]",&ptV[0],true);
  bn->add(50,0,500,"m_{V1} [GeV]",&mV[1],true);
  bn->add(50,0,2000,"p_{T,V1} [GeV]",&ptV[1],true);
  bn->add(50,0,500,"m_{V2} [GeV]",&mV[2],true);
  bn->add(50,0,2000,"p_{T,V2} [GeV]",&ptV[2],true);
  bn->add(40,0,2000,"m_{VV} [GeV]",&mHH,true);
  bn->add(50,0,2000,"p_{T,VV} [GeV]",&ptHH,true);
}

Float_t figSigSF = 1; //sig
Float_t figMax = 1.2;

Bool_t Pass(Int_t ich) {
  //Bool_t cut = mLL>0;
  // Bool_t cut = mLL/GeV>80 && mLL/GeV<100 && mHH>0 && reg==4;
  // Bool_t cut = mLL/GeV>80 && mLL/GeV<100 && mHH>0 && (reg==3 || reg==4);
  Bool_t cut = mLL/GeV>80 && mLL/GeV<100 && mHH>0;
  return cut;
}

void GetVars(Int_t ich) {
  
  TLorentzVector temp_V;
  TLorentzVector temp_VVV;
  temp_VVV.SetPtEtaPhiM(0,0,0,0);

  for(int i=0; i<truthV_n; i++){
    temp_V.SetPtEtaPhiM(0,0,0,0);
      temp_V.SetPtEtaPhiM(truthV_pt[i],truthV_eta[i],truthV_phi[i],truthV_m[i]);
      temp_VVV += temp_V;
  }
  if(truthV_n==3) EVVV = temp_VVV.E()/GeV;

  TLorentzVector vec_l1;
  TLorentzVector vec_l2;
  if(el_n==2 && el_charge[0]*el_charge[1]<0 && mu_n==0) {
    vec_l1.SetPtEtaPhiE(el_pt[0],el_eta[0],el_phi[0],el_E[0]);
    vec_l2.SetPtEtaPhiE(el_pt[1],el_eta[1],el_phi[1],el_E[1]);
    mLL = (vec_l1+vec_l2).M();
    ptLL = (vec_l1+vec_l2).Pt();
  }
  else if(mu_n==2 && mu_charge[0]*mu_charge[1]<0 && el_n==0) {
    vec_l1.SetPtEtaPhiE(mu_pt[0],mu_eta[0],mu_phi[0],mu_E[0]);
    vec_l2.SetPtEtaPhiE(mu_pt[1],mu_eta[1],mu_phi[1],mu_E[1]);
    mLL = (vec_l1+vec_l2).M();
    ptLL = (vec_l1+vec_l2).Pt();
  }
  else {
    mLL = 0;
    ptLL = 0;
  }
  v1 = vec_l1+vec_l2;

  Int_t ifj1(-1);
  Float_t maxPt = -999;
  for(Int_t i=0; i<fatjet_n; i++) {
    if(fatjet_pt[i]/GeV<50) continue;
    if(fatjet_pt[i]>maxPt) {
      maxPt = fatjet_pt[i];
      ifj1 = i;
    }
  }

  Int_t ifj2(-1);
  maxPt = -999;
  for(Int_t i=0; i<fatjet_n; i++) {
    if(fatjet_pt[i]/GeV<50) continue;
    if(i==ifj1) continue;
    if(fatjet_pt[i]>maxPt) {
      maxPt = fatjet_pt[i];
      ifj2 = i;
    }
  }

  Int_t ij[4] = {-1,-1,-1,-1};
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
  maxPt = -999;
  for(Int_t i=0; i<jet_n; i++) {
    if(jet_pt[i]/GeV<30) continue;
    if(i==ij[0] || i==ij[1]) continue;
    if(jet_pt[i]>maxPt) {
      maxPt = jet_pt[i];
      ij[2] = i;
    }
  }
  maxPt = -999;
  for(Int_t i=0; i<jet_n; i++) {
    if(jet_pt[i]/GeV<30) continue;
    if(i==ij[0] || i==ij[1] || i==ij[2]) continue;
    if(jet_pt[i]>maxPt) {
      maxPt = jet_pt[i];
      ij[3] = i;
    }
  }

  v2.SetPxPyPzE(0,0,0,0);
  v3.SetPxPyPzE(0,0,0,0);

  reg = 0;
  if(v2.Pt()==0) {
    if(ptLL/GeV>950 && ifj1>=0 && fatjet_pt[ifj1]/GeV>750 && fatjet_sj_n[ifj1]==2 && fatjet_m_sj1[ifj1]/GeV>70 && fatjet_m_sj1[ifj1]/GeV<150 && fatjet_m_sj2[ifj1]/GeV>70 && fatjet_m_sj2[ifj1]/GeV<150 && fatjet_tau2[ifj1]/fatjet_tau1[ifj1]<0.45) {
      v2.SetPtEtaPhiM(fatjet_pt_sj1[ifj1],fatjet_eta_sj1[ifj1],fatjet_phi_sj1[ifj1],fatjet_m_sj1[ifj1]);
      v3.SetPtEtaPhiM(fatjet_pt_sj2[ifj1],fatjet_eta_sj2[ifj1],fatjet_phi_sj2[ifj1],fatjet_m_sj2[ifj1]);
      tauXX = fatjet_tau2[ifj1]/fatjet_tau1[ifj1];
      reg = 1;
    }
  }
  if(v2.Pt()==0) {
    if(ptLL/GeV>550 && ij[0]>=0 && jet_pt[ij[0]]/GeV>300 && jet_m[ij[0]]/GeV>70 && jet_m[ij[0]]/GeV<150 && ij[1]>=0 && ij[2]>=0 && jet_tau2[ij[0]]/jet_tau1[ij[0]]<0.4) {
      TLorentzVector tmp1; TLorentzVector tmp2; TLorentzVector tmp3;
      tmp1.SetPtEtaPhiM(jet_pt[ij[0]],jet_eta[ij[0]],jet_phi[ij[0]],jet_m[ij[0]]);
      tmp2.SetPtEtaPhiM(jet_pt[ij[1]],jet_eta[ij[1]],jet_phi[ij[1]],jet_m[ij[1]]);
      tmp3.SetPtEtaPhiM(jet_pt[ij[2]],jet_eta[ij[2]],jet_phi[ij[2]],jet_m[ij[2]]);
      TLorentzVector tmp23 = tmp2+tmp3;
      if(tmp23.M()/GeV>70 && tmp23.M()/GeV<110 && tmp23.Pt()/GeV>150 && tmp23.DeltaR(tmp1)<tmp23.DeltaR(v1) && tmp23.DeltaR(tmp1)<tmp1.DeltaR(v1) && (tmp1+tmp23).Pt()/GeV>550) {
    	v2 = tmp1;
    	v3 = tmp23;
    	tauXX = jet_tau2[ij[0]]/jet_tau1[ij[0]];
  	reg = 2;
      }
    }
  }

  if(v2.Pt()==0) {
    if(ij[0]>=0 && jet_pt[ij[0]]/GeV>700 && jet_m[ij[0]]/GeV>70 && jet_m[ij[0]]/GeV<150 && jet_tau2[ij[0]]/jet_tau1[ij[0]]<0.60 && ij[1]>=0 && ij[2]>=0 && ptLL/GeV>300) {
      TLorentzVector tmp1; TLorentzVector tmp2; TLorentzVector tmp3;
      tmp1.SetPtEtaPhiM(jet_pt[ij[0]],jet_eta[ij[0]],jet_phi[ij[0]],jet_m[ij[0]]);
      tmp2.SetPtEtaPhiM(jet_pt[ij[1]],jet_eta[ij[1]],jet_phi[ij[1]],jet_m[ij[1]]);
      tmp3.SetPtEtaPhiM(jet_pt[ij[2]],jet_eta[ij[2]],jet_phi[ij[2]],jet_m[ij[2]]);
      TLorentzVector tmp23 = tmp2+tmp3;
      if(tmp23.M()/GeV>75 && tmp23.M()/GeV<115 && tmp23.Pt()/GeV>50 && tmp23.DeltaR(v1)<tmp23.DeltaR(tmp1) && tmp23.DeltaR(v1)<tmp1.DeltaR(v1) && (tmp23+v1).Pt()/GeV>700) {
    	v2 = tmp1;
    	v3 = tmp23;
    	tauXX = jet_tau2[ij[0]]/jet_tau1[ij[0]];
    	reg = 3;
      }
    }
  }
  if(v2.Pt()==0) {
    if(ij[0]>=0 && jet_pt[ij[0]]/GeV>700 && jet_m[ij[0]]/GeV>70 && jet_m[ij[0]]/GeV<150 && ij[1]>=0 && jet_pt[ij[1]]/GeV>250 && jet_m[ij[1]]/GeV>70 && jet_m[ij[1]]/GeV<150 && jet_tau2[ij[0]]/jet_tau1[ij[0]]<0.52 && jet_tau2[ij[1]]/jet_tau1[ij[1]]<0.52 && ptLL/GeV>300) {
      TLorentzVector tmp1; TLorentzVector tmp2;
      tmp1.SetPtEtaPhiM(jet_pt[ij[0]],jet_eta[ij[0]],jet_phi[ij[0]],jet_m[ij[0]]);
      tmp2.SetPtEtaPhiM(jet_pt[ij[1]],jet_eta[ij[1]],jet_phi[ij[1]],jet_m[ij[1]]);
      if(tmp2.DeltaR(v1)<tmp2.DeltaR(tmp1) && tmp2.DeltaR(v1)<v1.DeltaR(tmp1) && (v1+tmp2).Pt()/GeV>700) {
	v2 = tmp1;
	v3 = tmp2;
	tauXX = jet_tau2[ij[0]]/jet_tau1[ij[0]];
	reg = 4;
      }

    }
  }

  V0.SetPxPyPzE(0,0,0,0); V1.SetPxPyPzE(0,0,0,0); V2.SetPxPyPzE(0,0,0,0);
  isLL[0] = isLL[1] = isLL[2] = false;
  mV[0] = mV[1] = mV[2] = 0;
  pV[0] = pV[1] = pV[2] = 0;
  ptV[0] = ptV[1] = ptV[2] = 0;
  mHH = pHH = ptHH = 0;
  if(reg==1 || reg==2) {
    V0 = v1; V1 = v2; V2 = v3;
    isLL[0] = true;
    mV[0] = V0.M(); mV[1] = V1.M(); mV[2] = V2.M();
    pV[0] = V0.P(); pV[1] = V1.P(); pV[2] = V2.P();
    ptV[0] = V0.Pt(); ptV[1] = V1.Pt(); ptV[2] = V2.Pt();
    mHH = (V1+V2).M();
    pHH = (V1+V2).P();
    ptHH = (V1+V2).Pt();
  }
  else if(reg==3 || reg==4) {
    V0 = v2; V1 = v1; V2 = v3;
    isLL[1] = true;
    mV[0] = V0.M(); mV[1] = V1.M(); mV[2] = V2.M();
    pV[0] = V0.P(); pV[1] = V1.P(); pV[2] = V2.P();
    ptV[0] = V0.Pt(); ptV[1] = V1.Pt(); ptV[2] = V2.Pt();
    mHH = (V1+V2).M();
    pHH = (V1+V2).P();
    ptHH = (V1+V2).Pt();
  }

  else if(reg==6){
  	V0 = v1;
    mV[0] = V0.M(); pV[0] = V0.P(); ptV[0] = V0.Pt();
    isLL[0] = true;
    mHH = v2.M();
    pHH = v2.P();
    ptHH = v2.Pt();
  }
   if(ich==0) mc_event_weight *= lumi*0.00062418374e03/0.00062418374e04; //600,1000,6000
   /////////////////////////////quadrant_1                             //mass,fw,fww
   // if(ich==0)      mc_event_weight *= lumi*3.95981917e-02/3.95981917e-01; // 600,30,30
 
  else if(ich==1) mc_event_weight *= lumi*6.83279757e01/68.72184;
  else if(ich==2) mc_event_weight *= lumi*6.87452922e01/68.70894;
  else if(ich==3) mc_event_weight *= lumi*6.154065e02/3.03691e-07;
  else if(ich==4) mc_event_weight *= lumi*6.150717e02/3.03093e-07;
  else if(ich==5) mc_event_weight *= lumi*0.01454e03/1470.53;
}

void plot() {
  char str[200];

  initialize();

  TString title_sig = "hsig_fw_70_fww_m70";
  TChain* ch[6] = {0,0,0,0,0,0};

  ch[0] = new TChain("tau");
   // ch[0]->Add("../../../rootfiles/test600/ntuple_2lep/S2_red_15.root"); //600,1000,6000
   // ch[0]->Add("../../../rootfiles/test600/ntuple_2lep/S2_red_380.root"); //600,1000,6000
   ch[0]->Add("../../../rootfiles/rhoH_005/ntuple_2lep/S2_red_1.root"); //600,1500,0
  ////////////////////////quadrant_1                           //mass,fw,fww
  // ch[0]->Add("../../rootfiles/ntuple_2lep/S2_red_68.root"); //600,30,30
  Init(ch[0]);

  ch[1] = new TChain("tau");
  ch[1]->Add("../../../rootfiles/ntuple_merged2/zee_red.root");
  Init(ch[1]);

  ch[2] = new TChain("tau");
  ch[2]->Add("../../../rootfiles/ntuple_merged2/zmm_red.root");
  Init(ch[2]);

  ch[3] = new TChain("tau");
  ch[3]->Add("../../../rootfiles/ntuple_merged2/zee_vjj_red.root");
  Init(ch[3]);

  ch[4] = new TChain("tau");
  ch[4]->Add("../../../rootfiles/ntuple_merged2/zmm_vjj_red.root");
  Init(ch[4]);

  ch[5] = new TChain("tau");
  ch[5]->Add("../../../rootfiles/ntuple_merged2/zll_vv4j_red.root");
  Init(ch[5]);
  
  //////////////////////////////////
  ///////////////////////////////////////
  TH1F* ht1 = new TH1F("ht1","",100,0,1); 
  TH1F* ht2 = new TH1F("ht2","",100,0,1);

  TH1F* hsig = new TH1F(title_sig,"",40,0,2000);
  TH1F* hbkg1 = new TH1F("hbkg1","zll",40,0,2000);

  // sig,bkg
  TH1F* h1[50]; TH1F* h2[50];
  
  for(Int_t i=0; i<bn->nvar; i++) {
    sprintf(str,"h1_var%d",i+1);
    h1[i]  = new TH1F(str,"",bn->nbin[i],bn->xlo[i],bn->xhi[i]);
    sprintf(str,"h2_var%d",i+1);
    h2[i]  = new TH1F(str,"",bn->nbin[i],bn->xlo[i],bn->xhi[i]);
  }

  TH1F* hh1 = new TH1F("hh1","",200,0,1); 
  TH1F* hh2 = new TH1F("hh2","",200,0,1);
  
  TH1F* hz = new TH1F("hz","",1000,0,10); 
  TH1F* h_EVVV = new TH1F("hz","",100,0,7000); 

  TH1F* hdr1 = new TH1F("hdr1","",200,-5,5); 
  TH1F* hdr2 = new TH1F("hdr2","",200,-5,5);
  
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

      Double_t wt = mc_event_weight;
      n_eve[ich][0] += wt;

      if(truthV_n==3) h_EVVV->Fill(EVVV);

      if(Pass(ich)) {
	n_eve[ich][1] += wt;
     
     if(ich==0) hsig->Fill(mHH/GeV,wt);
    else  hbkg1->Fill(mHH/GeV,wt);

	for(Int_t i=0; i<bn->nvar; i++) {
	  if(ich==0) h1[i]->Fill(bn->getVal(i),wt);
	  else h2[i]->Fill(bn->getVal(i),wt);
	}

	if(ich==0) ht1->Fill(tauXX);
	else ht2->Fill(tauXX);

      }
    }
  }
  printf("Sig: %.1f %.3f\n",n_eve[0][0],n_eve[0][1]);
  printf("Bkg1: %.1f %.3f\n",n_eve[1][0],n_eve[1][1]);
  printf("Bkg2: %.1f %.3f\n",n_eve[2][0],n_eve[2][1]);
  printf("Bkg3: %.1f %.3f\n",n_eve[3][0],n_eve[3][1]);
  printf("Bkg4: %.1f %.3f\n",n_eve[4][0],n_eve[4][1]);
  printf("Bkg5: %.1f %.3f\n",n_eve[5][0],n_eve[5][1]);
  printf("significance: %0.3f\n", significance(n_eve[1][1]+n_eve[2][1]+n_eve[3][1]+n_eve[4][1]+n_eve[5][1], n_eve[0][1], 0) );

  TCanvas* c01 = new TCanvas("c01","c01",600,600);
  ht1->SetLineColor(kRed);
  ht2->SetLineColor(kBlue);
  ht1->Scale(ht2->Integral()/ht1->Integral());
  SetMax(ht1,ht2,1.1);
  ht1->Draw();
  ht2->Draw("same");

  TLegend* lg = new TLegend(0.62,0.74,0.90,0.89,"");
  if(figSigSF==1) sprintf(str,"Signal");
  else sprintf(str,"Signal#times%1.0f",figSigSF);
  lg->AddEntry(h1[0],str,"F");
  lg->AddEntry(h2[0],"Background","F");
  lg->SetBorderSize(0);
  lg->SetMargin(0.25);
  lg->SetFillColor(kWhite);

  THStack* hsk[50];
  TCanvas* cv[50];

  for(Int_t i=0; i<bn->nvar; i++) {
    h1[i]->Scale(figSigSF);
    h1[i]->SetLineColor(kRed);
    h1[i]->SetFillColor(kRed);
    h2[i]->SetLineColor(kBlue);
    h2[i]->SetFillColor(kBlue);

    sprintf(str,"hsk_var%d",i+1);
    hsk[i] = new THStack(str,"");
    if(bn->titleX[i]=="m_{VV} [GeV]") {
      hsk[i]->Add(h2[i]);
      hsk[i]->Add(h1[i]);
    }
    else {
      hsk[i]->Add(h1[i]);
      hsk[i]->Add(h2[i]);
    }
    
    sprintf(str,"c%d",i+1);
    cv[i] = new TCanvas(str,str,600,600);
    hsk[i]->Draw("hist");
    hsk[i]->GetXaxis()->SetTitle(bn->titleX[i].c_str());
    sprintf(str,"Events/( %4.2f )",h1[i]->GetBinWidth(1));
    hsk[i]->GetYaxis()->SetTitle(str);
    lg->Draw("same");

    // sprintf(str,"c%d.eps",i+1);
    // cv[i]->SaveAs(str);
  }
  TCanvas* cvvv = new TCanvas();
  h_EVVV->Draw();
// TFile f1("mHH_600_hist_bkg_2lep_lumi3000.root","recreate");
// TFile f1("mHH_600_hist_sig_2lep.root","recreate");
// TFile f1("mHH_600_hist_sig_2lep.root","update");
// f1.cd();
// hsig->Write();
// hbkg1->Write();
// f1.Close();
}
