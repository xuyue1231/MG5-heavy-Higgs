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
TLorentzVector v1, v3;
Bool_t isFatJet_v3;
Float_t tau21_v3;
Int_t reg;

Float_t mLL;
Float_t ptLL;
Float_t ptL1;
Float_t ptL2;
Float_t dphiLL;
Float_t mV3;
Float_t pV3;
Float_t ptV3;
Float_t d_eta_jj;
Bool_t hasBjet;

void initialize() {
  bn = new binning();

  bn->add(100,0,1000,"m_{LL} [GeV]",&mLL,true);
  bn->add(50,0,1500,"p_{T,LL} [GeV]",&ptLL,true);
  bn->add(50,0,1500,"p_{T,L1} [GeV]",&ptL1,true);
  bn->add(50,0,500,"p_{T,L2} [GeV]",&ptL2,true);
  bn->add(20,0,acos(-1),"#Delta#phi_{LL}",&dphiLL,false);
  //bn->add(50,0,500,"m_{V(jj)} [GeV]",&mV3,true);
  bn->add(20,0,250,"m_{V(jj)} [GeV]",&mV3,true);
  bn->add(100,0,2000,"m_{LL} [GeV]",&mLL,true);
  bn->add(50,0,2500,"p_{T,V(jj)} [GeV]",&ptV3,true);
  bn->add(50,0,1500,"MET [GeV]",&MET_et,true);
  bn->add(20,0,5,"#Delta#eta_{jj}",&d_eta_jj,true);
}

Float_t figSigSF = 1; //sig
Float_t figMax = 1.2;

Bool_t Pass(Int_t ich) { //cut1 ptL1>300 ,, previous ptL1>200 1,2
                        // change ptL1>400 to ptL2>450 2
  // cut on mLL performed well in reg 1
  Bool_t cut1 = mLL/GeV>300 && ptV3/GeV>0 && MET_et/GeV>100 && ptLL/GeV>100 && ptL1/GeV>300 && ptL2/GeV>50 && dphiLL>2.0 && !hasBjet && (reg==1 || reg==2); // without mass window cut 
  // Bool_t cut1 = mV3/GeV>60 && mV3/GeV<150 && mLL/GeV>300 && ptV3/GeV>0 && MET_et/GeV>100 && ptLL/GeV>100 && ptL1/GeV>300 && ptL2/GeV>50 && dphiLL>2.0 && !hasBjet && (reg==1 || reg==2); // add cut on mLL
  Bool_t cut3 = mLL/GeV>400 && ptV3/GeV>0 && MET_et/GeV>100 && ptLL/GeV>100 && ptL1/GeV>450 && ptL2/GeV>50 && dphiLL>1.6 && (reg==3);
  return (cut1 || cut3);
  // return (cut1);
}

void GetVars(Int_t ich) {
  TLorentzVector vec_l1;
  TLorentzVector vec_l2;
  //// two electrons ???
  if(el_n==1 && mu_n==1 && el_charge[0]*mu_charge[0]>0) {
    vec_l1.SetPtEtaPhiE(el_pt[0],el_eta[0],el_phi[0],el_E[0]);
    vec_l2.SetPtEtaPhiE(mu_pt[0],mu_eta[0],mu_phi[0],mu_E[0]);
    mLL = (vec_l1+vec_l2).M();
    ptLL = (vec_l1+vec_l2).Pt();
    dphiLL = vec_l1.DeltaR(vec_l2);
  }
  else if(mu_n==2 && mu_charge[0]*mu_charge[1]>0 && el_n==0) {
    vec_l1.SetPtEtaPhiE(mu_pt[0],mu_eta[0],mu_phi[0],mu_E[0]);
    vec_l2.SetPtEtaPhiE(mu_pt[1],mu_eta[1],mu_phi[1],mu_E[1]);
    mLL = (vec_l1+vec_l2).M();
    ptLL = (vec_l1+vec_l2).Pt();
    dphiLL = vec_l1.DeltaR(vec_l2);
  }
  else if(el_n==2 && el_charge[0]*el_charge[1]>0 && mu_n==0) {
    vec_l1.SetPtEtaPhiE(el_pt[0],el_eta[0],el_phi[0],el_E[0]);
    vec_l2.SetPtEtaPhiE(el_pt[1],el_eta[1],el_phi[1],el_E[1]);
    mLL = (vec_l1+vec_l2).M();
    ptLL = (vec_l1+vec_l2).Pt();
    dphiLL = vec_l1.DeltaR(vec_l2);
  }
  else {
    mLL = 0;
    ptLL = 0;
    dphiLL = 0;
  }
  ptL1 = ptL2 = 0;
  if(vec_l1.Pt()>vec_l2.Pt()) {
    ptL1 = vec_l1.Pt(); ptL2 = vec_l2.Pt();
  }
  else {
    ptL1 = vec_l2.Pt(); ptL2 = vec_l1.Pt();
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

  hasBjet = false;
  for(Int_t i=0; i<jet_n; i++) {
    if(jet_pt[i]/GeV>30 && fabs(jet_eta[i])<2.5 && jet_isBtagged[i]) {
      hasBjet = true;
      break;
    }
  } 
  
  isFatJet_v3 = false;
  tau21_v3 = -999;
  v3.SetXYZT(0,0,0,0);

  reg = 0;
  if(v1.Pt()>0) {
    if(v3.Pt()==0) {
      //if(ij[0]>=0 && jet_tau2[ij[0]]/jet_tau1[ij[0]]<0.7) { // zll_wlv_red.root has no tauXX
      if(ij[0]>=0) {
	TLorentzVector tmp;
	tmp.SetPtEtaPhiM(jet_pt[ij[0]],jet_eta[ij[0]],jet_phi[ij[0]],jet_m[ij[0]]); //main
	if(tmp.Pt()/GeV>400) {
	  v3 = tmp;
	  isFatJet_v3 = false;
	  reg = 1;
	}
      }
    }
    if(v3.Pt()==0) {
      if(ifj1>=0 && fatjet_tau2[ifj1]/fatjet_tau1[ifj1]<0.6) {
        TLorentzVector tmp;
      	tmp.SetPtEtaPhiM(fatjet_pt[ifj1],fatjet_eta[ifj1],fatjet_phi[ifj1],fatjet_m[ifj1]);
	if(tmp.Pt()/GeV>100) {
	  v3 = tmp;
	  isFatJet_v3 = true;
	  tau21_v3 = fatjet_tau2[ifj1]/fatjet_tau1[ifj1];
	  reg = 2;
        }
      }
    }
    if(v3.Pt()==0) {
      if(ij[0]>=0 && ij[1]>=0) {
        TLorentzVector tmp1; TLorentzVector tmp2;
	tmp1.SetPtEtaPhiM(jet_pt[ij[0]],jet_eta[ij[0]],jet_phi[ij[0]],jet_m[ij[0]]);
	tmp2.SetPtEtaPhiM(jet_pt[ij[1]],jet_eta[ij[1]],jet_phi[ij[1]],jet_m[ij[1]]);
	TLorentzVector tmp = tmp1+tmp2;
  d_eta_jj = fabs(jet_eta[ij[0]]-jet_eta[ij[1]])*1000.;
  // if((tmp1+tmp2).Pt()/GeV>250 && d_eta_jj<1.5*1000.) {
	if((tmp1+tmp2).Pt()/GeV>0) {
	  v3 = tmp;
	  isFatJet_v3 = false;
	  reg = 3;
        }
      }
    }
  }

  mV3 = v3.M();
  pV3 = v3.P();
  ptV3 = v3.Pt();

  if(ich==0)      mc_event_weight *= lumi*0.1638887/8.19444; // 600
  
  else if(ich==1) mc_event_weight *= lumi*3.134693e02/1.55753e-07; //zll_wlv
  else if(ich==2) mc_event_weight *= lumi*0.014e03/24624.98383; //ttV_3lep
  else if(ich==3) mc_event_weight *= lumi*0.010286e03/2057.61; //ttV_2l_ss
  else if(ich==4) mc_event_weight *= lumi*0.004726e03/467.61446;//zll_wlv_vjj
  else if(ich==5) mc_event_weight *= lumi*3.5841227/358.384; //www_2l_ss
}

void plot() {
  char str[200];

  initialize();
  TChain* ch[6] = {0,0,0,0,0,0};

  ch[0] = new TChain("tau");
  ch[0]->Add("../../../rootfiles/ntupleForWWW/sig_2lepWW_red_001.root"); // 1000，1000
  // ch[0]->Add("../../../rootfiles/rhoH_005/ntuple_www/S2_red_265.root");// 700，700
  Init(ch[0]);
  
  ch[1] = new TChain("tau");
  ch[1]->Add("../../../rootfiles/ntuple_merged/zll_wlv_red.root"); // no filter
  Init(ch[1]);

  ch[2] = new TChain("tau");
  ch[2]->Add("../../../rootfiles/ntuple_merged2/ttV_3lep_red.root");
  Init(ch[2]);

  ch[3] = new TChain("tau");
  ch[3]->Add("../../../rootfiles/ntuple_merged2/ttV_2l_ss_red.root");
  Init(ch[3]);

  ch[4] = new TChain("tau");
  ch[4]->Add("../../../rootfiles/ntuple_merged2/zll_wlv_vjj_red.root");
  Init(ch[4]);

  ch[5] = new TChain("tau");
  ch[5]->Add("../../../rootfiles/ntuple_merged2/www_2l_ss_red.root");
  Init(ch[5]);

  // sig,bkg
  TH1F* h1[50]; TH1F* h2[50]; TH1F* h3[50]; TH1F* h4[50];
  
  for(Int_t i=0; i<bn->nvar; i++) {
    sprintf(str,"h1_var%d",i+1);
    h1[i]  = new TH1F(str,"",bn->nbin[i],bn->xlo[i],bn->xhi[i]);
    sprintf(str,"h2_var%d",i+1);
    h2[i]  = new TH1F(str,"",bn->nbin[i],bn->xlo[i],bn->xhi[i]);
    sprintf(str,"h3_var%d",i+1);
    h3[i]  = new TH1F(str,"",bn->nbin[i],bn->xlo[i],bn->xhi[i]);
    sprintf(str,"h4_var%d",i+1);
    h4[i]  = new TH1F(str,"",bn->nbin[i],bn->xlo[i],bn->xhi[i]);
  }

  TH1F* hh1 = new TH1F("hh1","",20,0,1); 
  TH1F* hh2 = new TH1F("hh2","",20,0,1);

  TH1F* hsig = new TH1F("sig","",40,0,2000);
  TH1F* hbkg1 = new TH1F("hbkg1","",40,0,2000);
  TH1F* hbkg2 = new TH1F("hbkg2","",40,0,2000);
  TH1F* hbkg3 = new TH1F("hbkg3","",40,0,2000);

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

      if(Pass(ich)) {

	n_eve[ich][1] += wt;

	if(ich==0) hsig->Fill(mV3/GeV,wt);
	else if(ich==1) hbkg1->Fill(mV3/GeV,wt);
	else if(ich==2||ich==3) hbkg1->Fill(mV3/GeV,wt);
	else if(ich==4||ich==5) hbkg2->Fill(mV3/GeV,wt);

	for(Int_t i=0; i<bn->nvar; i++) {
	  if(ich==0) h1[i]->Fill(bn->getVal(i),wt);
	  else if(ich==1) h2[i]->Fill(bn->getVal(i),wt);
	  else if(ich==4||ich==5||ich==2||ich==3) h3[i]->Fill(bn->getVal(i),wt);
	  // else if(ich==4||ich==5) h4[i]->Fill(bn->getVal(i),wt);
	}

	if(ich==0) {
	  if(isFatJet_v3) hh1->Fill(tau21_v3,wt);
	  if(isFatJet_v3) n_eve[0][2] += wt;
	}
	else {
	  if(isFatJet_v3) hh2->Fill(tau21_v3,wt);
	  if(isFatJet_v3) n_eve[1][2] += wt;
	}
      }
    }
  }
  printf("Sig: %.5f %.3f\n",n_eve[0][0],n_eve[0][1]);
  printf("Bkg1: %.5f %.3f\n",n_eve[1][0],n_eve[1][1]);
  printf("Bkg2: %.5f %.3f\n",n_eve[2][0],n_eve[2][1]);
  printf("Bkg3: %.5f %.3f\n",n_eve[3][0],n_eve[3][1]);
  printf("Bkg4: %.5f %.3f\n",n_eve[4][0],n_eve[4][1]);
  printf("Bkg4: %.5f %.3f\n",n_eve[5][0],n_eve[5][1]);
  printf("HasFatJet: %.1f(sig) %.1f(bkg)\n",n_eve[0][2],n_eve[1][2]);
  printf("significance: %0.3f\n", significance(n_eve[1][1]+n_eve[2][1]+n_eve[3][1]+n_eve[4][1]+n_eve[5][1], n_eve[0][1], 0) );
  
  gStyle->SetLegendFont(62);
  TLegend* lg = new TLegend(0.6,0.65,0.88,0.87,"");
  if(figSigSF==1) sprintf(str,"Signal");
  else sprintf(str,"Sig.#times%.1f",figSigSF);
  lg->AddEntry(h1[0],"Signal","F");
  lg->AddEntry(h2[0],"WZ+jets","F");
  lg->AddEntry(h3[0],"Other","F");
  // lg->AddEntry(h3[0],"t#bar{t}V","F");
  // lg->AddEntry(h4[0],"VVV","F");
  lg->SetBorderSize(0);
  lg->SetMargin(0.25);
  lg->SetFillColor(kWhite);

  TCanvas* c01 = new TCanvas("c01","c01",600,600);
  hh1->SetLineColor(kRed);
  hh2->SetLineColor(kBlue);
  hh1->SetMarkerColor(kRed);
  hh2->SetMarkerColor(kBlue);
  hh1->SetMarkerSize(0.8);
  hh2->SetMarkerSize(0.8);
  hh1->Scale(hh2->Integral()/hh1->Integral());
  SetMax(hh1,hh2,1.1);
  hh1->Draw();
  hh2->Draw("same");
  //////////////////////////////////////////////
   TStyle *atlasStyle= new TStyle("ATLAS","Atlas style");

  // use plain black on white colors
  Int_t icol=0;
  atlasStyle->SetFrameBorderMode(icol);
  atlasStyle->SetCanvasBorderMode(icol);
  atlasStyle->SetPadBorderMode(icol);
  atlasStyle->SetPadColor(icol);
  atlasStyle->SetCanvasColor(icol);
  atlasStyle->SetStatColor(icol);
  //atlasStyle->SetFillColor(icol);

  // set the paper & margin sizes
  atlasStyle->SetPaperSize(20,26);
  atlasStyle->SetPadTopMargin(0.07);
  atlasStyle->SetPadRightMargin(0.07);
  atlasStyle->SetPadBottomMargin(0.13);
  atlasStyle->SetPadLeftMargin(0.13);

  // use large fonts
  //Int_t font=72;
  Int_t font=42;
  Double_t tsize=0.05;
  atlasStyle->SetTextFont(font);

  atlasStyle->SetTextSize(tsize);
  atlasStyle->SetLabelFont(font,"x");
  atlasStyle->SetTitleFont(font,"x");
  atlasStyle->SetLabelFont(font,"y");
  atlasStyle->SetTitleFont(font,"y");
  atlasStyle->SetLabelFont(font,"z");
  atlasStyle->SetTitleFont(font,"z");

  atlasStyle->SetLabelSize(tsize,"x");
  atlasStyle->SetTitleSize(tsize,"x");
  atlasStyle->SetLabelSize(tsize,"y");
  atlasStyle->SetTitleSize(tsize,"y");
  atlasStyle->SetLabelSize(tsize,"z");
  atlasStyle->SetTitleSize(tsize,"z");

  // Xin:
  atlasStyle->SetTitleOffset(1.1,"X");
  atlasStyle->SetTitleOffset(1.2,"Y");

  //use bold lines and markers
  atlasStyle->SetMarkerStyle(20);
  atlasStyle->SetMarkerSize(1.2);
  atlasStyle->SetHistLineWidth(2.);
  atlasStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  //get rid of X error bars and y error bar caps
  //atlasStyle->SetErrorX(0.001);

  //do not display any of the standard histogram decorations
  atlasStyle->SetOptTitle(0);
  //atlasStyle->SetOptStat(1111);
  atlasStyle->SetOptStat(0);
  //atlasStyle->SetOptFit(1111);
  atlasStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  atlasStyle->SetPadTickX(1);
  atlasStyle->SetPadTickY(1);

  //gROOT->SetStyle("Plain");
  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();
  //////////////////////////////////////////////
  THStack* hsk[50];
  TCanvas* cv[50];
  for(Int_t i=0; i<bn->nvar; i++) {
    h1[i]->Scale(figSigSF);
    h1[i]->SetLineColor(kRed);
    h1[i]->SetFillColor(kRed);
    h2[i]->SetLineColor(kBlue);
    h2[i]->SetFillColor(kBlue);
    h3[i]->SetLineColor(kGreen-3);
    h3[i]->SetFillColor(kGreen-3);
    // h4[i]->SetLineColor(kBlue-5);
    // h4[i]->SetFillColor(kBlue-5);

    sprintf(str,"hsk_var%d",i+1);
    hsk[i] = new THStack(str,"");
    if(bn->titleX[i]=="m_{V(jj)} [GeV]" || bn->titleX[i]=="#Delta#phi_{LL}") {
      hsk[i]->Add(h4[i]);
      hsk[i]->Add(h3[i]);
      hsk[i]->Add(h2[i]);
      hsk[i]->Add(h1[i]);
    }
    else {
      hsk[i]->Add(h1[i]);
      hsk[i]->Add(h4[i]);
      hsk[i]->Add(h3[i]);
      hsk[i]->Add(h2[i]);
    }
    
    sprintf(str,"c%d",i+1);
    cv[i] = new TCanvas(str,str,912,600);
    hsk[i]->Draw("hist");
    hsk[i]->GetXaxis()->SetTitle(bn->titleX[i].c_str());
    sprintf(str,"Events/( %4.2f )",h1[i]->GetBinWidth(1));
    hsk[i]->GetYaxis()->SetTitle(str);
    lg->Draw("same");

    // sprintf(str,"c%d.eps",i+1);
    // cv[i]->SaveAs(str);
  }
}
