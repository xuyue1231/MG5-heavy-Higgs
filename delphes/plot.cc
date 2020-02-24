#include "TROOT.h"
#include "TStyle.h"
#include "TH2F.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TMinuit.h"
#include "TFile.h"
#include "TVector2.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "reader.h"

#define GeV 1000

void plot(const char* infile="qqz_a1rho_I1_I2_red.root") {
  TChain chain("tau");
  chain.Add(infile);

  Init(&chain);
  
  TH1F* h3a = new TH1F("h3a","",100,0,200);
  TH1F* h3b = new TH1F("h3b","",100,0,200);

  Int_t numberOfEntries = (Int_t)chain.GetEntries();
  // Loop over all events
  printf(" %d entries to be processed\n",numberOfEntries);
  for(Int_t entry = 0; entry < numberOfEntries; entry++)
  {
    // Load selected branches with data from specified event
    chain.GetEntry(entry);

    TLorentzVector tau1;
    TLorentzVector tau2;
    for(Int_t i=0; i<truthTau_n; i++) {
      if(abs(truthTau_vis_pdgId[i])==0 && abs(truthTau_vis_pdgIdMoth[i])<100) {
	TLorentzVector tmp1;
	TLorentzVector tmp2;
	tmp1.SetPtEtaPhiM(truthTau_vis_pt[i],truthTau_vis_eta[i],truthTau_vis_phi[i],truthTau_vis_m[i]);
	tmp2.SetPtEtaPhiM(truthTau_inv_pt[i],truthTau_inv_eta[i],truthTau_inv_phi[i],truthTau_inv_m[i]);
	if(truthTau_vis_charge[i]<0) tau1 = tmp1+tmp2;
	else                         tau2 = tmp1+tmp2;
      }
    }

    if(mc_event_weight>0) h3a->Fill((tau1+tau2).M()/GeV,mc_event_weight);
    else h3b->Fill((tau1+tau2).M()/GeV,mc_event_weight);
  }
  
  TCanvas* c3a = new TCanvas("c3a","c3a",600,600);
  h3a->SetXTitle("m_{#tau#tau} [GeV]");
  //h3a->SetMinimum(0);
  h3a->Draw();
 
  TCanvas* c3b = new TCanvas("c3b","c3b",600,600);
  h3b->SetXTitle("m_{#tau#tau} [GeV]");
  //h3b->SetMinimum(0);
  h3b->Draw();

  printf("%f %f\n",h3a->Integral(h3a->GetXaxis()->FindBin(110),h3a->GetXaxis()->FindBin(123.9)),h3a->Integral(1,h3a->GetNbinsX()));
  printf("%f %f\n",h3b->Integral(h3b->GetXaxis()->FindBin(110),h3b->GetXaxis()->FindBin(120)),h3b->Integral(1,h3b->GetNbinsX()));
}
