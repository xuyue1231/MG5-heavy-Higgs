/*
Simple macro showing how to access branches from the delphes output root file,
loop over events, and plot simple quantities such as the jet pt and the di-electron invariant
mass.

root -l examples/Example1.C'("delphes_output.root")'
*/

#include "TROOT.h"
#include "TH2F.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

//const Double_t M_PI=3.1415926536;
#define GeV 1000
#define c_light 2.99792458e8

//------------------------------------------------------------------------------
Int_t           EventNumber;
Float_t         mc_event_weight;
Float_t         mc_event_ScalePDF;
Float_t         mc_event_AlphaQED;
Float_t         mc_event_AlphaQCD;
Int_t           mu_n;
Float_t         mu_E[10];   //[mu_n]
Float_t         mu_pt[10];   //[mu_n]
Float_t         mu_eta[10];   //[mu_n]
Float_t         mu_phi[10];   //[mu_n]
Int_t           mu_charge[10];   //[mu_n]
Int_t           el_n;
Float_t         el_E[10];   //[el_n]
Float_t         el_pt[10];   //[el_n]
Float_t         el_eta[10];   //[el_n]
Float_t         el_phi[10];   //[el_n]
Int_t           el_charge[10];   //[el_n]
Int_t           jet_n;
Float_t         jet_E[50];   //[jet_n]
Float_t         jet_m[50];   //[jet_n]
Float_t         jet_pt[50];   //[jet_n]
Float_t         jet_eta[50];   //[jet_n]
Float_t         jet_phi[50];   //[jet_n]
Float_t         jet_tau1[50];   //[jet_n]
Float_t         jet_tau2[50];   //[jet_n]
Float_t         jet_tau3[50];   //[jet_n]
Float_t         jet_tau4[50];   //[jet_n]
Bool_t          jet_isBtagged[50];   //[jet_n]
Int_t           fatjet_n;
Float_t         fatjet_E[50];   //[fatjet_n]
Float_t         fatjet_m[50];   //[fatjet_n]
Float_t         fatjet_pt[50];   //[fatjet_n]
Float_t         fatjet_eta[50];   //[fatjet_n]
Float_t         fatjet_phi[50];   //[fatjet_n]
Float_t         fatjet_tau1[50];   //[fatjet_n]
Float_t         fatjet_tau2[50];   //[fatjet_n]
Float_t         fatjet_tau3[50];   //[fatjet_n]
Float_t         fatjet_tau4[50];   //[fatjet_n]
Int_t           fatjet_sj_n[50];   //[fatjet_n]
Float_t         fatjet_m_sj1[50];   //[fatjet_n]
Float_t         fatjet_pt_sj1[50];   //[fatjet_n]
Float_t         fatjet_eta_sj1[50];   //[fatjet_n]
Float_t         fatjet_phi_sj1[50];   //[fatjet_n]
Float_t         fatjet_m_sj2[50];   //[fatjet_n]
Float_t         fatjet_pt_sj2[50];   //[fatjet_n]
Float_t         fatjet_eta_sj2[50];   //[fatjet_n]
Float_t         fatjet_phi_sj2[50];   //[fatjet_n]
Float_t         MET_ex;
Float_t         MET_ey;
Float_t         MET_et;
Float_t         MET_phi;
Float_t         MET_sumet;
Float_t         MET_Truth_NonInt_ex;
Float_t         MET_Truth_NonInt_ey;
Float_t         MET_Truth_NonInt_et;
Int_t           jet_truth_n;
Float_t         jet_truth_pt[50];   //[jet_truth_n]
Float_t         jet_truth_eta[50];   //[jet_truth_n]
Float_t         jet_truth_phi[50];   //[jet_truth_n]
Float_t         jet_truth_m[50];   //[jet_truth_n]
Int_t           fatjet_truth_n;
Float_t         fatjet_truth_pt[50];   //[fatjet_truth_n]
Float_t         fatjet_truth_eta[50];   //[fatjet_truth_n]
Float_t         fatjet_truth_phi[50];   //[fatjet_truth_n]
Float_t         fatjet_truth_m[50];   //[fatjet_truth_n]
Int_t           truthV_n;
Int_t           truthV_pdgId[10];   //[truthV_n]
Float_t         truthV_pt[10];   //[truthV_n]
Float_t         truthV_m[10];   //[truthV_n]
Float_t         truthV_eta[10];   //[truthV_n]
Float_t         truthV_phi[10];   //[truthV_n]
Int_t           truthV_charge[10];   //[truthV_n]
Int_t           truthV_dau1_pdgId[10];   //[truthV_n]
Float_t         truthV_dau1_pt[10];   //[truthV_n]
Float_t         truthV_dau1_m[10];   //[truthV_n]
Float_t         truthV_dau1_eta[10];   //[truthV_n]
Float_t         truthV_dau1_phi[10];   //[truthV_n]
Int_t           truthV_dau1_charge[10];   //[truthV_n]
Int_t           truthV_dau2_pdgId[10];   //[truthV_n]
Float_t         truthV_dau2_pt[10];   //[truthV_n]
Float_t         truthV_dau2_m[10];   //[truthV_n]
Float_t         truthV_dau2_eta[10];   //[truthV_n]
Float_t         truthV_dau2_phi[10];   //[truthV_n]
Int_t           truthV_dau2_charge[10];   //[truthV_n]

void SetBranch(TTree* tree) {
  tree->Branch("EventNumber", &EventNumber,"EventNumber/I");
  tree->Branch("mc_event_weight", &mc_event_weight,"mc_event_weight/F");
  //tree->Branch("mc_event_ScalePDF", &mc_event_ScalePDF,"mc_event_ScalePDF/F");
  //tree->Branch("mc_event_AlphaQED", &mc_event_AlphaQED,"mc_event_AlphaQED/F");
  //tree->Branch("mc_event_AlphaQCD", &mc_event_AlphaQCD,"mc_event_AlphaQCD/F");
  tree->Branch("mu_n", &mu_n,"mu_n/I");
  tree->Branch("mu_E", mu_E,"mu_E[mu_n]/F");
  tree->Branch("mu_pt", mu_pt,"mu_pt[mu_n]/F");
  tree->Branch("mu_eta", mu_eta,"mu_eta[mu_n]/F");
  tree->Branch("mu_phi", mu_phi,"mu_phi[mu_n]/F");
  tree->Branch("mu_charge", mu_charge,"mu_charge[mu_n]/I");
  tree->Branch("el_n", &el_n,"el_n/I");
  tree->Branch("el_E", el_E,"el_E[el_n]/F");
  tree->Branch("el_pt", el_pt,"el_pt[el_n]/F");
  tree->Branch("el_eta", el_eta,"el_eta[el_n]/F");
  tree->Branch("el_phi", el_phi,"el_phi[el_n]/F");
  tree->Branch("el_charge", el_charge,"el_charge[el_n]/I");
  tree->Branch("jet_n", &jet_n,"jet_n/I");
  tree->Branch("jet_E", jet_E,"jet_E[jet_n]/F");
  tree->Branch("jet_m", jet_m,"jet_m[jet_n]/F");
  tree->Branch("jet_pt", jet_pt,"jet_pt[jet_n]/F");
  tree->Branch("jet_eta", jet_eta,"jet_eta[jet_n]/F");
  tree->Branch("jet_phi", jet_phi,"jet_phi[jet_n]/F");
  tree->Branch("jet_tau1", jet_tau1,"jet_tau1[jet_n]/F");
  tree->Branch("jet_tau2", jet_tau2,"jet_tau2[jet_n]/F");
  tree->Branch("jet_tau3", jet_tau3,"jet_tau3[jet_n]/F");
  tree->Branch("jet_tau4", jet_tau4,"jet_tau4[jet_n]/F");
  tree->Branch("jet_isBtagged", jet_isBtagged,"jet_isBtagged[jet_n]/O");
  tree->Branch("fatjet_n", &fatjet_n,"fatjet_n/I");
  tree->Branch("fatjet_E", fatjet_E,"fatjet_E[fatjet_n]/F");
  tree->Branch("fatjet_m", fatjet_m,"fatjet_m[fatjet_n]/F");
  tree->Branch("fatjet_pt", fatjet_pt,"fatjet_pt[fatjet_n]/F");
  tree->Branch("fatjet_eta", fatjet_eta,"fatjet_eta[fatjet_n]/F");
  tree->Branch("fatjet_phi", fatjet_phi,"fatjet_phi[fatjet_n]/F");
  tree->Branch("fatjet_tau1", fatjet_tau1,"fatjet_tau1[fatjet_n]/F");
  tree->Branch("fatjet_tau2", fatjet_tau2,"fatjet_tau2[fatjet_n]/F");
  tree->Branch("fatjet_tau3", fatjet_tau3,"fatjet_tau3[fatjet_n]/F");
  tree->Branch("fatjet_tau4", fatjet_tau4,"fatjet_tau4[fatjet_n]/F");
  tree->Branch("fatjet_sj_n", fatjet_sj_n,"fatjet_sj_n[fatjet_n]/I");
  tree->Branch("fatjet_m_sj1", fatjet_m_sj1,"fatjet_m_sj1[fatjet_n]/F");
  tree->Branch("fatjet_pt_sj1", fatjet_pt_sj1,"fatjet_pt_sj1[fatjet_n]/F");
  tree->Branch("fatjet_eta_sj1", fatjet_eta_sj1,"fatjet_eta_sj1[fatjet_n]/F");
  tree->Branch("fatjet_phi_sj1", fatjet_phi_sj1,"fatjet_phi_sj1[fatjet_n]/F");
  tree->Branch("fatjet_m_sj2", fatjet_m_sj2,"fatjet_m_sj2[fatjet_n]/F");
  tree->Branch("fatjet_pt_sj2", fatjet_pt_sj2,"fatjet_pt_sj2[fatjet_n]/F");
  tree->Branch("fatjet_eta_sj2", fatjet_eta_sj2,"fatjet_eta_sj2[fatjet_n]/F");
  tree->Branch("fatjet_phi_sj2", fatjet_phi_sj2,"fatjet_phi_sj2[fatjet_n]/F");
  tree->Branch("MET_ex", &MET_ex,"MET_ex/F");
  tree->Branch("MET_ey", &MET_ey,"MET_ey/F");
  tree->Branch("MET_et", &MET_et,"MET_et/F");
  tree->Branch("MET_phi", &MET_phi,"MET_phi/F");
  tree->Branch("MET_sumet", &MET_sumet,"MET_sumet/F");
  tree->Branch("MET_Truth_NonInt_ex", &MET_Truth_NonInt_ex,"MET_Truth_NonInt_ex/F");
  tree->Branch("MET_Truth_NonInt_ey", &MET_Truth_NonInt_ey,"MET_Truth_NonInt_ey/F");
  tree->Branch("MET_Truth_NonInt_et", &MET_Truth_NonInt_et,"MET_Truth_NonInt_et/F");
  tree->Branch("jet_truth_n", &jet_truth_n,"jet_truth_n/I");
  tree->Branch("jet_truth_pt", jet_truth_pt,"jet_truth_pt[jet_truth_n]/F");
  tree->Branch("jet_truth_eta", jet_truth_eta,"jet_truth_eta[jet_truth_n]/F");
  tree->Branch("jet_truth_phi", jet_truth_phi,"jet_truth_phi[jet_truth_n]/F");
  tree->Branch("jet_truth_m", jet_truth_m,"jet_truth_m[jet_truth_n]/F");
  tree->Branch("fatjet_truth_n", &fatjet_truth_n,"fatjet_truth_n/I");
  tree->Branch("fatjet_truth_pt", fatjet_truth_pt,"fatjet_truth_pt[fatjet_truth_n]/F");
  tree->Branch("fatjet_truth_eta", fatjet_truth_eta,"fatjet_truth_eta[fatjet_truth_n]/F");
  tree->Branch("fatjet_truth_phi", fatjet_truth_phi,"fatjet_truth_phi[fatjet_truth_n]/F");
  tree->Branch("fatjet_truth_m", fatjet_truth_m,"fatjet_truth_m[fatjet_truth_n]/F");
  tree->Branch("truthV_n", &truthV_n,"truthV_n/I");
  tree->Branch("truthV_pdgId", truthV_pdgId,"truthV_pdgId[truthV_n]/I");
  tree->Branch("truthV_pt", truthV_pt,"truthV_pt[truthV_n]/F");
  tree->Branch("truthV_m", truthV_m,"truthV_m[truthV_n]/F");
  tree->Branch("truthV_eta", truthV_eta,"truthV_eta[truthV_n]/F");
  tree->Branch("truthV_phi", truthV_phi,"truthV_phi[truthV_n]/F");
  tree->Branch("truthV_charge", truthV_charge,"truthV_charge[truthV_n]/I");
  tree->Branch("truthV_dau1_pdgId", truthV_dau1_pdgId,"truthV_dau1_pdgId[truthV_n]/I");
  tree->Branch("truthV_dau1_pt", truthV_dau1_pt,"truthV_dau1_pt[truthV_n]/F");
  tree->Branch("truthV_dau1_m", truthV_dau1_m,"truthV_dau1_m[truthV_n]/F");
  tree->Branch("truthV_dau1_eta", truthV_dau1_eta,"truthV_dau1_eta[truthV_n]/F");
  tree->Branch("truthV_dau1_phi", truthV_dau1_phi,"truthV_dau1_phi[truthV_n]/F");
  tree->Branch("truthV_dau1_charge", truthV_dau1_charge,"truthV_dau1_charge[truthV_n]/I");
  tree->Branch("truthV_dau2_pdgId", truthV_dau2_pdgId,"truthV_dau2_pdgId[truthV_n]/I");
  tree->Branch("truthV_dau2_pt", truthV_dau2_pt,"truthV_dau2_pt[truthV_n]/F");
  tree->Branch("truthV_dau2_m", truthV_dau2_m,"truthV_dau2_m[truthV_n]/F");
  tree->Branch("truthV_dau2_eta", truthV_dau2_eta,"truthV_dau2_eta[truthV_n]/F");
  tree->Branch("truthV_dau2_phi", truthV_dau2_phi,"truthV_dau2_phi[truthV_n]/F");
  tree->Branch("truthV_dau2_charge", truthV_dau2_charge,"truthV_dau2_charge[truthV_n]/I");
}

Double_t Max(Double_t a, Double_t b) {
  if(a>b) return a;
  else return b;
}

void SetMax(TH1* h1, TH1* h2, Double_t scale=1.0) {
  h1->SetMaximum(scale*Max(h1->GetMaximum(),h2->GetMaximum()));
  h2->SetMaximum(scale*Max(h1->GetMaximum(),h2->GetMaximum()));
}

Float_t DR(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2) {
  Float_t dphi = fabs(phi1-phi2);
  if(dphi>M_PI) dphi = 2*M_PI-dphi;
  Float_t deta = fabs(eta1-eta2);
  return sqrt(dphi*dphi+deta*deta);
}

Bool_t isNotSelf(TClonesArray *branchGenParticle, GenParticle *particle) {
  if(particle->D1>=0 && particle->D2>=0) {
    for(Int_t i=particle->D1; i<=particle->D2; i++) {
      GenParticle *part = (GenParticle*) branchGenParticle->At(i);
      if(part->PID==particle->PID) return false;
    }
  }
  return true;
}

Int_t GetMotherID(TClonesArray *branchGenParticle, GenParticle *particle, bool removeSelf=false) {
  if(particle->M1>=0 && particle->M2>=0) {
    for(Int_t i=particle->M1; i<=particle->M2; i++) {
      GenParticle *mon = (GenParticle*) branchGenParticle->At(i);
      if(mon->PID==particle->PID && removeSelf) return GetMotherID(branchGenParticle,mon,removeSelf);
    }
    GenParticle *part = (GenParticle*) branchGenParticle->At(particle->M1);
    return part->PID;
  }
  else if(particle->M1>=0) {
    GenParticle *part = (GenParticle*) branchGenParticle->At(particle->M1);
    if(part->PID==particle->PID && removeSelf) return GetMotherID(branchGenParticle,part,removeSelf);
    return part->PID;
  }
  return 0;
}

void Check(TClonesArray *branchGenParticle, GenParticle *particle) {
  GenParticle *self(0);
  if(particle->D1>=0 && particle->D2>=0) {
    for(Int_t i=particle->D1; i<=particle->D2; i++) {
      GenParticle *part = (GenParticle*) branchGenParticle->At(i);
      printf(" %d (%d %f/%f/%f)",part->PID,part->Status,part->E,part->Eta,part->Phi);
      if(particle->PID==part->PID) self = part;
    }
    printf("\n");
  }
  if(self && self->Status!=1) Check(branchGenParticle,self);
  return;
}

GenParticle* GetFirstSelf(TClonesArray *branchGenParticle, GenParticle *particle) {
  if(particle->M1>=0) {
    for(Int_t i=particle->M1; i<=TMath::Max(particle->M1,particle->M2); i++) {
      GenParticle *part = (GenParticle*) branchGenParticle->At(i);
      if(particle->PID==part->PID) {
	Bool_t hasGam = false;
	for(Int_t j=part->D1; j<=part->D2; j++) {
	  GenParticle *dau = (GenParticle*) branchGenParticle->At(j);
	  if(dau->PID==22) hasGam = true;
	  if(dau->PID!=part->PID && dau->PID!=22) printf("Error: %d has daughter %d\n",part->PID,dau->PID);
	}
	if(hasGam && part->E>particle->E) return GetFirstSelf(branchGenParticle,part);
      }
    }
  }
  return particle;
}

GenParticle* GetFinalSelf(TClonesArray *branchGenParticle, GenParticle *particle) {
  if(particle->Status==1) return particle;
  if(particle->D1>=0 && particle->D2>=0) {
    for(Int_t i=particle->D1; i<=particle->D2; i++) {
      GenParticle *part = (GenParticle*) branchGenParticle->At(i);
      if(particle->PID==part->PID) {
	return GetFinalSelf(branchGenParticle,part);
      }
    }
  }
  return particle;
}

GenParticle* GetFinalParticle(TClonesArray *branchGenParticle, GenParticle *particle) {
  if(particle->Status==1) return particle;
  if(particle->D1>=0 && particle->D2>=0) {
    for(Int_t i=particle->D1; i<=particle->D2; i++) {
      GenParticle *part = (GenParticle*) branchGenParticle->At(i);
      if(particle->PID==part->PID) {
	if(part->Status==1) return part;
	else return GetFinalParticle(branchGenParticle,part);
      }
    }
  }
  return particle;
}

GenParticle* GetFirstParticle(TClonesArray *branchGenParticle, GenParticle *particle, Int_t aID) {
  if(abs(particle->PID)==aID) return particle;
  if(particle->D1>=0 && particle->D2>=0) {
    for(Int_t i=particle->D1; i<=particle->D2; i++) {
      GenParticle *part = (GenParticle*) branchGenParticle->At(i);
      if(abs(part->PID)==aID) return part;
      else if(part->Status!=1) {
	GenParticle* tmp = GetFirstParticle(branchGenParticle,part,aID);
	if(tmp) return tmp;
      }
    }
  }
  return 0;
}

void reduce(const char* input="../Events/run_01/tag_1_delphes_events.root", const char* out="ntuple/sig.root")
{
  //gSystem->Load("libDelphes");
  char str[200];

  // split by ','
  std::string argStr = input;
  std::vector<std::string> fileList;
  for(size_t i=0,n; i <= argStr.length(); i=n+1) {
    n = argStr.find_first_of(',',i);
    if (n == std::string::npos) n = argStr.length();
    std::string tmp = argStr.substr(i,n-i);
    fileList.push_back(tmp);
  }
  
  // Create chain of root trees
  TChain chain("Delphes");
  for(std::vector<std::string>::size_type i=0; i<fileList.size(); i++) {
    chain.Add(fileList[i].c_str());
  }

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Int_t numberOfEntries = (Int_t)treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchFatJet = treeReader->UseBranch("FatJet");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchGenFatJet = treeReader->UseBranch("GenFatJet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
  TClonesArray *branchGenParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchScalarHT = treeReader->UseBranch("ScalarHT");
  TClonesArray *branchTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchNeutral = treeReader->UseBranch("EFlowPhoton");
  TClonesArray *branchNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");

  Int_t eve1 = 0;
  Int_t eve2 = 0;
  Int_t eve3 = 0;

  Float_t n_eve[] = {0,0,0};

  TFile* fout = new TFile(out,"RECREATE"); //output file
  TTree* tree = new TTree("tau",str); //output tree
  SetBranch(tree);

  // Loop over all events
  printf(" %d entries to be processed\n",numberOfEntries);
  for(Int_t entry = 0; entry < numberOfEntries; entry++)
  //for(Int_t entry = 0; entry < 10000; entry++)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    if((entry+1)%10000==0) printf(" -> %d entries processed...\n",entry+1);

    HepMCEvent *m_Event = (HepMCEvent*) branchEvent->At(0);
    EventNumber = m_Event->Number;
    mc_event_weight = m_Event->Weight;
    mc_event_ScalePDF = m_Event->ScalePDF;
    mc_event_AlphaQED = m_Event->AlphaQED;
    mc_event_AlphaQCD = m_Event->AlphaQCD;

    eve1++;
    n_eve[0] += mc_event_weight;

    el_n = 0;
    for(Int_t i=0; i<branchElectron->GetEntries(); i++) {
      Electron *ele = (Electron*) branchElectron->At(i);
      if(ele->PT>15 && fabs(ele->Eta)<2.5) {
	el_E[el_n] = ele->PT*cosh(ele->Eta)*GeV;
	el_pt[el_n] = ele->PT*GeV;
	el_eta[el_n] = ele->Eta;
	el_phi[el_n] = ele->Phi;
	el_charge[el_n] = ele->Charge;
	el_n++;
      }
    }

    mu_n = 0;
    for(Int_t i=0; i<branchMuon->GetEntries(); i++) {
      Muon *muo = (Muon*) branchMuon->At(i);
      if(muo->PT>10 && fabs(muo->Eta)<2.5) {
	mu_E[mu_n] = muo->PT*cosh(muo->Eta)*GeV;
	mu_pt[mu_n] = muo->PT*GeV;
	mu_eta[mu_n] = muo->Eta;
	mu_phi[mu_n] = muo->Phi;
	mu_charge[mu_n] = muo->Charge;
	mu_n++;
      }
    }

    truthV_n = 0;

    for(Int_t i=0; i<branchGenParticle->GetEntries(); i++) {
      GenParticle *part = (GenParticle*) branchGenParticle->At(i);
      if(!part) continue;

      if((part->PID==23 || abs(part->PID)==24) && isNotSelf(branchGenParticle,part)) {
        truthV_pdgId[truthV_n] = part->PID;
        truthV_pt[truthV_n] = part->PT*GeV;
        truthV_m[truthV_n] = part->Mass*GeV;
        truthV_eta[truthV_n] = part->Eta;
        truthV_phi[truthV_n] = part->Phi;
        truthV_charge[truthV_n] = part->Charge;
	
	truthV_dau1_pdgId[truthV_n] = 0;
        truthV_dau1_pt[truthV_n] = 0;
        truthV_dau1_m[truthV_n] = 0;
        truthV_dau1_eta[truthV_n] = 0;
        truthV_dau1_phi[truthV_n] = 0;
        truthV_dau1_charge[truthV_n] = 0;
	truthV_dau2_pdgId[truthV_n] = 0;
        truthV_dau2_pt[truthV_n] = 0;
        truthV_dau2_m[truthV_n] = 0;
        truthV_dau2_eta[truthV_n] = 0;
        truthV_dau2_phi[truthV_n] = 0;
        truthV_dau2_charge[truthV_n] = 0;

	for(Int_t j=part->D1; j<=part->D2; j++) {
	  GenParticle *dau = (GenParticle*) branchGenParticle->At(j);
	  if(truthV_dau1_pdgId[truthV_n]==0) {
	    truthV_dau1_pdgId[truthV_n] = dau->PID;
	    truthV_dau1_pt[truthV_n] = dau->PT*GeV;
	    truthV_dau1_m[truthV_n] = dau->Mass*GeV;
	    truthV_dau1_eta[truthV_n] = dau->Eta;
	    truthV_dau1_phi[truthV_n] = dau->Phi;
	    truthV_dau1_charge[truthV_n] = dau->Charge;
	  }
	  else if(truthV_dau2_pdgId[truthV_n]==0) {
	    truthV_dau2_pdgId[truthV_n] = dau->PID;
	    truthV_dau2_pt[truthV_n] = dau->PT*GeV;
	    truthV_dau2_m[truthV_n] = dau->Mass*GeV;
	    truthV_dau2_eta[truthV_n] = dau->Eta;
	    truthV_dau2_phi[truthV_n] = dau->Phi;
	    truthV_dau2_charge[truthV_n] = dau->Charge;
	  }
	}
	truthV_n++;
      }
    }

    jet_n = 0;
    for(Int_t i=0; i<branchJet->GetEntries(); i++) {
      Jet *jet = (Jet*) branchJet->At(i);
      
      Int_t iEl = -1;
      for(Int_t j=0; j<el_n; j++) {
	if(DR(el_eta[j],el_phi[j],jet->Eta,jet->Phi)<0.4) {
	  iEl = j;
	}
      }
      if(iEl>=0) continue;

      Int_t iMu = -1;
      for(Int_t j=0; j<mu_n; j++) {
	if(DR(mu_eta[j],mu_phi[j],jet->Eta,jet->Phi)<0.4) {
	  iMu = j;
	}
      }
      if(iMu>=0) continue;

      if(jet->PT>30 && fabs(jet->Eta)<4.0) {
	jet_E[jet_n] = sqrt(pow(jet->PT*cosh(jet->Eta),2)+pow(jet->Mass,2))*GeV;
	jet_m[jet_n] = jet->Mass*GeV;
	jet_pt[jet_n] = jet->PT*GeV;
	jet_eta[jet_n] = jet->Eta;
	jet_phi[jet_n] = jet->Phi;
	jet_tau1[jet_n] = jet->Tau[0];
        jet_tau2[jet_n] = jet->Tau[1];
        jet_tau3[jet_n] = jet->Tau[2];
        jet_tau4[jet_n] = jet->Tau[3];
	jet_isBtagged[jet_n] = fabs(jet->Eta)<2.5 ? (jet->BTag % 2)==1 : false;
	jet_n++;
      }
    }

    fatjet_n = 0;
    for(Int_t i=0; i<branchFatJet->GetEntries(); i++) {
      Jet *jet = (Jet*) branchFatJet->At(i);
      
      Int_t iEl = -1;
      for(Int_t j=0; j<el_n; j++) {
	if(DR(el_eta[j],el_phi[j],jet->Eta,jet->Phi)<1.0) {
	  iEl = j;
	}
      }
      if(iEl>=0) continue;

      Int_t iMu = -1;
      for(Int_t j=0; j<mu_n; j++) {
	if(DR(mu_eta[j],mu_phi[j],jet->Eta,jet->Phi)<1.0) {
	  iMu = j;
	}
      }
      if(iMu>=0) continue;

      if(jet->PT>40 && fabs(jet->Eta)<4.0) {
	fatjet_E[fatjet_n] = sqrt(pow(jet->PT*cosh(jet->Eta),2)+pow(jet->Mass,2))*GeV;
	fatjet_m[fatjet_n] = jet->Mass*GeV;
	fatjet_pt[fatjet_n] = jet->PT*GeV;
	fatjet_eta[fatjet_n] = jet->Eta;
	fatjet_phi[fatjet_n] = jet->Phi;
	fatjet_tau1[fatjet_n] = jet->Tau[0];
        fatjet_tau2[fatjet_n] = jet->Tau[1];
        fatjet_tau3[fatjet_n] = jet->Tau[2];
        fatjet_tau4[fatjet_n] = jet->Tau[3];
	if(jet->SoftDroppedP4[0].Pt() && jet->SoftDroppedP4[1].Pt()>0) fatjet_sj_n[fatjet_n] = 2;
	else if(jet->SoftDroppedP4[0].Pt() || jet->SoftDroppedP4[1].Pt()>0) fatjet_sj_n[fatjet_n] = 1;
	else fatjet_sj_n[fatjet_n] = 0;
        fatjet_m_sj1[fatjet_n] = jet->SoftDroppedP4[0].M()*GeV;
        fatjet_pt_sj1[fatjet_n] = jet->SoftDroppedP4[0].Pt()*GeV;
        fatjet_eta_sj1[fatjet_n] = jet->SoftDroppedP4[0].Eta();
        fatjet_phi_sj1[fatjet_n] = jet->SoftDroppedP4[0].Phi();
        fatjet_m_sj2[fatjet_n] = jet->SoftDroppedP4[1].M()*GeV;
        fatjet_pt_sj2[fatjet_n] = jet->SoftDroppedP4[1].Pt()*GeV;
        fatjet_eta_sj2[fatjet_n] = jet->SoftDroppedP4[1].Eta();
        fatjet_phi_sj2[fatjet_n] = jet->SoftDroppedP4[1].Phi();
	fatjet_n++;
      }
    }

    jet_truth_n = 0;
    for(Int_t i=0; i<branchGenJet->GetEntries(); i++) {
      Jet *jet = (Jet*) branchGenJet->At(i);
      jet_truth_pt[jet_truth_n] = jet->PT*GeV;
      jet_truth_eta[jet_truth_n] = jet->Eta;
      jet_truth_phi[jet_truth_n] = jet->Phi;
      jet_truth_m[jet_truth_n] = jet->Mass*GeV;
      jet_truth_n++;
    }

    fatjet_truth_n = 0;
    for(Int_t i=0; i<branchGenFatJet->GetEntries(); i++) {
      Jet *jet = (Jet*) branchGenFatJet->At(i);
      fatjet_truth_pt[fatjet_truth_n] = jet->PT*GeV;
      fatjet_truth_eta[fatjet_truth_n] = jet->Eta;
      fatjet_truth_phi[fatjet_truth_n] = jet->Phi;
      fatjet_truth_m[fatjet_truth_n] = jet->Mass*GeV;
      fatjet_truth_n++;
    }

    MissingET *m_MET = (MissingET*) branchMissingET->At(0);
    MET_ex = m_MET->MET*GeV*cos(m_MET->Phi);
    MET_ey = m_MET->MET*GeV*sin(m_MET->Phi);
    MET_et = m_MET->MET*GeV;
    MET_phi = m_MET->Phi;

    ScalarHT *m_HT = (ScalarHT*) branchScalarHT->At(0);
    MET_sumet = m_HT->HT*GeV;

    MET_Truth_NonInt_ex = 0;
    MET_Truth_NonInt_ey = 0;
    for(Int_t i=0; i<branchGenParticle->GetEntries(); i++) {
      GenParticle *part = (GenParticle*) branchGenParticle->At(i);
      if(!part) continue;
      if(abs(part->PID)==12 || abs(part->PID)==14 || abs(part->PID)==16) {
	MET_Truth_NonInt_ex += part->Px*GeV;
	MET_Truth_NonInt_ey += part->Py*GeV;
      }
    }
    MET_Truth_NonInt_et = sqrt(pow(MET_Truth_NonInt_ex,2)+pow(MET_Truth_NonInt_ey,2));

    //if((el_n==2 && mu_n==0) || (el_n==0 && mu_n==2)) {
    if(true) {
      eve2++;
      n_eve[1] += mc_event_weight;

      tree->Fill();
    }
  }
  printf(" Tot %d, Sel %d\n",eve1,eve2);
  printf(" Tot %g, Sel %g\n",n_eve[0],n_eve[1]);

  fout->cd();
  tree->AutoSave();
}
