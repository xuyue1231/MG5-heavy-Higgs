   // Declaration of leaf types
   Int_t           EventNumber;
   Float_t         mc_event_weight;
   Int_t           mu_n;
   Float_t         mu_E[4];   //[mu_n]
   Float_t         mu_pt[4];   //[mu_n]
   Float_t         mu_eta[4];   //[mu_n]
   Float_t         mu_phi[4];   //[mu_n]
   Int_t           mu_charge[4];   //[mu_n]
   Int_t           el_n;
   Float_t         el_E[4];   //[el_n]
   Float_t         el_pt[4];   //[el_n]
   Float_t         el_eta[4];   //[el_n]
   Float_t         el_phi[4];   //[el_n]
   Int_t           el_charge[4];   //[el_n]
   Int_t           jet_n;
   Float_t         jet_E[10];   //[jet_n]
   Float_t         jet_m[10];   //[jet_n]
   Float_t         jet_pt[10];   //[jet_n]
   Float_t         jet_eta[10];   //[jet_n]
   Float_t         jet_phi[10];   //[jet_n]
   Float_t         jet_tau1[10];  
   Float_t         jet_tau2[10];   
   Float_t         jet_tau3[10];   
   Float_t         jet_tau4[10];  
   Bool_t          jet_isBtagged[11]; 
   Int_t           fatjet_n;
   Float_t         fatjet_E[8];   //[fatjet_n]
   Float_t         fatjet_m[8];   //[fatjet_n]
   Float_t         fatjet_pt[8];   //[fatjet_n]
   Float_t         fatjet_eta[8];   //[fatjet_n]
   Float_t         fatjet_phi[8];   //[fatjet_n]
   Float_t         fatjet_tau1[8];   //[fatjet_n]
   Float_t         fatjet_tau2[8];   //[fatjet_n]
   Float_t         fatjet_tau3[8];   //[fatjet_n]
   Float_t         fatjet_tau4[8];   //[fatjet_n]
   Int_t           fatjet_sj_n[8];   //[fatjet_n]
   Float_t         fatjet_m_sj1[8];   //[fatjet_n]
   Float_t         fatjet_pt_sj1[8];   //[fatjet_n]
   Float_t         fatjet_eta_sj1[8];   //[fatjet_n]
   Float_t         fatjet_phi_sj1[8];   //[fatjet_n]
   Float_t         fatjet_m_sj2[8];   //[fatjet_n]
   Float_t         fatjet_pt_sj2[8];   //[fatjet_n]
   Float_t         fatjet_eta_sj2[8];   //[fatjet_n]
   Float_t         fatjet_phi_sj2[8];   //[fatjet_n]
   Float_t         MET_ex;
   Float_t         MET_ey;
   Float_t         MET_et;
   Float_t         MET_phi;
   Float_t         MET_sumet;
   Float_t         MET_Truth_NonInt_ex;
   Float_t         MET_Truth_NonInt_ey;
   Float_t         MET_Truth_NonInt_et;
   Int_t           jet_truth_n;
   Float_t         jet_truth_pt[15];   //[jet_truth_n]
   Float_t         jet_truth_eta[15];   //[jet_truth_n]
   Float_t         jet_truth_phi[15];   //[jet_truth_n]
   Float_t         jet_truth_m[15];   //[jet_truth_n]
   Int_t           fatjet_truth_n;
   Float_t         fatjet_truth_pt[10];   //[fatjet_truth_n]
   Float_t         fatjet_truth_eta[10];   //[fatjet_truth_n]
   Float_t         fatjet_truth_phi[10];   //[fatjet_truth_n]
   Float_t         fatjet_truth_m[10];   //[fatjet_truth_n]
   Int_t           truthV_n;
   Int_t           truthV_pdgId[3];   //[truthV_n]
   Float_t         truthV_pt[3];   //[truthV_n]
   Float_t         truthV_m[3];   //[truthV_n]
   Float_t         truthV_eta[3];   //[truthV_n]
   Float_t         truthV_phi[3];   //[truthV_n]
   Int_t           truthV_charge[3];   //[truthV_n]
   Int_t           truthV_dau1_pdgId[3];   //[truthV_n]
   Float_t         truthV_dau1_pt[3];   //[truthV_n]
   Float_t         truthV_dau1_m[3];   //[truthV_n]
   Float_t         truthV_dau1_eta[3];   //[truthV_n]
   Float_t         truthV_dau1_phi[3];   //[truthV_n]
   Int_t           truthV_dau1_charge[3];   //[truthV_n]
   Int_t           truthV_dau2_pdgId[3];   //[truthV_n]
   Float_t         truthV_dau2_pt[3];   //[truthV_n]
   Float_t         truthV_dau2_m[3];   //[truthV_n]
   Float_t         truthV_dau2_eta[3];   //[truthV_n]
   Float_t         truthV_dau2_phi[3];   //[truthV_n]
   Int_t           truthV_dau2_charge[3];   //[truthV_n]

   // List of branches
   TBranch        *b_EventNumber;   //!
   TBranch        *b_mc_event_weight;   //!
   TBranch        *b_mu_n;   //!
   TBranch        *b_mu_E;   //!
   TBranch        *b_mu_pt;   //!
   TBranch        *b_mu_eta;   //!
   TBranch        *b_mu_phi;   //!
   TBranch        *b_mu_charge;   //!
   TBranch        *b_el_n;   //!
   TBranch        *b_el_E;   //!
   TBranch        *b_el_pt;   //!
   TBranch        *b_el_eta;   //!
   TBranch        *b_el_phi;   //!
   TBranch        *b_el_charge;   //!
   TBranch        *b_jet_n;   //!
   TBranch        *b_jet_E;   //!
   TBranch        *b_jet_m;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_tau2;   //!
   TBranch        *b_jet_tau1;   //!
   TBranch        *b_jet_tau3;   //!
   TBranch        *b_jet_tau4;   //!
   TBranch        *b_jet_isBtagged;   //!
   TBranch        *b_fatjet_n;   //!
   TBranch        *b_fatjet_E;   //!
   TBranch        *b_fatjet_m;   //!
   TBranch        *b_fatjet_pt;   //!
   TBranch        *b_fatjet_eta;   //!
   TBranch        *b_fatjet_phi;   //!
   TBranch        *b_fatjet_tau1;   //!
   TBranch        *b_fatjet_tau2;   //!
   TBranch        *b_fatjet_tau3;   //!
   TBranch        *b_fatjet_tau4;   //!
   TBranch        *b_fatjet_sj_n;   //!
   TBranch        *b_fatjet_m_sj1;   //!
   TBranch        *b_fatjet_pt_sj1;   //!
   TBranch        *b_fatjet_eta_sj1;   //!
   TBranch        *b_fatjet_phi_sj1;   //!
   TBranch        *b_fatjet_m_sj2;   //!
   TBranch        *b_fatjet_pt_sj2;   //!
   TBranch        *b_fatjet_eta_sj2;   //!
   TBranch        *b_fatjet_phi_sj2;   //!
   TBranch        *b_MET_ex;   //!
   TBranch        *b_MET_ey;   //!
   TBranch        *b_MET_et;   //!
   TBranch        *b_MET_phi;   //!
   TBranch        *b_MET_sumet;   //!
   TBranch        *b_MET_Truth_NonInt_ex;   //!
   TBranch        *b_MET_Truth_NonInt_ey;   //!
   TBranch        *b_MET_Truth_NonInt_et;   //!
   TBranch        *b_jet_truth_n;   //!
   TBranch        *b_jet_truth_pt;   //!
   TBranch        *b_jet_truth_eta;   //!
   TBranch        *b_jet_truth_phi;   //!
   TBranch        *b_jet_truth_m;   //!
   TBranch        *b_fatjet_truth_n;   //!
   TBranch        *b_fatjet_truth_pt;   //!
   TBranch        *b_fatjet_truth_eta;   //!
   TBranch        *b_fatjet_truth_phi;   //!
   TBranch        *b_fatjet_truth_m;   //!
   TBranch        *b_truthV_n;   //!
   TBranch        *b_truthV_pdgId;   //!
   TBranch        *b_truthV_pt;   //!
   TBranch        *b_truthV_m;   //!
   TBranch        *b_truthV_eta;   //!
   TBranch        *b_truthV_phi;   //!
   TBranch        *b_truthV_charge;   //!
   TBranch        *b_truthV_dau1_pdgId;   //!
   TBranch        *b_truthV_dau1_pt;   //!
   TBranch        *b_truthV_dau1_m;   //!
   TBranch        *b_truthV_dau1_eta;   //!
   TBranch        *b_truthV_dau1_phi;   //!
   TBranch        *b_truthV_dau1_charge;   //!
   TBranch        *b_truthV_dau2_pdgId;   //!
   TBranch        *b_truthV_dau2_pt;   //!
   TBranch        *b_truthV_dau2_m;   //!
   TBranch        *b_truthV_dau2_eta;   //!
   TBranch        *b_truthV_dau2_phi;   //!
   TBranch        *b_truthV_dau2_charge;   //!

void Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   TTree* fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("mc_event_weight", &mc_event_weight, &b_mc_event_weight);
   fChain->SetBranchAddress("mu_n", &mu_n, &b_mu_n);
   fChain->SetBranchAddress("mu_E", mu_E, &b_mu_E);
   fChain->SetBranchAddress("mu_pt", mu_pt, &b_mu_pt);
   fChain->SetBranchAddress("mu_eta", mu_eta, &b_mu_eta);
   fChain->SetBranchAddress("mu_phi", mu_phi, &b_mu_phi);
   fChain->SetBranchAddress("mu_charge", mu_charge, &b_mu_charge);
   fChain->SetBranchAddress("el_n", &el_n, &b_el_n);
   fChain->SetBranchAddress("el_E", el_E, &b_el_E);
   fChain->SetBranchAddress("el_pt", el_pt, &b_el_pt);
   fChain->SetBranchAddress("el_eta", el_eta, &b_el_eta);
   fChain->SetBranchAddress("el_phi", el_phi, &b_el_phi);
   fChain->SetBranchAddress("el_charge", el_charge, &b_el_charge);
   fChain->SetBranchAddress("jet_n", &jet_n, &b_jet_n);
   fChain->SetBranchAddress("jet_E", jet_E, &b_jet_E);
   fChain->SetBranchAddress("jet_m", jet_m, &b_jet_m);
   fChain->SetBranchAddress("jet_pt", jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_tau2", jet_tau2, &b_jet_tau2);
   fChain->SetBranchAddress("jet_tau1", jet_tau1, &b_jet_tau1);
   fChain->SetBranchAddress("jet_tau3", jet_tau3, &b_jet_tau3);
   fChain->SetBranchAddress("jet_tau4", jet_tau4, &b_jet_tau4);
   fChain->SetBranchAddress("jet_isBtagged", jet_isBtagged, &b_jet_isBtagged);
   fChain->SetBranchAddress("fatjet_n", &fatjet_n, &b_fatjet_n);
   fChain->SetBranchAddress("fatjet_E", fatjet_E, &b_fatjet_E);
   fChain->SetBranchAddress("fatjet_m", fatjet_m, &b_fatjet_m);
   fChain->SetBranchAddress("fatjet_pt", fatjet_pt, &b_fatjet_pt);
   fChain->SetBranchAddress("fatjet_eta", fatjet_eta, &b_fatjet_eta);
   fChain->SetBranchAddress("fatjet_phi", fatjet_phi, &b_fatjet_phi);
   fChain->SetBranchAddress("fatjet_tau1", fatjet_tau1, &b_fatjet_tau1);
   fChain->SetBranchAddress("fatjet_tau2", fatjet_tau2, &b_fatjet_tau2);
   fChain->SetBranchAddress("fatjet_tau3", fatjet_tau3, &b_fatjet_tau3);
   fChain->SetBranchAddress("fatjet_tau4", fatjet_tau4, &b_fatjet_tau4);
   fChain->SetBranchAddress("fatjet_sj_n", fatjet_sj_n, &b_fatjet_sj_n);
   fChain->SetBranchAddress("fatjet_m_sj1", fatjet_m_sj1, &b_fatjet_m_sj1);
   fChain->SetBranchAddress("fatjet_pt_sj1", fatjet_pt_sj1, &b_fatjet_pt_sj1);
   fChain->SetBranchAddress("fatjet_eta_sj1", fatjet_eta_sj1, &b_fatjet_eta_sj1);
   fChain->SetBranchAddress("fatjet_phi_sj1", fatjet_phi_sj1, &b_fatjet_phi_sj1);
   fChain->SetBranchAddress("fatjet_m_sj2", fatjet_m_sj2, &b_fatjet_m_sj2);
   fChain->SetBranchAddress("fatjet_pt_sj2", fatjet_pt_sj2, &b_fatjet_pt_sj2);
   fChain->SetBranchAddress("fatjet_eta_sj2", fatjet_eta_sj2, &b_fatjet_eta_sj2);
   fChain->SetBranchAddress("fatjet_phi_sj2", fatjet_phi_sj2, &b_fatjet_phi_sj2);
   fChain->SetBranchAddress("MET_ex", &MET_ex, &b_MET_ex);
   fChain->SetBranchAddress("MET_ey", &MET_ey, &b_MET_ey);
   fChain->SetBranchAddress("MET_et", &MET_et, &b_MET_et);
   fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   fChain->SetBranchAddress("MET_sumet", &MET_sumet, &b_MET_sumet);
   fChain->SetBranchAddress("MET_Truth_NonInt_ex", &MET_Truth_NonInt_ex, &b_MET_Truth_NonInt_ex);
   fChain->SetBranchAddress("MET_Truth_NonInt_ey", &MET_Truth_NonInt_ey, &b_MET_Truth_NonInt_ey);
   fChain->SetBranchAddress("MET_Truth_NonInt_et", &MET_Truth_NonInt_et, &b_MET_Truth_NonInt_et);
   fChain->SetBranchAddress("jet_truth_n", &jet_truth_n, &b_jet_truth_n);
   fChain->SetBranchAddress("jet_truth_pt", jet_truth_pt, &b_jet_truth_pt);
   fChain->SetBranchAddress("jet_truth_eta", jet_truth_eta, &b_jet_truth_eta);
   fChain->SetBranchAddress("jet_truth_phi", jet_truth_phi, &b_jet_truth_phi);
   fChain->SetBranchAddress("jet_truth_m", jet_truth_m, &b_jet_truth_m);
   fChain->SetBranchAddress("fatjet_truth_n", &fatjet_truth_n, &b_fatjet_truth_n);
   fChain->SetBranchAddress("fatjet_truth_pt", fatjet_truth_pt, &b_fatjet_truth_pt);
   fChain->SetBranchAddress("fatjet_truth_eta", fatjet_truth_eta, &b_fatjet_truth_eta);
   fChain->SetBranchAddress("fatjet_truth_phi", fatjet_truth_phi, &b_fatjet_truth_phi);
   fChain->SetBranchAddress("fatjet_truth_m", fatjet_truth_m, &b_fatjet_truth_m);
   fChain->SetBranchAddress("truthV_n", &truthV_n, &b_truthV_n);
   fChain->SetBranchAddress("truthV_pdgId", truthV_pdgId, &b_truthV_pdgId);
   fChain->SetBranchAddress("truthV_pt", truthV_pt, &b_truthV_pt);
   fChain->SetBranchAddress("truthV_m", truthV_m, &b_truthV_m);
   fChain->SetBranchAddress("truthV_eta", truthV_eta, &b_truthV_eta);
   fChain->SetBranchAddress("truthV_phi", truthV_phi, &b_truthV_phi);
   fChain->SetBranchAddress("truthV_charge", truthV_charge, &b_truthV_charge);
   fChain->SetBranchAddress("truthV_dau1_pdgId", truthV_dau1_pdgId, &b_truthV_dau1_pdgId);
   fChain->SetBranchAddress("truthV_dau1_pt", truthV_dau1_pt, &b_truthV_dau1_pt);
   fChain->SetBranchAddress("truthV_dau1_m", truthV_dau1_m, &b_truthV_dau1_m);
   fChain->SetBranchAddress("truthV_dau1_eta", truthV_dau1_eta, &b_truthV_dau1_eta);
   fChain->SetBranchAddress("truthV_dau1_phi", truthV_dau1_phi, &b_truthV_dau1_phi);
   fChain->SetBranchAddress("truthV_dau1_charge", truthV_dau1_charge, &b_truthV_dau1_charge);
   fChain->SetBranchAddress("truthV_dau2_pdgId", truthV_dau2_pdgId, &b_truthV_dau2_pdgId);
   fChain->SetBranchAddress("truthV_dau2_pt", truthV_dau2_pt, &b_truthV_dau2_pt);
   fChain->SetBranchAddress("truthV_dau2_m", truthV_dau2_m, &b_truthV_dau2_m);
   fChain->SetBranchAddress("truthV_dau2_eta", truthV_dau2_eta, &b_truthV_dau2_eta);
   fChain->SetBranchAddress("truthV_dau2_phi", truthV_dau2_phi, &b_truthV_dau2_phi);
   fChain->SetBranchAddress("truthV_dau2_charge", truthV_dau2_charge, &b_truthV_dau2_charge);
}
