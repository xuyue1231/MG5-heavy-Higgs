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
Int_t           ph_n;
Float_t         ph_pt[10];   //[mu_n]
Float_t         ph_eta[10];   //[mu_n]
Float_t         ph_phi[10];   //[mu_n]
Int_t           mu_n;
Float_t         mu_E[10];   //[mu_n]
Float_t         mu_pt[10];   //[mu_n]
Float_t         mu_eta[10];   //[mu_n]
Float_t         mu_phi[10];   //[mu_n]
Float_t         mu_E_cor[10];   //[mu_n]
Float_t         mu_pt_cor[10];   //[mu_n]
Float_t         mu_eta_cor[10];   //[mu_n]
Float_t         mu_phi_cor[10];   //[mu_n]
Int_t           mu_charge[10];   //[mu_n]
Float_t         mu_x[10];   //[mu_n]
Float_t         mu_y[10];   //[mu_n]
Float_t         mu_z[10];   //[mu_n]
Float_t         mu_xd[10];   //[mu_n]
Float_t         mu_yd[10];   //[mu_n]
Float_t         mu_zd[10];   //[mu_n]
Float_t  	mu_truth_dr[10];   //[mu_n]
Float_t  	mu_truth_E_cor[10];   //[mu_n]
Float_t  	mu_truth_pt_cor[10];   //[mu_n]
Float_t  	mu_truth_eta_cor[10];   //[mu_n]
Float_t  	mu_truth_phi_cor[10];   //[mu_n]
Float_t  	mu_truth_E[10];   //[mu_n]
Float_t  	mu_truth_pt[10];   //[mu_n]
Float_t  	mu_truth_eta[10];   //[mu_n]
Float_t  	mu_truth_phi[10];   //[mu_n]
Int_t    	mu_truth_pdgIdMoth[10];   //[mu_n]
Int_t           el_n;
Float_t         el_E[10];   //[el_n]
Float_t         el_pt[10];   //[el_n]
Float_t         el_eta[10];   //[el_n]
Float_t         el_phi[10];   //[el_n]
Float_t         el_E_cor[10];   //[el_n]
Float_t         el_pt_cor[10];   //[el_n]
Float_t         el_eta_cor[10];   //[el_n]
Float_t         el_phi_cor[10];   //[el_n]
Int_t           el_charge[10];   //[el_n]
Float_t         el_x[10];   //[el_n]
Float_t         el_y[10];   //[el_n]
Float_t         el_z[10];   //[el_n]
Float_t         el_xd[10];   //[el_n]
Float_t         el_yd[10];   //[el_n]
Float_t         el_zd[10];   //[el_n]
Float_t  	el_truth_dr[10];   //[el_n]
Float_t  	el_truth_E_cor[10];   //[el_n]
Float_t  	el_truth_pt_cor[10];   //[el_n]
Float_t  	el_truth_eta_cor[10];   //[el_n]
Float_t  	el_truth_phi_cor[10];   //[el_n]
Float_t  	el_truth_E[10];   //[el_n]
Float_t  	el_truth_pt[10];   //[el_n]
Float_t  	el_truth_eta[10];   //[el_n]
Float_t  	el_truth_phi[10];   //[el_n]
Int_t   	el_truth_pdgIdMoth[10];   //[el_n]
Int_t           tau_n;
Float_t         tau_E[10];   //[tau_n]
Float_t         tau_pt[10];   //[tau_n]
Float_t         tau_eta[10];   //[tau_n]
Float_t         tau_phi[10];   //[tau_n]
Int_t           tau_charge[10];   //[tau_n]
Int_t           tau_numTrack[10];   //[tau_n]
Int_t           tau_numNeut[10];   //[tau_n]
Float_t         tau_trk1_pt[10];   //[tau_n]
Float_t         tau_trk1_eta[10];   //[tau_n]
Float_t         tau_trk1_phi[10];   //[tau_n]
Float_t         tau_trk1_m[10];   //[tau_n]
Float_t         tau_trk1_x[10];   //[tau_n]
Float_t         tau_trk1_y[10];   //[tau_n]
Float_t         tau_trk1_z[10];   //[tau_n]
Float_t         tau_trk1_xd[10];   //[tau_n]
Float_t         tau_trk1_yd[10];   //[tau_n]
Float_t         tau_trk1_zd[10];   //[tau_n]
Int_t           tau_trk1_charge[10];   //[tau_n]
Int_t           tau_trk1_pdgId[10];   //[tau_n]
Float_t         tau_trk2_pt[10];   //[tau_n]
Float_t         tau_trk2_eta[10];   //[tau_n]
Float_t         tau_trk2_phi[10];   //[tau_n]
Float_t         tau_trk2_m[10];   //[tau_n]
Float_t         tau_trk2_x[10];   //[tau_n]
Float_t         tau_trk2_y[10];   //[tau_n]
Float_t         tau_trk2_z[10];   //[tau_n]
Float_t         tau_trk2_xd[10];   //[tau_n]
Float_t         tau_trk2_yd[10];   //[tau_n]
Float_t         tau_trk2_zd[10];   //[tau_n]
Int_t           tau_trk2_charge[10];   //[tau_n]
Int_t           tau_trk2_pdgId[10];   //[tau_n]
Float_t         tau_trk3_pt[10];   //[tau_n]
Float_t         tau_trk3_eta[10];   //[tau_n]
Float_t         tau_trk3_phi[10];   //[tau_n]
Float_t         tau_trk3_m[10];   //[tau_n]
Float_t         tau_trk3_x[10];   //[tau_n]
Float_t         tau_trk3_y[10];   //[tau_n]
Float_t         tau_trk3_z[10];   //[tau_n]
Float_t         tau_trk3_xd[10];   //[tau_n]
Float_t         tau_trk3_yd[10];   //[tau_n]
Float_t         tau_trk3_zd[10];   //[tau_n]
Int_t           tau_trk3_charge[10];   //[tau_n]
Int_t           tau_trk3_pdgId[10];   //[tau_n]
Int_t           jet_n;
Float_t         jet_E[50];   //[jet_n]
Float_t         jet_m[50];   //[jet_n]
Float_t         jet_pt[50];   //[jet_n]
Float_t         jet_eta[50];   //[jet_n]
Float_t         jet_phi[50];   //[jet_n]
Bool_t          jet_isTau[50];   //[jet_n]
Bool_t          jet_isBtagged[50];   //[jet_n]
Float_t         MET_ex;
Float_t         MET_ey;
Float_t         MET_ez;
Float_t         MET_et;
Float_t         MET_phi;
Float_t         MET_sumet;
Float_t         MET_Truth_NonInt_ex;
Float_t         MET_Truth_NonInt_ey;
Float_t         MET_Truth_NonInt_ez;
Float_t         MET_Truth_NonInt_et;
Int_t           mc_n;
Float_t         mc_pt[20];   //[mc_n]
Float_t         mc_m[20];   //[mc_n]
Float_t         mc_eta[20];   //[mc_n]
Float_t         mc_phi[20];   //[mc_n]
Int_t           mc_status[20];   //[mc_n]
Int_t           mc_pdgId[20];   //[mc_n]
Int_t           mc_pdgIdMoth[20];   //[mc_n]
Int_t           mc_charge[20];   //[mc_n]
Int_t           jet_truth_n;
Float_t         jet_truth_pt[50];   //[jet_truth_n]
Float_t         jet_truth_eta[50];   //[jet_truth_n]
Float_t         jet_truth_phi[50];   //[jet_truth_n]
Float_t         jet_truth_m[50];   //[jet_truth_n]
Int_t           truthTau_n;
Int_t           truthTau_vis_pdgId[10];   //[truthTau_n]
Int_t           truthTau_vis_pdgIdMoth[10];   //[truthTau_n]
Float_t         truthTau_vis_pt[10];   //[truthTau_n]
Float_t         truthTau_vis_m[10];   //[truthTau_n]
Float_t         truthTau_vis_eta[10];   //[truthTau_n]
Float_t         truthTau_vis_phi[10];   //[truthTau_n]
Int_t           truthTau_vis_charge[10];   //[truthTau_n]
Int_t           truthTau_vis_numTrack[10];   //[truthTau_n]
Int_t           truthTau_vis_numNeut[10];   //[truthTau_n]
Int_t           truthTau_vis_numPi0[10];   //[truthTau_n]
Int_t           truthTau_vis_numGam[10];   //[truthTau_n]
Float_t         truthTau_nu3_pt[10];   //[truthTau_n]
Float_t         truthTau_nu3_m[10];   //[truthTau_n]
Float_t         truthTau_nu3_eta[10];   //[truthTau_n]
Float_t         truthTau_nu3_phi[10];   //[truthTau_n]
Float_t         truthTau_inv_pt[10];   //[truthTau_n]
Float_t         truthTau_inv_m[10];   //[truthTau_n]
Float_t         truthTau_inv_eta[10];   //[truthTau_n]
Float_t         truthTau_inv_phi[10];   //[truthTau_n]
Float_t         truthTau_trk1_pt[10];   //[truthTau_n]
Float_t         truthTau_trk1_eta[10];   //[truthTau_n]
Float_t         truthTau_trk1_phi[10];   //[truthTau_n]
Float_t         truthTau_trk1_m[10];   //[truthTau_n]
Float_t         truthTau_trk1_x[10];   //[truthTau_n]
Float_t         truthTau_trk1_y[10];   //[truthTau_n]
Float_t         truthTau_trk1_z[10];   //[truthTau_n]
Float_t         truthTau_trk1_xd[10];   //[truthTau_n]
Float_t         truthTau_trk1_yd[10];   //[truthTau_n]
Float_t         truthTau_trk1_zd[10];   //[truthTau_n]
Int_t           truthTau_trk1_charge[10];   //[tau_n]
Int_t           truthTau_trk1_pdgId[10];   //[truthTau_n]
Float_t         truthTau_trk2_pt[10];   //[truthTau_n]
Float_t         truthTau_trk2_eta[10];   //[truthTau_n]
Float_t         truthTau_trk2_phi[10];   //[truthTau_n]
Float_t         truthTau_trk2_m[10];   //[truthTau_n]
Float_t         truthTau_trk2_x[10];   //[truthTau_n]
Float_t         truthTau_trk2_y[10];   //[truthTau_n]
Float_t         truthTau_trk2_z[10];   //[truthTau_n]
Float_t         truthTau_trk2_xd[10];   //[truthTau_n]
Float_t         truthTau_trk2_yd[10];   //[truthTau_n]
Float_t         truthTau_trk2_zd[10];   //[truthTau_n]
Int_t           truthTau_trk2_charge[10];   //[tau_n]
Int_t           truthTau_trk2_pdgId[10];   //[truthTau_n]
Float_t         truthTau_trk3_pt[10];   //[truthTau_n]
Float_t         truthTau_trk3_eta[10];   //[truthTau_n]
Float_t         truthTau_trk3_phi[10];   //[truthTau_n]
Float_t         truthTau_trk3_m[10];   //[truthTau_n]
Float_t         truthTau_trk3_x[10];   //[truthTau_n]
Float_t         truthTau_trk3_y[10];   //[truthTau_n]
Float_t         truthTau_trk3_z[10];   //[truthTau_n]
Float_t         truthTau_trk3_xd[10];   //[truthTau_n]
Float_t         truthTau_trk3_yd[10];   //[truthTau_n]
Float_t         truthTau_trk3_zd[10];   //[truthTau_n]
Int_t           truthTau_trk3_charge[10];   //[tau_n]
Int_t           truthTau_trk3_pdgId[10];   //[truthTau_n]
Float_t         truthTau_neut1_e[10];   //[truthTau_n]
Float_t         truthTau_neut1_pt[10];   //[truthTau_n]
Float_t         truthTau_neut1_eta[10];   //[truthTau_n]
Float_t         truthTau_neut1_phi[10];   //[truthTau_n]
Int_t           truthTau_neut1_pdgId[10];   //[truthTau_n]
Float_t         truthTau_neut2_e[10];   //[truthTau_n]
Float_t         truthTau_neut2_pt[10];   //[truthTau_n]
Float_t         truthTau_neut2_eta[10];   //[truthTau_n]
Float_t         truthTau_neut2_phi[10];   //[truthTau_n]
Int_t           truthTau_neut2_pdgId[10];   //[truthTau_n]
Float_t         truthTau_neut3_e[10];   //[truthTau_n]
Float_t         truthTau_neut3_pt[10];   //[truthTau_n]
Float_t         truthTau_neut3_eta[10];   //[truthTau_n]
Float_t         truthTau_neut3_phi[10];   //[truthTau_n]
Int_t           truthTau_neut3_pdgId[10];   //[truthTau_n]

void SetBranch(TTree* tree) {
  tree->Branch("EventNumber", &EventNumber,"EventNumber/I");
  tree->Branch("mc_event_weight", &mc_event_weight,"mc_event_weight/F");
  //tree->Branch("mc_event_ScalePDF", &mc_event_ScalePDF,"mc_event_ScalePDF/F");
  //tree->Branch("mc_event_AlphaQED", &mc_event_AlphaQED,"mc_event_AlphaQED/F");
  //tree->Branch("mc_event_AlphaQCD", &mc_event_AlphaQCD,"mc_event_AlphaQCD/F");
  tree->Branch("ph_n", &ph_n,"ph_n/I");
  tree->Branch("ph_pt", ph_pt,"ph_pt[ph_n]/F");
  tree->Branch("ph_eta", ph_eta,"ph_eta[ph_n]/F");
  tree->Branch("ph_phi", ph_phi,"ph_phi[ph_n]/F");
  tree->Branch("mu_n", &mu_n,"mu_n/I");
  tree->Branch("mu_E", mu_E,"mu_E[mu_n]/F");
  tree->Branch("mu_pt", mu_pt,"mu_pt[mu_n]/F");
  tree->Branch("mu_eta", mu_eta,"mu_eta[mu_n]/F");
  tree->Branch("mu_phi", mu_phi,"mu_phi[mu_n]/F");
  tree->Branch("mu_E_cor", mu_E_cor,"mu_E_cor[mu_n]/F");
  tree->Branch("mu_pt_cor", mu_pt_cor,"mu_pt_cor[mu_n]/F");
  tree->Branch("mu_eta_cor", mu_eta_cor,"mu_eta_cor[mu_n]/F");
  tree->Branch("mu_phi_cor", mu_phi_cor,"mu_phi_cor[mu_n]/F");
  tree->Branch("mu_charge", mu_charge,"mu_charge[mu_n]/I");
  tree->Branch("mu_x", mu_x,"mu_x[mu_n]/F");
  tree->Branch("mu_y", mu_y,"mu_y[mu_n]/F");
  tree->Branch("mu_z", mu_z,"mu_z[mu_n]/F");
  tree->Branch("mu_xd", mu_xd,"mu_xd[mu_n]/F");
  tree->Branch("mu_yd", mu_yd,"mu_yd[mu_n]/F");
  tree->Branch("mu_zd", mu_zd,"mu_zd[mu_n]/F");
  tree->Branch("mu_truth_dr", mu_truth_dr,"mu_truth_dr[mu_n]/F");
  tree->Branch("mu_truth_E_cor", mu_truth_E_cor,"mu_truth_E_cor[mu_n]/F");
  tree->Branch("mu_truth_pt_cor", mu_truth_pt_cor,"mu_truth_pt_cor[mu_n]/F");
  tree->Branch("mu_truth_eta_cor", mu_truth_eta_cor,"mu_truth_eta_cor[mu_n]/F");
  tree->Branch("mu_truth_phi_cor", mu_truth_phi_cor,"mu_truth_phi_cor[mu_n]/F");
  tree->Branch("mu_truth_E", mu_truth_E,"mu_truth_E[mu_n]/F");
  tree->Branch("mu_truth_pt", mu_truth_pt,"mu_truth_pt[mu_n]/F");
  tree->Branch("mu_truth_eta", mu_truth_eta,"mu_truth_eta[mu_n]/F");
  tree->Branch("mu_truth_phi", mu_truth_phi,"mu_truth_phi[mu_n]/F");
  tree->Branch("mu_truth_pdgIdMoth", mu_truth_pdgIdMoth,"mu_truth_pdgIdMoth[mu_n]/I");
  tree->Branch("el_n", &el_n,"el_n/I");
  tree->Branch("el_E", el_E,"el_E[el_n]/F");
  tree->Branch("el_pt", el_pt,"el_pt[el_n]/F");
  tree->Branch("el_eta", el_eta,"el_eta[el_n]/F");
  tree->Branch("el_phi", el_phi,"el_phi[el_n]/F");
  tree->Branch("el_E_cor", el_E_cor,"el_E_cor[el_n]/F");
  tree->Branch("el_pt_cor", el_pt_cor,"el_pt_cor[el_n]/F");
  tree->Branch("el_eta_cor", el_eta_cor,"el_eta_cor[el_n]/F");
  tree->Branch("el_phi_cor", el_phi_cor,"el_phi_cor[el_n]/F");
  tree->Branch("el_charge", el_charge,"el_charge[el_n]/I");
  tree->Branch("el_x", el_x,"el_x[el_n]/F");
  tree->Branch("el_y", el_y,"el_y[el_n]/F");
  tree->Branch("el_z", el_z,"el_z[el_n]/F");
  tree->Branch("el_xd", el_xd,"el_xd[el_n]/F");
  tree->Branch("el_yd", el_yd,"el_yd[el_n]/F");
  tree->Branch("el_zd", el_zd,"el_zd[el_n]/F");
  tree->Branch("el_truth_dr", el_truth_dr,"el_truth_dr[el_n]/F");
  tree->Branch("el_truth_E_cor", el_truth_E_cor,"el_truth_E_cor[el_n]/F");
  tree->Branch("el_truth_pt_cor", el_truth_pt_cor,"el_truth_pt_cor[el_n]/F");
  tree->Branch("el_truth_eta_cor", el_truth_eta_cor,"el_truth_eta_cor[el_n]/F");
  tree->Branch("el_truth_phi_cor", el_truth_phi_cor,"el_truth_phi_cor[el_n]/F");
  tree->Branch("el_truth_E", el_truth_E,"el_truth_E[el_n]/F");
  tree->Branch("el_truth_pt", el_truth_pt,"el_truth_pt[el_n]/F");
  tree->Branch("el_truth_eta", el_truth_eta,"el_truth_eta[el_n]/F");
  tree->Branch("el_truth_phi", el_truth_phi,"el_truth_phi[el_n]/F");
  tree->Branch("el_truth_pdgIdMoth", el_truth_pdgIdMoth,"el_truth_pdgIdMoth[el_n]/I");
  tree->Branch("tau_n", &tau_n,"tau_n/I");
  tree->Branch("tau_E", tau_E,"tau_E[tau_n]/F");
  tree->Branch("tau_pt", tau_pt,"tau_pt[tau_n]/F");
  tree->Branch("tau_eta", tau_eta,"tau_eta[tau_n]/F");
  tree->Branch("tau_phi", tau_phi,"tau_phi[tau_n]/F");
  tree->Branch("tau_charge", tau_charge,"tau_charge[tau_n]/I");
  tree->Branch("tau_numTrack", tau_numTrack,"tau_numTrack[tau_n]/I");
  tree->Branch("tau_numNeut", tau_numNeut,"tau_numNeut[tau_n]/I");
  tree->Branch("tau_trk1_pt", tau_trk1_pt,"tau_trk1_pt[tau_n]/F");
  tree->Branch("tau_trk1_eta", tau_trk1_eta,"tau_trk1_eta[tau_n]/F");
  tree->Branch("tau_trk1_phi", tau_trk1_phi,"tau_trk1_phi[tau_n]/F");
  tree->Branch("tau_trk1_m", tau_trk1_m,"tau_trk1_m[tau_n]/F");
  tree->Branch("tau_trk1_x", tau_trk1_x,"tau_trk1_x[tau_n]/F");
  tree->Branch("tau_trk1_y", tau_trk1_y,"tau_trk1_y[tau_n]/F");
  tree->Branch("tau_trk1_z", tau_trk1_z,"tau_trk1_z[tau_n]/F");
  tree->Branch("tau_trk1_xd", tau_trk1_xd,"tau_trk1_xd[tau_n]/F");
  tree->Branch("tau_trk1_yd", tau_trk1_yd,"tau_trk1_yd[tau_n]/F");
  tree->Branch("tau_trk1_zd", tau_trk1_zd,"tau_trk1_zd[tau_n]/F");
  tree->Branch("tau_trk1_charge", tau_trk1_charge,"tau_trk1_charge[tau_n]/I");
  tree->Branch("tau_trk1_pdgId", tau_trk1_pdgId,"tau_trk1_pdgId[tau_n]/I");
  tree->Branch("tau_trk2_pt", tau_trk2_pt,"tau_trk2_pt[tau_n]/F");
  tree->Branch("tau_trk2_eta", tau_trk2_eta,"tau_trk2_eta[tau_n]/F");
  tree->Branch("tau_trk2_phi", tau_trk2_phi,"tau_trk2_phi[tau_n]/F");
  tree->Branch("tau_trk2_m", tau_trk2_m,"tau_trk2_m[tau_n]/F");
  tree->Branch("tau_trk2_x", tau_trk2_x,"tau_trk2_x[tau_n]/F");
  tree->Branch("tau_trk2_y", tau_trk2_y,"tau_trk2_y[tau_n]/F");
  tree->Branch("tau_trk2_z", tau_trk2_z,"tau_trk2_z[tau_n]/F");
  tree->Branch("tau_trk2_xd", tau_trk2_xd,"tau_trk2_xd[tau_n]/F");
  tree->Branch("tau_trk2_yd", tau_trk2_yd,"tau_trk2_yd[tau_n]/F");
  tree->Branch("tau_trk2_zd", tau_trk2_zd,"tau_trk2_zd[tau_n]/F");
  tree->Branch("tau_trk2_charge", tau_trk2_charge,"tau_trk2_charge[tau_n]/I");
  tree->Branch("tau_trk2_pdgId", tau_trk2_pdgId,"tau_trk2_pdgId[tau_n]/I");
  tree->Branch("tau_trk3_pt", tau_trk3_pt,"tau_trk3_pt[tau_n]/F");
  tree->Branch("tau_trk3_eta", tau_trk3_eta,"tau_trk3_eta[tau_n]/F");
  tree->Branch("tau_trk3_phi", tau_trk3_phi,"tau_trk3_phi[tau_n]/F");
  tree->Branch("tau_trk3_m", tau_trk3_m,"tau_trk3_m[tau_n]/F");
  tree->Branch("tau_trk3_x", tau_trk3_x,"tau_trk3_x[tau_n]/F");
  tree->Branch("tau_trk3_y", tau_trk3_y,"tau_trk3_y[tau_n]/F");
  tree->Branch("tau_trk3_z", tau_trk3_z,"tau_trk3_z[tau_n]/F");
  tree->Branch("tau_trk3_xd", tau_trk3_xd,"tau_trk3_xd[tau_n]/F");
  tree->Branch("tau_trk3_yd", tau_trk3_yd,"tau_trk3_yd[tau_n]/F");
  tree->Branch("tau_trk3_zd", tau_trk3_zd,"tau_trk3_zd[tau_n]/F");
  tree->Branch("tau_trk3_charge", tau_trk3_charge,"tau_trk3_charge[tau_n]/I");
  tree->Branch("tau_trk3_pdgId", tau_trk3_pdgId,"tau_trk3_pdgId[tau_n]/I");
  tree->Branch("jet_n", &jet_n,"jet_n/I");
  tree->Branch("jet_E", jet_E,"jet_E[jet_n]/F");
  tree->Branch("jet_m", jet_m,"jet_m[jet_n]/F");
  tree->Branch("jet_pt", jet_pt,"jet_pt[jet_n]/F");
  tree->Branch("jet_eta", jet_eta,"jet_eta[jet_n]/F");
  tree->Branch("jet_phi", jet_phi,"jet_phi[jet_n]/F");
  tree->Branch("jet_isTau", jet_isTau,"jet_isTau[jet_n]/O");
  tree->Branch("jet_isBtagged", jet_isBtagged,"jet_isBtagged[jet_n]/O");
  tree->Branch("MET_ex", &MET_ex,"MET_ex/F");
  tree->Branch("MET_ey", &MET_ey,"MET_ey/F");
  tree->Branch("MET_ez", &MET_ez,"MET_ez/F");
  tree->Branch("MET_et", &MET_et,"MET_et/F");
  tree->Branch("MET_phi", &MET_phi,"MET_phi/F");
  tree->Branch("MET_sumet", &MET_sumet,"MET_sumet/F");
  tree->Branch("MET_Truth_NonInt_ex", &MET_Truth_NonInt_ex,"MET_Truth_NonInt_ex/F");
  tree->Branch("MET_Truth_NonInt_ey", &MET_Truth_NonInt_ey,"MET_Truth_NonInt_ey/F");
  tree->Branch("MET_Truth_NonInt_ez", &MET_Truth_NonInt_ez,"MET_Truth_NonInt_ez/F");
  tree->Branch("MET_Truth_NonInt_et", &MET_Truth_NonInt_et,"MET_Truth_NonInt_et/F");
  tree->Branch("mc_n", &mc_n,"mc_n/I");
  tree->Branch("mc_pt", mc_pt,"mc_pt[mc_n]/F");
  tree->Branch("mc_m", mc_m,"mc_m[mc_n]/F");
  tree->Branch("mc_eta", mc_eta,"mc_eta[mc_n]/F");
  tree->Branch("mc_phi", mc_phi,"mc_phi[mc_n]/F");
  tree->Branch("mc_status", mc_status,"mc_status[mc_n]/I");
  tree->Branch("mc_pdgId", mc_pdgId,"mc_pdgId[mc_n]/I");
  tree->Branch("mc_pdgIdMoth", mc_pdgIdMoth,"mc_pdgIdMoth[mc_n]/I");
  tree->Branch("mc_charge", mc_charge,"mc_charge[mc_n]/I");
  tree->Branch("jet_truth_n", &jet_truth_n,"jet_truth_n/I");
  tree->Branch("jet_truth_pt", jet_truth_pt,"jet_truth_pt[jet_truth_n]/F");
  tree->Branch("jet_truth_eta", jet_truth_eta,"jet_truth_eta[jet_truth_n]/F");
  tree->Branch("jet_truth_phi", jet_truth_phi,"jet_truth_phi[jet_truth_n]/F");
  tree->Branch("jet_truth_m", jet_truth_m,"jet_truth_m[jet_truth_n]/F");
  tree->Branch("truthTau_n", &truthTau_n,"truthTau_n/I");
  tree->Branch("truthTau_vis_pdgId", truthTau_vis_pdgId,"truthTau_vis_pdgId[truthTau_n]/I");
  tree->Branch("truthTau_vis_pdgIdMoth", truthTau_vis_pdgIdMoth,"truthTau_vis_pdgIdMoth[truthTau_n]/I");
  tree->Branch("truthTau_vis_pt", truthTau_vis_pt,"truthTau_vis_pt[truthTau_n]/F");
  tree->Branch("truthTau_vis_m", truthTau_vis_m,"truthTau_vis_m[truthTau_n]/F");
  tree->Branch("truthTau_vis_eta", truthTau_vis_eta,"truthTau_vis_eta[truthTau_n]/F");
  tree->Branch("truthTau_vis_phi", truthTau_vis_phi,"truthTau_vis_phi[truthTau_n]/F");
  tree->Branch("truthTau_vis_charge", truthTau_vis_charge,"truthTau_vis_charge[truthTau_n]/I");
  tree->Branch("truthTau_vis_numTrack", truthTau_vis_numTrack,"truthTau_vis_numTrack[truthTau_n]/I");
  tree->Branch("truthTau_vis_numNeut", truthTau_vis_numNeut,"truthTau_vis_numNeut[truthTau_n]/I");
  tree->Branch("truthTau_vis_numPi0", truthTau_vis_numPi0,"truthTau_vis_numPi0[truthTau_n]/I");
  tree->Branch("truthTau_vis_numGam", truthTau_vis_numGam,"truthTau_vis_numGam[truthTau_n]/I");
  tree->Branch("truthTau_nu3_pt", truthTau_nu3_pt,"truthTau_nu3_pt[truthTau_n]/F");
  tree->Branch("truthTau_nu3_m", truthTau_nu3_m,"truthTau_nu3_m[truthTau_n]/F");
  tree->Branch("truthTau_nu3_eta", truthTau_nu3_eta,"truthTau_nu3_eta[truthTau_n]/F");
  tree->Branch("truthTau_nu3_phi", truthTau_nu3_phi,"truthTau_nu3_phi[truthTau_n]/F");
  tree->Branch("truthTau_inv_pt", truthTau_inv_pt,"truthTau_inv_pt[truthTau_n]/F");
  tree->Branch("truthTau_inv_m", truthTau_inv_m,"truthTau_inv_m[truthTau_n]/F");
  tree->Branch("truthTau_inv_eta", truthTau_inv_eta,"truthTau_inv_eta[truthTau_n]/F");
  tree->Branch("truthTau_inv_phi", truthTau_inv_phi,"truthTau_inv_phi[truthTau_n]/F");
  tree->Branch("truthTau_trk1_pt", truthTau_trk1_pt,"truthTau_trk1_pt[truthTau_n]/F");
  tree->Branch("truthTau_trk1_eta", truthTau_trk1_eta,"truthTau_trk1_eta[truthTau_n]/F");
  tree->Branch("truthTau_trk1_phi", truthTau_trk1_phi,"truthTau_trk1_phi[truthTau_n]/F");
  tree->Branch("truthTau_trk1_m", truthTau_trk1_m,"truthTau_trk1_m[truthTau_n]/F");
  tree->Branch("truthTau_trk1_x", truthTau_trk1_x,"truthTau_trk1_x[truthTau_n]/F");
  tree->Branch("truthTau_trk1_y", truthTau_trk1_y,"truthTau_trk1_y[truthTau_n]/F");
  tree->Branch("truthTau_trk1_z", truthTau_trk1_z,"truthTau_trk1_z[truthTau_n]/F");
  tree->Branch("truthTau_trk1_xd", truthTau_trk1_xd,"truthTau_trk1_xd[truthTau_n]/F");
  tree->Branch("truthTau_trk1_yd", truthTau_trk1_yd,"truthTau_trk1_yd[truthTau_n]/F");
  tree->Branch("truthTau_trk1_zd", truthTau_trk1_zd,"truthTau_trk1_zd[truthTau_n]/F");
  tree->Branch("truthTau_trk1_charge", truthTau_trk1_charge,"truthTau_trk1_charge[truthTau_n]/I");
  tree->Branch("truthTau_trk1_pdgId", truthTau_trk1_pdgId,"truthTau_trk1_pdgId[truthTau_n]/I");
  tree->Branch("truthTau_trk2_pt", truthTau_trk2_pt,"truthTau_trk2_pt[truthTau_n]/F");
  tree->Branch("truthTau_trk2_eta", truthTau_trk2_eta,"truthTau_trk2_eta[truthTau_n]/F");
  tree->Branch("truthTau_trk2_phi", truthTau_trk2_phi,"truthTau_trk2_phi[truthTau_n]/F");
  tree->Branch("truthTau_trk2_m", truthTau_trk2_m,"truthTau_trk2_m[truthTau_n]/F");
  tree->Branch("truthTau_trk2_x", truthTau_trk2_x,"truthTau_trk2_x[truthTau_n]/F");
  tree->Branch("truthTau_trk2_y", truthTau_trk2_y,"truthTau_trk2_y[truthTau_n]/F");
  tree->Branch("truthTau_trk2_z", truthTau_trk2_z,"truthTau_trk2_z[truthTau_n]/F");
  tree->Branch("truthTau_trk2_xd", truthTau_trk2_xd,"truthTau_trk2_xd[truthTau_n]/F");
  tree->Branch("truthTau_trk2_yd", truthTau_trk2_yd,"truthTau_trk2_yd[truthTau_n]/F");
  tree->Branch("truthTau_trk2_zd", truthTau_trk2_zd,"truthTau_trk2_zd[truthTau_n]/F");
  tree->Branch("truthTau_trk2_charge", truthTau_trk2_charge,"truthTau_trk2_charge[truthTau_n]/I");
  tree->Branch("truthTau_trk2_pdgId", truthTau_trk2_pdgId,"truthTau_trk2_pdgId[truthTau_n]/I");
  tree->Branch("truthTau_trk3_pt", truthTau_trk3_pt,"truthTau_trk3_pt[truthTau_n]/F");
  tree->Branch("truthTau_trk3_eta", truthTau_trk3_eta,"truthTau_trk3_eta[truthTau_n]/F");
  tree->Branch("truthTau_trk3_phi", truthTau_trk3_phi,"truthTau_trk3_phi[truthTau_n]/F");
  tree->Branch("truthTau_trk3_m", truthTau_trk3_m,"truthTau_trk3_m[truthTau_n]/F");
  tree->Branch("truthTau_trk3_x", truthTau_trk3_x,"truthTau_trk3_x[truthTau_n]/F");
  tree->Branch("truthTau_trk3_y", truthTau_trk3_y,"truthTau_trk3_y[truthTau_n]/F");
  tree->Branch("truthTau_trk3_z", truthTau_trk3_z,"truthTau_trk3_z[truthTau_n]/F");
  tree->Branch("truthTau_trk3_xd", truthTau_trk3_xd,"truthTau_trk3_xd[truthTau_n]/F");
  tree->Branch("truthTau_trk3_yd", truthTau_trk3_yd,"truthTau_trk3_yd[truthTau_n]/F");
  tree->Branch("truthTau_trk3_zd", truthTau_trk3_zd,"truthTau_trk3_zd[truthTau_n]/F");
  tree->Branch("truthTau_trk3_charge", truthTau_trk3_charge,"truthTau_trk3_charge[truthTau_n]/I");
  tree->Branch("truthTau_trk3_pdgId", truthTau_trk3_pdgId,"truthTau_trk3_pdgId[truthTau_n]/I");
  tree->Branch("truthTau_neut1_e", truthTau_neut1_e,"truthTau_neut1_e[truthTau_n]/F");
  tree->Branch("truthTau_neut1_pt", truthTau_neut1_pt,"truthTau_neut1_pt[truthTau_n]/F");
  tree->Branch("truthTau_neut1_eta", truthTau_neut1_eta,"truthTau_neut1_eta[truthTau_n]/F");
  tree->Branch("truthTau_neut1_phi", truthTau_neut1_phi,"truthTau_neut1_phi[truthTau_n]/F");
  tree->Branch("truthTau_neut1_pdgId", truthTau_neut1_pdgId,"truthTau_neut1_pdgId[truthTau_n]/I");
  tree->Branch("truthTau_neut2_e", truthTau_neut2_e,"truthTau_neut2_e[truthTau_n]/F");
  tree->Branch("truthTau_neut2_pt", truthTau_neut2_pt,"truthTau_neut2_pt[truthTau_n]/F");
  tree->Branch("truthTau_neut2_eta", truthTau_neut2_eta,"truthTau_neut2_eta[truthTau_n]/F");
  tree->Branch("truthTau_neut2_phi", truthTau_neut2_phi,"truthTau_neut2_phi[truthTau_n]/F");
  tree->Branch("truthTau_neut2_pdgId", truthTau_neut2_pdgId,"truthTau_neut2_pdgId[truthTau_n]/I");
  tree->Branch("truthTau_neut3_e", truthTau_neut3_e,"truthTau_neut3_e[truthTau_n]/F");
  tree->Branch("truthTau_neut3_pt", truthTau_neut3_pt,"truthTau_neut3_pt[truthTau_n]/F");
  tree->Branch("truthTau_neut3_eta", truthTau_neut3_eta,"truthTau_neut3_eta[truthTau_n]/F");
  tree->Branch("truthTau_neut3_phi", truthTau_neut3_phi,"truthTau_neut3_phi[truthTau_n]/F");
  tree->Branch("truthTau_neut3_pdgId", truthTau_neut3_pdgId,"truthTau_neut3_pdgId[truthTau_n]/I");
}

void RemoveEle(Int_t i) {
  for(Int_t j=i; j<el_n-1; j++) {
    el_E[j] = el_E[j+1];
    el_pt[j] = el_pt[j+1];
    el_eta[j] = el_eta[j+1];
    el_phi[j] = el_phi[j+1];
    el_E_cor[j] = el_E_cor[j+1];
    el_pt_cor[j] = el_pt_cor[j+1];
    el_eta_cor[j] = el_eta_cor[j+1];
    el_phi_cor[j] = el_phi_cor[j+1];
    el_charge[j] = el_charge[j+1];
    el_x[j] = el_x[j+1];
    el_y[j] = el_y[j+1];
    el_z[j] = el_z[j+1];
    el_xd[j] = el_xd[j+1];
    el_yd[j] = el_yd[j+1];
    el_zd[j] = el_zd[j+1];
    el_truth_dr[j] = el_truth_dr[j+1];
    el_truth_E_cor[j] = el_truth_E_cor[j+1];
    el_truth_pt_cor[j] = el_truth_pt_cor[j+1];
    el_truth_eta_cor[j] = el_truth_eta_cor[j+1];
    el_truth_phi_cor[j] = el_truth_phi_cor[j+1];
    el_truth_E[j] = el_truth_E[j+1];
    el_truth_pt[j] = el_truth_pt[j+1];
    el_truth_eta[j] = el_truth_eta[j+1];
    el_truth_phi[j] = el_truth_phi[j+1];
    el_truth_pdgIdMoth[j] = el_truth_pdgIdMoth[j+1];
  }
  el_n--;
}
    
void RemoveMuo(Int_t i) {
  for(Int_t j=i; j<mu_n-1; j++) {
    mu_E[j] = mu_E[j+1];
    mu_pt[j] = mu_pt[j+1];
    mu_eta[j] = mu_eta[j+1];
    mu_phi[j] = mu_phi[j+1];
    mu_E_cor[j] = mu_E_cor[j+1];
    mu_pt_cor[j] = mu_pt_cor[j+1];
    mu_eta_cor[j] = mu_eta_cor[j+1];
    mu_phi_cor[j] = mu_phi_cor[j+1];
    mu_charge[j] = mu_charge[j+1];
    mu_x[j] = mu_x[j+1];
    mu_y[j] = mu_y[j+1];
    mu_z[j] = mu_z[j+1];
    mu_xd[j] = mu_xd[j+1];
    mu_yd[j] = mu_yd[j+1];
    mu_zd[j] = mu_zd[j+1];
    mu_truth_dr[j] = mu_truth_dr[j+1];
    mu_truth_E_cor[j] = mu_truth_E_cor[j+1];
    mu_truth_pt_cor[j] = mu_truth_pt_cor[j+1];
    mu_truth_eta_cor[j] = mu_truth_eta_cor[j+1];
    mu_truth_phi_cor[j] = mu_truth_phi_cor[j+1];
    mu_truth_E[j] = mu_truth_E[j+1];
    mu_truth_pt[j] = mu_truth_pt[j+1];
    mu_truth_eta[j] = mu_truth_eta[j+1];
    mu_truth_phi[j] = mu_truth_phi[j+1];
    mu_truth_pdgIdMoth[j] = mu_truth_pdgIdMoth[j+1];
  }
  mu_n--;
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

// _p in MeV, _x,_y in mm, q is track charge (truth style)
void GetImpactPar(Double_t& xd, Double_t& yd, Double_t& zd, Double_t _px, Double_t _py, Double_t _eta, Double_t _x, Double_t _y, Double_t _z, Double_t q) {
  // to GeV and m:
  Double_t px = _px/GeV;
  Double_t py = _py/GeV;
  Double_t x = _x/1000;
  Double_t y = _y/1000;
  Double_t z = _z/1000;
  
  Double_t fBz = 3.5;
  Double_t pt = sqrt(px*px + py*py);
  
  Double_t r = pt / (q * fBz) * 1.0E9/c_light;        // in [m]

  Double_t phi_0 = TMath::ATan2(py, px); // [rad] in [-pi, pi]

  // 2. helix axis coordinates
  Double_t x_c = x + r*TMath::Sin(phi_0);
  Double_t y_c = y - r*TMath::Cos(phi_0);

  Double_t r_c = TMath::Hypot(x_c, y_c);
  //Double_t phi_c = TMath::ATan2(y_c, x_c);
  //Double_t phi = phi_c;
  //if(x_c < 0.0) phi += TMath::Pi();

  Double_t rcu = TMath::Abs(r);
  Double_t rc2 = r_c*r_c;
  // calculate coordinates of closest approach to track circle in transverse plane xd, yd, zd
  xd = x_c*x_c*x_c - x_c*rcu*r_c + x_c*y_c*y_c;
  xd = (rc2 > 0.0) ? xd / rc2 : -999;
  yd = y_c*(-rcu*r_c + rc2);
  yd = (rc2 > 0.0) ? yd / rc2 : -999;

  Double_t dphi = TMath::ATan2(y-y_c, x-x_c) - TMath::ATan2(yd-y_c, xd-x_c);
  if(dphi>=M_PI) dphi -= 2*M_PI;
  else if(dphi<-M_PI) dphi += 2*M_PI;
  Double_t arc = 0;
  if(q<0) arc = dphi*fabs(r);
  else arc = -dphi*fabs(r);
  zd = z - arc*sinh(_eta);
      
  // m to mm:
  xd *= 1000;
  yd *= 1000;
  zd *= 1000;

  //Double_t dxy = (xd*py - yd*px)/pt;
}

void AnalyzeTau(TClonesArray *branchGenParticle, GenParticle *particle, TLorentzVector& vis, TLorentzVector& invis, TLorentzVector& nu3, int& id, std::vector<TLorentzVector>& trk, std::vector<TVector3>& trk_pos, std::vector<TVector3>& trk_d, std::vector<int>& trk_charge, std::vector<int>& trk_id, std::vector<TLorentzVector>& neut, std::vector<TVector3>& neut_pos, std::vector<int>& neut_id) {
  TLorentzVector _vis;
  TLorentzVector _invis;
  TLorentzVector _nu3;
  TVector3 _pos;
  TVector3 _d;
  int _id = 0;
  Double_t _xd, _yd, _zd;

  for(Int_t i=particle->D1; i<=particle->D2; i++) {
    GenParticle *dau = (GenParticle*) branchGenParticle->At(i);
    int kidId = dau->PID;
    int kidIda = abs(dau->PID);
    int kidStat = dau->Status;
    
    if(kidId==particle->PID) {
      AnalyzeTau(branchGenParticle,dau,vis,invis,nu3,_id,trk,trk_pos,trk_d,trk_charge,trk_id,neut,neut_pos,neut_id);
    }
    else if(kidIda==12 || kidIda==14 || kidIda==16) {
      TLorentzVector tmp;
      tmp.SetPxPyPzE(dau->Px,dau->Py,dau->Pz,dau->E);
      _invis += tmp;
      if(kidIda==16) {
	tmp.SetPxPyPzE(dau->Px,dau->Py,dau->Pz,dau->E);
	_nu3 += tmp;
      }
    }
    else if(
      kidIda==11 || //e+/-
      kidIda==13 || //mu+/-
      kidIda==211 || //pi+/-
      kidIda==321 || //K+/-
      kidIda==323 || //K*+/-
      kidIda==2212) { //p+/-
      _vis.SetPxPyPzE(dau->Px,dau->Py,dau->Pz,dau->E);
      _pos.SetXYZ(dau->X,dau->Y,dau->Z);
      GetImpactPar(_xd,_yd,_zd,dau->Px*GeV,dau->Py*GeV,dau->Eta,dau->X,dau->Y,dau->Z,dau->Charge);
      _d.SetXYZ(_xd,_yd,_zd);
      trk.push_back(_vis);
      trk_pos.push_back(_pos);
      trk_d.push_back(_d);
      trk_charge.push_back(dau->Charge);
      trk_id.push_back(kidId);
      if((kidIda==11 || kidIda==13) && _id==0) _id = kidId;
    }
    else if(
      (kidId==22 && kidStat==1) || //pi0
      kidId==111 || //pi0
      kidId==130 || //K_L0
      kidId==221 || //eta
      kidId==223 || //omega
      kidId==310 || //K_S0
      kidIda==311 || //K0
      kidIda==2112) { //n
      _vis.SetPxPyPzE(dau->Px,dau->Py,dau->Pz,dau->E);
      _pos.SetXYZ(dau->X,dau->Y,dau->Z);
      neut.push_back(_vis);
      neut_pos.push_back(_pos);
      neut_id.push_back(kidId);
    }
    else if(kidStat!=1) {
      AnalyzeTau(branchGenParticle,dau,vis,invis,nu3,id,trk,trk_pos,trk_d,trk_charge,trk_id,neut,neut_pos,neut_id);
    }
    else {
      printf("AnalyzeTau: unknown tau daughter PID = %d, Status = %d\n",kidId,kidStat);
    }
  }

  TLorentzVector _all;
  _all.SetPxPyPzE(particle->Px,particle->Py,particle->Pz,particle->E);
  _invis += invis;
  _nu3 += nu3;
  _vis = _all - _invis;

  vis = _vis;
  invis = _invis;
  nu3 = _nu3;
  id = _id;
}

void findLeading(const std::vector<TLorentzVector>& vec, int* array, int dim, bool usePt=true) {
  float max = 0;
  for(int i=0; i<dim; i++) {
    max = 0;
    for(int j=0; j<(int)vec.size(); j++) {
      bool skip = false;
      for(int k=0; k<i; k++) {
	if(j==array[k]) {
	  skip = true;
	  break;
	}
      }
      if(skip) continue;
      if(usePt) {
	if(vec[j].Pt()>max) {
	  max = vec[j].Pt();
	  array[i] = j;
	}
      }
      else {
	if(vec[j].E()>max) {
	  max = vec[j].E();
	  array[i] = j;
	}
      }
    }
  }
}

void AnalyzeTauReco(Jet* tau, TClonesArray *branchTrack, TClonesArray *branchNeutralHadron, TClonesArray *branchNeutral, std::vector<TLorentzVector>& trk, std::vector<TVector3>& trk_pos, std::vector<TVector3>& trk_d, std::vector<int>& trk_charge, std::vector<int>& trk_id, std::vector<TLorentzVector>& neut) {
  TLorentzVector cand;
  TVector3 pos;
  TVector3 imp;
  TLorentzVector vtau = tau->P4();

  for(Int_t i=0; i<branchTrack->GetEntries(); i++) {
    Track *track = (Track*) branchTrack->At(i);
    cand.SetPtEtaPhiM(track->PT, track->Eta, track->Phi, 0.139570); //pi+/-
    if(vtau.DeltaR(cand)<0.2) {
      pos.SetXYZ(track->X, track->Y, track->Z);
      imp.SetXYZ(track->Xd, track->Yd, track->Zd);
      trk.push_back(cand);
      trk_pos.push_back(pos);
      trk_d.push_back(imp);
      trk_charge.push_back(track->Charge);
      trk_id.push_back(track->PID);
    }
  }

  for(Int_t i=0; i<branchNeutralHadron->GetEntries(); i++) {
    Tower *tower = (Tower*) branchNeutralHadron->At(i);
    cand.SetPtEtaPhiM(tower->ET, tower->Eta, tower->Phi, 0);
    if(vtau.DeltaR(cand)<0.2) {
      neut.push_back(cand);
    }
  }

  for(Int_t i=0; i<branchNeutral->GetEntries(); i++) {
    Tower *gam = (Tower*) branchNeutral->At(i);
    cand.SetPtEtaPhiM(gam->ET, gam->Eta, gam->Phi, 0);
    if(vtau.DeltaR(cand)<0.2) {
      neut.push_back(cand);
    }
  }
}

Double_t LogNormal(Double_t mean, Double_t sigma)
{
  Double_t a, b;

  if(mean > 0.0)
  {
    b = TMath::Sqrt(TMath::Log((1.0 + (sigma*sigma)/(mean*mean))));
    a = TMath::Log(mean) - 0.5*b*b;

    return TMath::Exp(a + b*gRandom->Gaus(0.0, 1.0));
  }
  else
  {
    return 0.0;
  }
}
    
void reduce(const char* input="../Events/run_01/tag_1_delphes_events.root", const char* out="ntuple/sig.root")
{
  //gSystem->Load("libDelphes");
  char str[200];

  TH1F* h1a[3];
  TH1F* h2a[3];
  TH1F* h1b[3];
  TH1F* h2b[3];
  
  h1a[0] = new TH1F("h1a_1","",40,-4,4);
  h2a[0] = new TH1F("h2a_1","",20,0,100);
  h1b[0] = new TH1F("h1b_1","",40,-4,4);
  h2b[0] = new TH1F("h2b_1","",20,0,100);
  
  h1a[1] = new TH1F("h1a_2","",40,-4,4);
  h2a[1] = new TH1F("h2a_2","",20,0,100);
  h1b[1] = new TH1F("h1b_2","",40,-4,4);
  h2b[1] = new TH1F("h2b_2","",20,0,100);
  
  h1a[2] = new TH1F("h1a_3","",40,-4,4);
  h2a[2] = new TH1F("h2a_3","",20,0,100);
  h1b[2] = new TH1F("h1b_3","",40,-4,4);
  h2b[2] = new TH1F("h2b_3","",20,0,100);

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
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
  TClonesArray *branchGenParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchScalarHT = treeReader->UseBranch("ScalarHT");
  TClonesArray *branchTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchNeutral = treeReader->UseBranch("EFlowPhoton");
  TClonesArray *branchNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");

  Int_t eve1 = 0;
  Int_t eve2 = 0;

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

    el_n = 0;
    for(Int_t i=0; i<branchElectron->GetEntries(); i++) {
      Electron *ele = (Electron*) branchElectron->At(i);
      if(ele->PT>15 && fabs(ele->Eta)<2.5) {
	el_E[el_n] = ele->PT*cosh(ele->Eta)*GeV;
	el_pt[el_n] = ele->PT*GeV;
	el_eta[el_n] = ele->Eta;
	el_phi[el_n] = ele->Phi;
	el_charge[el_n] = ele->Charge;
	el_x[el_n] = ele->X;
	el_y[el_n] = ele->Y;
	el_z[el_n] = ele->Z;
	el_xd[el_n] = ele->Xd;
	el_yd[el_n] = ele->Yd;
	el_zd[el_n] = ele->Zd;
	TLorentzVector clus0;
	clus0.SetPtEtaPhiM(ele->PT,ele->Eta,ele->Phi,0);
	TLorentzVector clus;
	for(Int_t j=0; j<branchGenParticle->GetEntries(); j++) {
	  GenParticle *part = (GenParticle*) branchGenParticle->At(j);
	  if(part->PID==22 && part->Status==1 && DR(part->Eta,part->Phi,clus0.Eta(),clus0.Phi())<0.1) {
	    TLorentzVector cand;
	    cand.SetPtEtaPhiM(part->PT, part->Eta, part->Phi, 0);
	    if(clus0.DeltaR(cand)<0.1) {
	      clus += cand;
	    }
	  }
	}
	Float_t sigma = sqrt(pow(clus.E()*0.01,2) + clus.E()*pow(0.15,2));
	Float_t energy = LogNormal(clus.E(), sigma);
	if(energy>0.1) {
	  clus.SetPtEtaPhiE(energy/cosh(clus.Eta()),clus.Eta(),clus.Phi(),energy);
	  clus0 += clus;
	}
	el_E_cor[el_n] = clus0.E()*GeV;
	el_pt_cor[el_n] = clus0.Pt()*GeV;
	el_eta_cor[el_n] = clus0.Eta();
	el_phi_cor[el_n] = clus0.Phi();
	el_truth_dr[el_n] = 999;
	el_truth_E_cor[el_n] = 0;
	el_truth_pt_cor[el_n] = 0;
	el_truth_eta_cor[el_n] = 0;
	el_truth_phi_cor[el_n] = 0;
	el_truth_E[el_n] = 0;
	el_truth_pt[el_n] = 0;
	el_truth_eta[el_n] = 0;
	el_truth_phi[el_n] = 0;
	el_truth_pdgIdMoth[el_n] = 0;
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
	mu_x[mu_n] = muo->X;
	mu_y[mu_n] = muo->Y;
	mu_z[mu_n] = muo->Z;
	mu_xd[mu_n] = muo->Xd;
	mu_yd[mu_n] = muo->Yd;
	mu_zd[mu_n] = muo->Zd;
	TLorentzVector clus0;
	clus0.SetPtEtaPhiM(muo->PT,muo->Eta,muo->Phi,0);
	TLorentzVector clus;
	for(Int_t j=0; j<branchGenParticle->GetEntries(); j++) {
	  GenParticle *part = (GenParticle*) branchGenParticle->At(j);
	  if(part->PID==22 && part->Status==1 && DR(part->Eta,part->Phi,clus0.Eta(),clus0.Phi())<0.1) {
	    TLorentzVector cand;
	    cand.SetPtEtaPhiM(part->PT, part->Eta, part->Phi, 0);
	    if(clus0.DeltaR(cand)<0.1) {
	      clus += cand;
	    }
	  }
	}
	Float_t sigma = sqrt(pow(clus.E()*0.01,2) + clus.E()*pow(0.15,2));
	Float_t energy = LogNormal(clus.E(), sigma);
	if(energy>0.1) {
	  clus.SetPtEtaPhiE(energy/cosh(clus.Eta()),clus.Eta(),clus.Phi(),energy);
	  clus0 += clus;
	}
	mu_E_cor[mu_n] = clus0.E()*GeV;
	mu_pt_cor[mu_n] = clus0.Pt()*GeV;
	mu_eta_cor[mu_n] = clus0.Eta();
	mu_phi_cor[mu_n] = clus0.Phi();
	mu_truth_dr[mu_n] = 999;
	mu_truth_E_cor[mu_n] = 0;
	mu_truth_pt_cor[mu_n] = 0;
	mu_truth_eta_cor[mu_n] = 0;
	mu_truth_phi_cor[mu_n] = 0;
	mu_truth_E[mu_n] = 0;
	mu_truth_pt[mu_n] = 0;
	mu_truth_eta[mu_n] = 0;
	mu_truth_phi[mu_n] = 0;
	mu_truth_pdgIdMoth[mu_n] = 0;
	mu_n++;
      }
    }

    ph_n = 0;
    for(Int_t i=0; i<branchPhoton->GetEntries(); i++) {
      Photon *pho = (Photon*) branchPhoton->At(i);
      if(pho->E>15 && fabs(pho->Eta)<2.5) {
	Bool_t OL = false;
	for(Int_t j=0; j<el_n; j++) {
	  if(DR(el_eta[j],el_phi[j],pho->Eta,pho->Phi)<0.2) {
	    OL = true;
	    break;
	  }
	}
	if(OL) continue;
      
	OL = false;
	for(Int_t j=0; j<mu_n; j++) {
	  if(DR(mu_eta[j],mu_phi[j],pho->Eta,pho->Phi)<0.2) {
	    OL = true;
	    break;
	  }
	}
	if(OL) continue;

	Bool_t isPhoton = false;
	Bool_t isPi0 = false;
	for(Int_t j=0; j<branchGenParticle->GetEntries(); j++) {
	  GenParticle *part = (GenParticle*) branchGenParticle->At(j);
	  if(part->PID==22 && part->Status==1 && DR(part->Eta,part->Phi,pho->Eta,pho->Phi)<0.2) {
	    isPhoton = true;
	  }
	  else if(part->PID==111 && part->Status==2 && DR(part->Eta,part->Phi,pho->Eta,pho->Phi)<0.2) {
	    isPi0 = true;
	    break;
	  }
	}
	if(!isPhoton || isPi0) continue;

	ph_pt[ph_n]  = pho->PT*GeV;
	ph_eta[ph_n] = pho->Eta;
	ph_phi[ph_n] = pho->Phi;
	ph_n++;
      }
    }

    mc_n = 0;
    truthTau_n = 0;

    for(Int_t i=0; i<branchGenParticle->GetEntries(); i++) {
      GenParticle *part = (GenParticle*) branchGenParticle->At(i);
      if(!part) continue;

      // tautau:
      if(abs(part->PID)==15 && isNotSelf(branchGenParticle,part)) {	
	GenParticle *dau = part;

	TLorentzVector vis;
	TLorentzVector invis;
	TLorentzVector nu3;
	int id;
	std::vector<TLorentzVector> trk;
	std::vector<TVector3> trk_pos;
	std::vector<TVector3> trk_d;
	std::vector<int> trk_charge;
	std::vector<int> trk_id;
	std::vector<TLorentzVector> neut;
	std::vector<TVector3> neut_pos;
	std::vector<int> neut_id;
	AnalyzeTau(branchGenParticle,dau,vis,invis,nu3,id,trk,trk_pos,trk_d,trk_charge,trk_id,neut,neut_pos,neut_id);
	// bool hasGam = false;
	// for(std::vector<int>::size_type j=0; j<neut_id.size(); j++){
	//   if(neut_id[j]==22) { hasGam = true; break; }
	// }
	// if(hasGam) {
	//   for(std::vector<int>::size_type j=0; j<trk_id.size(); j++) {
	// 	printf(" %d(%f %f %f)",trk_id[j],trk[j].E(),trk[j].Eta(),trk[j].Phi());
	//   }
	//   for(std::vector<int>::size_type j=0; j<neut_id.size(); j++) {
	// 	printf(" %d(%f %f %f)",neut_id[j],neut[j].E(),neut[j].Eta(),neut[j].Phi());
	//   }
	//   printf("  %d\n",id);
	// }
	int itrk[3] = {-1,-1,-1};
	findLeading(trk,itrk,3);

	truthTau_vis_pdgId[truthTau_n] = id;
	truthTau_vis_pdgIdMoth[truthTau_n] = GetMotherID(branchGenParticle,dau,true);
	truthTau_vis_pt[truthTau_n] = vis.Pt()*GeV;
	truthTau_vis_m[truthTau_n] = vis.M()*GeV;
	truthTau_vis_eta[truthTau_n] = vis.Eta();
	truthTau_vis_phi[truthTau_n] = vis.Phi();
	truthTau_vis_charge[truthTau_n] = dau->Charge;
	truthTau_vis_numTrack[truthTau_n] = (int)trk.size();
	truthTau_vis_numNeut[truthTau_n] = (int)neut.size();
	truthTau_vis_numPi0[truthTau_n] = 0;
	truthTau_vis_numGam[truthTau_n] = 0;
	for(std::vector<int>::size_type j=0; j<neut_id.size(); j++) {
	  if(neut_id[j]==111) truthTau_vis_numPi0[truthTau_n]++;
	  if(neut_id[j]==22) truthTau_vis_numGam[truthTau_n]++;
	}
	truthTau_nu3_pt[truthTau_n] = nu3.Pt()*GeV;
	truthTau_nu3_m[truthTau_n] = nu3.M()*GeV;
	truthTau_nu3_eta[truthTau_n] = nu3.Eta();
	truthTau_nu3_phi[truthTau_n] = nu3.Phi();
	truthTau_inv_pt[truthTau_n] = invis.Pt()*GeV;
	truthTau_inv_m[truthTau_n] = invis.M()*GeV;
	truthTau_inv_eta[truthTau_n] = invis.Eta();
	truthTau_inv_phi[truthTau_n] = invis.Phi();
	if(itrk[0]>=0) {
	  truthTau_trk1_pt[truthTau_n] = trk[itrk[0]].Pt()*GeV;
	  truthTau_trk1_eta[truthTau_n] = trk[itrk[0]].Eta();
	  truthTau_trk1_phi[truthTau_n] = trk[itrk[0]].Phi();
	  truthTau_trk1_m[truthTau_n] = trk[itrk[0]].M()*GeV;
	  truthTau_trk1_x[truthTau_n] = trk_pos[itrk[0]].X();
	  truthTau_trk1_y[truthTau_n] = trk_pos[itrk[0]].Y();
	  truthTau_trk1_z[truthTau_n] = trk_pos[itrk[0]].Z();
	  truthTau_trk1_xd[truthTau_n] = trk_d[itrk[0]].X();
	  truthTau_trk1_yd[truthTau_n] = trk_d[itrk[0]].Y();
	  truthTau_trk1_zd[truthTau_n] = trk_d[itrk[0]].Z();
	  truthTau_trk1_charge[truthTau_n] = trk_charge[itrk[0]];
	  truthTau_trk1_pdgId[truthTau_n] = trk_id[itrk[0]];
	}
	else {
	  truthTau_trk1_pt[truthTau_n] = 0;
	  truthTau_trk1_eta[truthTau_n] = 0;
	  truthTau_trk1_phi[truthTau_n] = 0;
	  truthTau_trk1_m[truthTau_n] = 0;
	  truthTau_trk1_x[truthTau_n] = 0;
	  truthTau_trk1_y[truthTau_n] = 0;
	  truthTau_trk1_z[truthTau_n] = 0;
	  truthTau_trk1_xd[truthTau_n] = 0;
	  truthTau_trk1_yd[truthTau_n] = 0;
	  truthTau_trk1_zd[truthTau_n] = 0;
	  truthTau_trk1_charge[truthTau_n] = 0;
	  truthTau_trk1_pdgId[truthTau_n] = 0;
	}
	if(itrk[1]>=0) {
	  truthTau_trk2_pt[truthTau_n] = trk[itrk[1]].Pt()*GeV;
	  truthTau_trk2_eta[truthTau_n] = trk[itrk[1]].Eta();
	  truthTau_trk2_phi[truthTau_n] = trk[itrk[1]].Phi();
	  truthTau_trk2_m[truthTau_n] = trk[itrk[1]].M()*GeV;
	  truthTau_trk2_x[truthTau_n] = trk_pos[itrk[1]].X();
	  truthTau_trk2_y[truthTau_n] = trk_pos[itrk[1]].Y();
	  truthTau_trk2_z[truthTau_n] = trk_pos[itrk[1]].Z();
	  truthTau_trk2_xd[truthTau_n] = trk_d[itrk[1]].X();
	  truthTau_trk2_yd[truthTau_n] = trk_d[itrk[1]].Y();
	  truthTau_trk2_zd[truthTau_n] = trk_d[itrk[1]].Z();
	  truthTau_trk2_charge[truthTau_n] = trk_charge[itrk[1]];
	  truthTau_trk2_pdgId[truthTau_n] = trk_id[itrk[1]];
	}
	else {
	  truthTau_trk2_pt[truthTau_n] = 0;
	  truthTau_trk2_eta[truthTau_n] = 0;
	  truthTau_trk2_phi[truthTau_n] = 0;
	  truthTau_trk2_m[truthTau_n] = 0;
	  truthTau_trk2_x[truthTau_n] = 0;
	  truthTau_trk2_y[truthTau_n] = 0;
	  truthTau_trk2_z[truthTau_n] = 0;
	  truthTau_trk2_xd[truthTau_n] = 0;
	  truthTau_trk2_yd[truthTau_n] = 0;
	  truthTau_trk2_zd[truthTau_n] = 0;
	  truthTau_trk2_charge[truthTau_n] = 0;
	  truthTau_trk2_pdgId[truthTau_n] = 0;
	}
	if(itrk[2]>=0) {
	  truthTau_trk3_pt[truthTau_n] = trk[itrk[2]].Pt()*GeV;
	  truthTau_trk3_eta[truthTau_n] = trk[itrk[2]].Eta();
	  truthTau_trk3_phi[truthTau_n] = trk[itrk[2]].Phi();
	  truthTau_trk3_m[truthTau_n] = trk[itrk[2]].M()*GeV;
	  truthTau_trk3_x[truthTau_n] = trk_pos[itrk[2]].X();
	  truthTau_trk3_y[truthTau_n] = trk_pos[itrk[2]].Y();
	  truthTau_trk3_z[truthTau_n] = trk_pos[itrk[2]].Z();
	  truthTau_trk3_xd[truthTau_n] = trk_d[itrk[2]].X();
	  truthTau_trk3_yd[truthTau_n] = trk_d[itrk[2]].Y();
	  truthTau_trk3_zd[truthTau_n] = trk_d[itrk[2]].Z();
	  truthTau_trk3_charge[truthTau_n] = trk_charge[itrk[2]];
	  truthTau_trk3_pdgId[truthTau_n] = trk_id[itrk[2]];
	}
	else {
	  truthTau_trk3_pt[truthTau_n] = 0;
	  truthTau_trk3_eta[truthTau_n] = 0;
	  truthTau_trk3_phi[truthTau_n] = 0;
	  truthTau_trk3_m[truthTau_n] = 0;
	  truthTau_trk3_x[truthTau_n] = 0;
	  truthTau_trk3_y[truthTau_n] = 0;
	  truthTau_trk3_z[truthTau_n] = 0;
	  truthTau_trk3_xd[truthTau_n] = 0;
	  truthTau_trk3_yd[truthTau_n] = 0;
	  truthTau_trk3_zd[truthTau_n] = 0;
	  truthTau_trk3_charge[truthTau_n] = 0;
	  truthTau_trk3_pdgId[truthTau_n] = 0;
	}

	int ineut[3] = {-1,-1,-1};
	findLeading(neut,ineut,3,false);

	if(ineut[0]>=0) {
	  truthTau_neut1_e[truthTau_n] = neut[ineut[0]].E()*GeV;
	  truthTau_neut1_pt[truthTau_n] = neut[ineut[0]].Pt()*GeV;
	  truthTau_neut1_eta[truthTau_n] = neut[ineut[0]].Eta();
	  truthTau_neut1_phi[truthTau_n] = neut[ineut[0]].Phi();
	  truthTau_neut1_pdgId[truthTau_n] = neut_id[ineut[0]];
	}
	else {
	  truthTau_neut1_e[truthTau_n] = 0;
	  truthTau_neut1_pt[truthTau_n] = 0;
	  truthTau_neut1_eta[truthTau_n] = 0;
	  truthTau_neut1_phi[truthTau_n] = 0;
	  truthTau_neut1_pdgId[truthTau_n] = 0;
	} 
	if(ineut[1]>=0) {
	  truthTau_neut2_e[truthTau_n] = neut[ineut[1]].E()*GeV;
	  truthTau_neut2_pt[truthTau_n] = neut[ineut[1]].Pt()*GeV;
	  truthTau_neut2_eta[truthTau_n] = neut[ineut[1]].Eta();
	  truthTau_neut2_phi[truthTau_n] = neut[ineut[1]].Phi();
	  truthTau_neut2_pdgId[truthTau_n] = neut_id[ineut[1]];
	}
	else {
	  truthTau_neut2_e[truthTau_n] = 0;
	  truthTau_neut2_pt[truthTau_n] = 0;
	  truthTau_neut2_eta[truthTau_n] = 0;
	  truthTau_neut2_phi[truthTau_n] = 0;
	  truthTau_neut2_pdgId[truthTau_n] = 0;
	} 
	if(ineut[2]>=0) {
	  truthTau_neut3_e[truthTau_n] = neut[ineut[2]].E()*GeV;
	  truthTau_neut3_pt[truthTau_n] = neut[ineut[2]].Pt()*GeV;
	  truthTau_neut3_eta[truthTau_n] = neut[ineut[2]].Eta();
	  truthTau_neut3_phi[truthTau_n] = neut[ineut[2]].Phi();
	  truthTau_neut3_pdgId[truthTau_n] = neut_id[ineut[2]];
	}
	else {
	  truthTau_neut3_e[truthTau_n] = 0;
	  truthTau_neut3_pt[truthTau_n] = 0;
	  truthTau_neut3_eta[truthTau_n] = 0;
	  truthTau_neut3_phi[truthTau_n] = 0;
	  truthTau_neut3_pdgId[truthTau_n] = 0;
	} 

	truthTau_n++;

	if(abs(id)==11) {
	  GenParticle* firstPart = GetFirstParticle(branchGenParticle,dau,11);
	  GenParticle* finalPart = firstPart ? GetFinalParticle(branchGenParticle,firstPart) : 0;
	  if(firstPart && finalPart) {
	    for(Int_t j=0; j<el_n; j++) {
	      if(DR(el_eta[j],el_phi[j],finalPart->Eta,finalPart->Phi)<el_truth_dr[j]) {
		el_truth_dr[j] = DR(el_eta[j],el_phi[j],finalPart->Eta,finalPart->Phi);
		el_truth_E_cor[j] = firstPart->PT*cosh(firstPart->Eta)*GeV;
		el_truth_pt_cor[j] = firstPart->PT*GeV;
		el_truth_eta_cor[j] = firstPart->Eta;
		el_truth_phi_cor[j] = firstPart->Phi;
		el_truth_E[j] = finalPart->PT*cosh(finalPart->Eta)*GeV;
		el_truth_pt[j] = finalPart->PT*GeV;
		el_truth_eta[j] = finalPart->Eta;
		el_truth_phi[j] = finalPart->Phi;
		el_truth_pdgIdMoth[j] = dau->PID;
	      }
	    }
	  }
	}
	else if(abs(id)==13) {
	  GenParticle* firstPart = GetFirstParticle(branchGenParticle,dau,13);
	  GenParticle* finalPart = firstPart ? GetFinalParticle(branchGenParticle,firstPart) : 0;
	  if(firstPart && finalPart) {
	    for(Int_t j=0; j<mu_n; j++) {
	      if(DR(mu_eta[j],mu_phi[j],finalPart->Eta,finalPart->Phi)<mu_truth_dr[j]) {
		mu_truth_dr[j] = DR(mu_eta[j],mu_phi[j],finalPart->Eta,finalPart->Phi);
		mu_truth_E_cor[j] = firstPart->PT*cosh(firstPart->Eta)*GeV;
		mu_truth_pt_cor[j] = firstPart->PT*GeV;
		mu_truth_eta_cor[j] = firstPart->Eta;
		mu_truth_phi_cor[j] = firstPart->Phi;
		mu_truth_E[j] = finalPart->PT*cosh(finalPart->Eta)*GeV;
		mu_truth_pt[j] = finalPart->PT*GeV;
		mu_truth_eta[j] = finalPart->Eta;
		mu_truth_phi[j] = finalPart->Phi;
		mu_truth_pdgIdMoth[j] = dau->PID;
	      }
	    }
	  }
	}
      }

      if(abs(part->PID)==11) {
	GenParticle* firstPart = GetFirstSelf(branchGenParticle,part);
	GenParticle* finalPart = GetFinalSelf(branchGenParticle,part);
	for(Int_t j=0; j<el_n; j++) {
	  if(DR(el_eta[j],el_phi[j],finalPart->Eta,finalPart->Phi)<el_truth_dr[j]) {
	    el_truth_dr[j] = DR(el_eta[j],el_phi[j],finalPart->Eta,finalPart->Phi);
	    el_truth_E_cor[j] = firstPart->PT*cosh(part->Eta)*GeV;
	    el_truth_pt_cor[j] = firstPart->PT*GeV;
	    el_truth_eta_cor[j] = firstPart->Eta;
	    el_truth_phi_cor[j] = firstPart->Phi;
	    el_truth_E[j] = finalPart->PT*cosh(part->Eta)*GeV;
	    el_truth_pt[j] = finalPart->PT*GeV;
	    el_truth_eta[j] = finalPart->Eta;
	    el_truth_phi[j] = finalPart->Phi;
	    el_truth_pdgIdMoth[j] = part->PID;
	  }
	}
      }
      else if(abs(part->PID)==13) {
	GenParticle* firstPart = GetFirstSelf(branchGenParticle,part);
	GenParticle* finalPart = GetFinalSelf(branchGenParticle,part);
	for(Int_t j=0; j<mu_n; j++) {
	  if(DR(mu_eta[j],mu_phi[j],finalPart->Eta,finalPart->Phi)<mu_truth_dr[j]) {
	    mu_truth_dr[j] = DR(mu_eta[j],mu_phi[j],finalPart->Eta,finalPart->Phi);
	    mu_truth_E_cor[j] = firstPart->PT*cosh(part->Eta)*GeV;
	    mu_truth_pt_cor[j] = firstPart->PT*GeV;
	    mu_truth_eta_cor[j] = firstPart->Eta;
	    mu_truth_phi_cor[j] = firstPart->Phi;
	    mu_truth_E[j] = finalPart->PT*cosh(part->Eta)*GeV;
	    mu_truth_pt[j] = finalPart->PT*GeV;
	    mu_truth_eta[j] = finalPart->Eta;
	    mu_truth_phi[j] = finalPart->Phi;
	    mu_truth_pdgIdMoth[j] = part->PID;
	  }
	}
      }
    
      if(part->PID==23 && isNotSelf(branchGenParticle,part)) {
        mc_pt[mc_n] = part->PT*GeV;
        mc_m[mc_n] = part->Mass*GeV;
        mc_eta[mc_n] = part->Eta;
        mc_phi[mc_n] = part->Phi;
        mc_status[mc_n] = part->Status;
        mc_pdgId[mc_n] = part->PID;
        mc_pdgIdMoth[mc_n] = GetMotherID(branchGenParticle,part,true);
        mc_charge[mc_n] = part->Charge;
        mc_n++;

	for(Int_t it=part->D1; it<=part->D2; it++) {
          GenParticle *dau_ = (GenParticle*) branchGenParticle->At(it);
	  GenParticle *dau = GetFinalSelf(branchGenParticle,dau_);
	  mc_pt[mc_n] = dau->PT*GeV;
	  mc_m[mc_n] = dau->Mass*GeV;
	  mc_eta[mc_n] = dau->Eta;
	  mc_phi[mc_n] = dau->Phi;
	  mc_status[mc_n] = dau->Status;
	  mc_pdgId[mc_n] = dau->PID;
	  mc_pdgIdMoth[mc_n] = part->PID;
	  mc_charge[mc_n] = dau->Charge;
	  mc_n++;
	}
      }
    }
    //if(entry<10) printf("\n");

    tau_n = 0;
    for(Int_t i=0; i<branchJet->GetEntries(); i++) {
      Jet *jet = (Jet*) branchJet->At(i);
      
      Int_t iEl = -1;
      for(Int_t j=0; j<el_n; j++) {
	if(DR(el_eta[j],el_phi[j],jet->Eta,jet->Phi)<0.2) {
	  iEl = j;
	}
      }

      Int_t iMu = -1;
      for(Int_t j=0; j<mu_n; j++) {
	if(DR(mu_eta[j],mu_phi[j],jet->Eta,jet->Phi)<0.2) {
	  iMu = j;
	}
      }

      // OL = false;
      // for(Int_t j=0; j<ph_n; j++) {
      // 	if(DR(ph_eta[j],ph_phi[j],jet->Eta,jet->Phi)<0.2) {
      // 	  OL = true;
      // 	  break;
      // 	}
      // }
      // if(OL) continue;

      if(jet->PT>25 && fabs(jet->Eta)<2.5 && jet->TauTag==1) {
	tau_E[tau_n] = sqrt(pow(jet->PT*cosh(jet->Eta),2)+pow(jet->Mass,2))*GeV;
        tau_pt[tau_n] = jet->PT*GeV;
        tau_eta[tau_n] = jet->Eta;
        tau_phi[tau_n] = jet->Phi;

	std::vector<TLorentzVector> trk;
	std::vector<TVector3> trk_pos;
	std::vector<TVector3> trk_d;
	std::vector<int> trk_charge;
	std::vector<int> trk_id;
	std::vector<TLorentzVector> neut;
	AnalyzeTauReco(jet,branchTrack,branchNeutralHadron,branchNeutral,trk,trk_pos,trk_d,trk_charge,trk_id,neut);
	int itrk[3] = {-1,-1,-1};
	findLeading(trk,itrk,3);

	tau_numTrack[tau_n] = (int)trk.size();
	tau_numNeut[tau_n] = (int)neut.size();
        tau_charge[tau_n] = 0;
	
	if(itrk[0]>=0) {
	  tau_trk1_pt[tau_n] = trk[itrk[0]].Pt()*GeV;
	  tau_trk1_eta[tau_n] = trk[itrk[0]].Eta();
	  tau_trk1_phi[tau_n] = trk[itrk[0]].Phi();
	  tau_trk1_m[tau_n] = trk[itrk[0]].M()*GeV;
	  tau_trk1_x[tau_n] = trk_pos[itrk[0]].X();
	  tau_trk1_y[tau_n] = trk_pos[itrk[0]].Y();
	  tau_trk1_z[tau_n] = trk_pos[itrk[0]].Z();
	  tau_trk1_xd[tau_n] = trk_d[itrk[0]].X();
	  tau_trk1_yd[tau_n] = trk_d[itrk[0]].Y();
	  tau_trk1_zd[tau_n] = trk_d[itrk[0]].Z();
	  tau_trk1_charge[tau_n] = trk_charge[itrk[0]];
	  tau_trk1_pdgId[tau_n] = trk_id[itrk[0]];
	  tau_charge[tau_n] += trk_charge[itrk[0]];
	}
	else {
	  tau_trk1_pt[tau_n] = 0;
	  tau_trk1_eta[tau_n] = 0;
	  tau_trk1_phi[tau_n] = 0;
	  tau_trk1_m[tau_n] = 0;
	  tau_trk1_x[tau_n] = 0;
	  tau_trk1_y[tau_n] = 0;
	  tau_trk1_z[tau_n] = 0;
	  tau_trk1_xd[tau_n] = 0;
	  tau_trk1_yd[tau_n] = 0;
	  tau_trk1_zd[tau_n] = 0;
	  tau_trk1_charge[tau_n] = 0;
	  tau_trk1_pdgId[tau_n] = 0;
	} 
	
	if(itrk[1]>=0) {
	  tau_trk2_pt[tau_n] = trk[itrk[1]].Pt()*GeV;
	  tau_trk2_eta[tau_n] = trk[itrk[1]].Eta();
	  tau_trk2_phi[tau_n] = trk[itrk[1]].Phi();
	  tau_trk2_m[tau_n] = trk[itrk[1]].M()*GeV;
	  tau_trk2_x[tau_n] = trk_pos[itrk[1]].X();
	  tau_trk2_y[tau_n] = trk_pos[itrk[1]].Y();
	  tau_trk2_z[tau_n] = trk_pos[itrk[1]].Z();
	  tau_trk2_xd[tau_n] = trk_d[itrk[1]].X();
	  tau_trk2_yd[tau_n] = trk_d[itrk[1]].Y();
	  tau_trk2_zd[tau_n] = trk_d[itrk[1]].Z();
	  tau_trk2_charge[tau_n] = trk_charge[itrk[1]];
	  tau_trk2_pdgId[tau_n] = trk_id[itrk[1]];
	  tau_charge[tau_n] += trk_charge[itrk[1]];
	}
	else {
	  tau_trk2_pt[tau_n] = 0;
	  tau_trk2_eta[tau_n] = 0;
	  tau_trk2_phi[tau_n] = 0;
	  tau_trk2_m[tau_n] = 0;
	  tau_trk2_x[tau_n] = 0;
	  tau_trk2_y[tau_n] = 0;
	  tau_trk2_z[tau_n] = 0;
	  tau_trk2_xd[tau_n] = 0;
	  tau_trk2_yd[tau_n] = 0;
	  tau_trk2_zd[tau_n] = 0;
	  tau_trk2_charge[tau_n] = 0;
	  tau_trk2_pdgId[tau_n] = 0;
	} 

	if(itrk[2]>=0) {
	  tau_trk3_pt[tau_n] = trk[itrk[2]].Pt()*GeV;
	  tau_trk3_eta[tau_n] = trk[itrk[2]].Eta();
	  tau_trk3_phi[tau_n] = trk[itrk[2]].Phi();
	  tau_trk3_m[tau_n] = trk[itrk[2]].M()*GeV;
	  tau_trk3_x[tau_n] = trk_pos[itrk[2]].X();
	  tau_trk3_y[tau_n] = trk_pos[itrk[2]].Y();
	  tau_trk3_z[tau_n] = trk_pos[itrk[2]].Z();
	  tau_trk3_xd[tau_n] = trk_d[itrk[2]].X();
	  tau_trk3_yd[tau_n] = trk_d[itrk[2]].Y();
	  tau_trk3_zd[tau_n] = trk_d[itrk[2]].Z();
	  tau_trk3_charge[tau_n] = trk_charge[itrk[2]];
	  tau_trk3_pdgId[tau_n] = trk_id[itrk[2]];
	  tau_charge[tau_n] += trk_charge[itrk[2]];
	}
	else {
	  tau_trk3_pt[tau_n] = 0;
	  tau_trk3_eta[tau_n] = 0;
	  tau_trk3_phi[tau_n] = 0;
	  tau_trk3_m[tau_n] = 0;
	  tau_trk3_x[tau_n] = 0;
	  tau_trk3_y[tau_n] = 0;
	  tau_trk3_z[tau_n] = 0;
	  tau_trk3_xd[tau_n] = 0;
	  tau_trk3_yd[tau_n] = 0;
	  tau_trk3_zd[tau_n] = 0;
	  tau_trk3_charge[tau_n] = 0;
	  tau_trk3_pdgId[tau_n] = 0;
	}

        tau_n++;
      }
    }

    jet_n = 0;
    for(Int_t i=0; i<branchJet->GetEntries(); i++) {
      Jet *jet = (Jet*) branchJet->At(i);
      
      Int_t iEl = -1;
      for(Int_t j=0; j<el_n; j++) {
	if(DR(el_eta[j],el_phi[j],jet->Eta,jet->Phi)<0.2) {
	  iEl = j;
	}
      }
      if(iEl>=0) continue;

      Int_t iMu = -1;
      for(Int_t j=0; j<mu_n; j++) {
	if(DR(mu_eta[j],mu_phi[j],jet->Eta,jet->Phi)<0.2) {
	  iMu = j;
	}
      }
      if(iMu>=0) continue;

      Int_t iTau = -1;
      for(Int_t j=0; j<tau_n; j++) {
	if(DR(tau_eta[j],tau_phi[j],jet->Eta,jet->Phi)<0.2) {
	  iTau = j;
	}
      }
      if(iTau>=0) continue;

      if((jet->PT>25 && fabs(jet->Eta)<2.5) || (jet->PT>30 && fabs(jet->Eta)<4.5)) {
	jet_E[jet_n] = sqrt(pow(jet->PT*cosh(jet->Eta),2)+pow(jet->Mass,2))*GeV;
	jet_m[jet_n] = jet->Mass*GeV;
	jet_pt[jet_n] = jet->PT*GeV;
	jet_eta[jet_n] = jet->Eta;
	jet_phi[jet_n] = jet->Phi;
	jet_isBtagged[jet_n] = fabs(jet->Eta)<2.5 ? (jet->BTag % 2)==1 : false;
	jet_n++;
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

    MissingET *m_MET = (MissingET*) branchMissingET->At(0);
    MET_ex = m_MET->MET*GeV*cos(m_MET->Phi);
    MET_ey = m_MET->MET*GeV*sin(m_MET->Phi);
    MET_et = m_MET->MET*GeV;
    MET_phi = m_MET->Phi;

    ScalarHT *m_HT = (ScalarHT*) branchScalarHT->At(0);
    MET_sumet = m_HT->HT*GeV;

    MET_Truth_NonInt_ex = 0;
    MET_Truth_NonInt_ey = 0;
    MET_Truth_NonInt_ez = 0;

    for(Int_t i=0; i<truthTau_n; i++) {
      MET_Truth_NonInt_ex += truthTau_inv_pt[i]*cos(truthTau_inv_phi[i]);
      MET_Truth_NonInt_ey += truthTau_inv_pt[i]*sin(truthTau_inv_phi[i]);
      MET_Truth_NonInt_ez += truthTau_inv_pt[i]*sinh(truthTau_inv_eta[i]);
    }
    for(Int_t i=0; i<mc_n; i++) {
      if((abs(mc_pdgId[i])==12 || abs(mc_pdgId[i])==14 || abs(mc_pdgId[i])==16) && (abs(mc_pdgIdMoth[i])==24 || abs(mc_pdgIdMoth[i])==37 || mc_pdgIdMoth[i]==23)) {
	MET_Truth_NonInt_ex += mc_pt[i]*cos(mc_phi[i]);
	MET_Truth_NonInt_ey += mc_pt[i]*sin(mc_phi[i]);
	MET_Truth_NonInt_ez += mc_pt[i]*sinh(mc_eta[i]);
      }
    }
    MET_Truth_NonInt_et = sqrt(pow(MET_Truth_NonInt_ex,2)+pow(MET_Truth_NonInt_ey,2));

    // eff study
    for(Int_t i=0; i<mc_n; i++) {
      if(abs(mc_pdgId[i])==11 && mc_pdgIdMoth[i]==23) {
	h1a[0]->Fill(mc_eta[i]);
	h2a[0]->Fill(mc_pt[i]/GeV);

	Int_t idx(-1);
	for(Int_t j=0; j<el_n; j++) {
	  if(DR(mc_eta[i],mc_phi[i],el_eta[j],el_phi[j])<0.2) {
	    idx = j; break;
	  }
	}
	if(idx>=0) {
	  h1b[0]->Fill(el_eta[idx]);
	  h2b[0]->Fill(el_pt[idx]/GeV);
	}
      }
      else if(abs(mc_pdgId[i])==13 && mc_pdgIdMoth[i]==23) {
	h1a[1]->Fill(mc_eta[i]);
	h2a[1]->Fill(mc_pt[i]/GeV);
	
	Int_t idx(-1);
	for(Int_t j=0; j<mu_n; j++) {
	  if(DR(mc_eta[i],mc_phi[i],mu_eta[j],mu_phi[j])<0.2) {
	    idx = j; break;
	  }
	}
	if(idx>=0) {
	  h1b[1]->Fill(mu_eta[idx]);
	  h2b[1]->Fill(mu_pt[idx]/GeV);
	}
      }
    }
    
    for(Int_t i=0; i<truthTau_n; i++) {
      if(abs(truthTau_vis_pdgId[i])==0) {
	h1a[2]->Fill(truthTau_vis_eta[i]);
	//h2a[2]->Fill(truthTau_vis_pt[i]/GeV);
	h2a[2]->Fill(truthTau_vis_pt[i]*cosh(truthTau_vis_eta[i])/GeV);

	Int_t idx(-1);
	for(Int_t j=0; j<tau_n; j++) {
	  if(DR(truthTau_vis_eta[i],truthTau_vis_phi[i],tau_eta[j],tau_phi[j])<0.2) {
	    idx = j; break;
	  }
	}
	if(idx>=0) {
	  h1b[2]->Fill(tau_eta[idx]);
	  //h2b[2]->Fill(tau_pt[idx]/GeV);
	  h2b[2]->Fill(tau_pt[idx]*cosh(tau_eta[idx])/GeV);
	}
      }
    }

    if(tau_n>=2 && el_n+mu_n>=2) {
    //if(el_n+mu_n>=4) {
      eve2++;
      tree->Fill();
    }
  }
  printf(" Tot %d, Sel %d\n",eve1,eve2);

  TCanvas* c1a = new TCanvas("c1a","c1a",600,600);
  h1a[0]->SetXTitle("#eta_{e}");
  h1a[0]->SetMinimum(0);
  h1a[0]->Draw();
  h1b[0]->SetLineColor(kRed);
  h1b[0]->Draw("same");

  TCanvas* c1b = new TCanvas("c1b","c1b",600,600);
  h2a[0]->SetXTitle("p_{T,e} [GeV]");
  h2a[0]->SetMinimum(0);
  h2a[0]->Draw();
  h2b[0]->SetLineColor(kRed);
  h2b[0]->Draw("same");

  TCanvas* c2a = new TCanvas("c2a","c2a",600,600);
  h1a[1]->SetXTitle("#eta_{#mu}");
  h1a[1]->SetMinimum(0);
  h1a[1]->Draw();
  h1b[1]->SetLineColor(kRed);
  h1b[1]->Draw("same");

  TCanvas* c2b = new TCanvas("c2b","c2b",600,600);
  h2a[1]->SetXTitle("p_{T,#mu} [GeV]");
  h2a[1]->SetMinimum(0);
  h2a[1]->Draw();
  h2b[1]->SetLineColor(kRed);
  h2b[1]->Draw("same");

  TCanvas* c3a = new TCanvas("c3a","c3a",600,600);
  h1a[2]->SetXTitle("#eta_{#tau}");
  h1a[2]->SetMinimum(0);
  h1a[2]->Draw();
  h1b[2]->SetLineColor(kRed);
  h1b[2]->Draw("same");

  TCanvas* c3b = new TCanvas("c3b","c3b",600,600);
  //h2a[2]->SetXTitle("p_{T,#tau} [GeV]");
  h2a[2]->SetXTitle("p_{#tau} [GeV]");
  h2a[2]->SetMinimum(0);
  h2a[2]->Draw();
  h2b[2]->SetLineColor(kRed);
  h2b[2]->Draw("same");
  
  fout->cd();
  tree->AutoSave();
}
