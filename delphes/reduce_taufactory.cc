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
#include "thrust.h"

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
Int_t           trk_n;
Float_t         trk_E[10];   //[trk_n]
Float_t         trk_pt[10];   //[trk_n]
Float_t         trk_eta[10];   //[trk_n]
Float_t         trk_phi[10];   //[trk_n]
Int_t           trk_charge[10];   //[trk_n]
Float_t         trk_x[10];   //[trk_n]
Float_t         trk_y[10];   //[trk_n]
Float_t         trk_z[10];   //[trk_n]
Float_t         trk_xd[10];   //[trk_n]
Float_t         trk_yd[10];   //[trk_n]
Float_t         trk_zd[10];   //[trk_n]
Int_t           pi0_n;
Float_t         pi0_E[20];   //[pi0_n]
Float_t         pi0_pt[20];   //[pi0_n]
Float_t         pi0_eta[20];   //[pi0_n]
Float_t         pi0_phi[20];   //[pi0_n]
Int_t           gam_n;
Float_t         gam_E[100];   //[gam_n]
Float_t         gam_pt[100];   //[gam_n]
Float_t         gam_eta[100];   //[gam_n]
Float_t         gam_phi[100];   //[gam_n]
Int_t           tau_n;
Float_t         tau_vis_E[10];   //[tau_n]
Float_t         tau_vis_pt[10];   //[tau_n]
Float_t         tau_vis_eta[10];   //[tau_n]
Float_t         tau_vis_phi[10];   //[tau_n]
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
Float_t         tau_neut1_pt[10];   //[tau_n]
Float_t         tau_neut1_eta[10];   //[tau_n]
Float_t         tau_neut1_phi[10];   //[tau_n]
Float_t         tau_neut1_m[10];   //[tau_n]
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
  tree->Branch("thrust", &eVal1,"thrust/D");
  tree->Branch("major", &eVal2,"major/D");
  tree->Branch("minor", &eVal3,"minor/D");
  //tree->Branch("mc_event_ScalePDF", &mc_event_ScalePDF,"mc_event_ScalePDF/F");
  //tree->Branch("mc_event_AlphaQED", &mc_event_AlphaQED,"mc_event_AlphaQED/F");
  //tree->Branch("mc_event_AlphaQCD", &mc_event_AlphaQCD,"mc_event_AlphaQCD/F");
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
  tree->Branch("pi0_n", &pi0_n,"pi0_n/I");
  tree->Branch("pi0_E", pi0_E,"pi0_E[pi0_n]/F");
  tree->Branch("pi0_pt", pi0_pt,"pi0_pt[pi0_n]/F");
  tree->Branch("pi0_eta", pi0_eta,"pi0_eta[pi0_n]/F");
  tree->Branch("pi0_phi", pi0_phi,"pi0_phi[pi0_n]/F");
  tree->Branch("gam_n", &gam_n,"gam_n/I");
  tree->Branch("gam_E", gam_E,"gam_E[gam_n]/F");
  tree->Branch("gam_pt", gam_pt,"gam_pt[gam_n]/F");
  tree->Branch("gam_eta", gam_eta,"gam_eta[gam_n]/F");
  tree->Branch("gam_phi", gam_phi,"gam_phi[gam_n]/F");
  tree->Branch("tau_n", &tau_n,"tau_n/I");
  tree->Branch("tau_vis_E", tau_vis_E,"tau_vis_E[tau_n]/F");
  tree->Branch("tau_vis_pt", tau_vis_pt,"tau_vis_pt[tau_n]/F");
  tree->Branch("tau_vis_eta", tau_vis_eta,"tau_vis_eta[tau_n]/F");
  tree->Branch("tau_vis_phi", tau_vis_phi,"tau_vis_phi[tau_n]/F");
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
  tree->Branch("tau_neut1_pt", tau_neut1_pt,"tau_neut1_pt[tau_n]/F");
  tree->Branch("tau_neut1_eta", tau_neut1_eta,"tau_neut1_eta[tau_n]/F");
  tree->Branch("tau_neut1_phi", tau_neut1_phi,"tau_neut1_phi[tau_n]/F");
  tree->Branch("tau_neut1_m", tau_neut1_m,"tau_neut1_m[tau_n]/F");
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

void SortTauTrack(Int_t i) {
  if(tau_trk1_pt[i]<tau_trk2_pt[i]) {
    Float_t tau_trk_pt = tau_trk2_pt[i];
    Float_t tau_trk_eta = tau_trk2_eta[i];
    Float_t tau_trk_phi = tau_trk2_phi[i];
    Float_t tau_trk_m = tau_trk2_m[i];
    Float_t tau_trk_x = tau_trk2_x[i];
    Float_t tau_trk_y = tau_trk2_y[i];
    Float_t tau_trk_z = tau_trk2_z[i];
    Float_t tau_trk_xd = tau_trk2_xd[i];
    Float_t tau_trk_yd = tau_trk2_yd[i];
    Float_t tau_trk_zd = tau_trk2_zd[i];
    Int_t   tau_trk_charge = tau_trk2_charge[i];
    Int_t   tau_trk_pdgId = tau_trk2_pdgId[i];
   
    tau_trk2_pt[i] = tau_trk1_pt[i];
    tau_trk2_eta[i] = tau_trk1_eta[i];
    tau_trk2_phi[i] = tau_trk1_phi[i];
    tau_trk2_m[i] = tau_trk1_m[i];
    tau_trk2_x[i] = tau_trk1_x[i];
    tau_trk2_y[i] = tau_trk1_y[i];
    tau_trk2_z[i] = tau_trk1_z[i];
    tau_trk2_xd[i] = tau_trk1_xd[i];
    tau_trk2_yd[i] = tau_trk1_yd[i];
    tau_trk2_zd[i] = tau_trk1_zd[i];
    tau_trk2_charge[i] = tau_trk1_charge[i];
    tau_trk2_pdgId[i] = tau_trk1_pdgId[i];

    tau_trk1_pt[i] = tau_trk_pt;
    tau_trk1_eta[i] = tau_trk_eta;
    tau_trk1_phi[i] = tau_trk_phi;
    tau_trk1_m[i] = tau_trk_m;
    tau_trk1_x[i] = tau_trk_x;
    tau_trk1_y[i] = tau_trk_y;
    tau_trk1_z[i] = tau_trk_z;
    tau_trk1_xd[i] = tau_trk_xd;
    tau_trk1_yd[i] = tau_trk_yd;
    tau_trk1_zd[i] = tau_trk_zd;
    tau_trk1_charge[i] = tau_trk_charge;
    tau_trk1_pdgId[i] = tau_trk_pdgId;
  }

  if(tau_trk1_pt[i]<tau_trk3_pt[i]) {
    Float_t tau_trk_pt = tau_trk3_pt[i];
    Float_t tau_trk_eta = tau_trk3_eta[i];
    Float_t tau_trk_phi = tau_trk3_phi[i];
    Float_t tau_trk_m = tau_trk3_m[i];
    Float_t tau_trk_x = tau_trk3_x[i];
    Float_t tau_trk_y = tau_trk3_y[i];
    Float_t tau_trk_z = tau_trk3_z[i];
    Float_t tau_trk_xd = tau_trk3_xd[i];
    Float_t tau_trk_yd = tau_trk3_yd[i];
    Float_t tau_trk_zd = tau_trk3_zd[i];
    Int_t   tau_trk_charge = tau_trk3_charge[i];
    Int_t   tau_trk_pdgId = tau_trk3_pdgId[i];
   
    tau_trk3_pt[i] = tau_trk1_pt[i];
    tau_trk3_eta[i] = tau_trk1_eta[i];
    tau_trk3_phi[i] = tau_trk1_phi[i];
    tau_trk3_m[i] = tau_trk1_m[i];
    tau_trk3_x[i] = tau_trk1_x[i];
    tau_trk3_y[i] = tau_trk1_y[i];
    tau_trk3_z[i] = tau_trk1_z[i];
    tau_trk3_xd[i] = tau_trk1_xd[i];
    tau_trk3_yd[i] = tau_trk1_yd[i];
    tau_trk3_zd[i] = tau_trk1_zd[i];
    tau_trk3_charge[i] = tau_trk1_charge[i];
    tau_trk3_pdgId[i] = tau_trk1_pdgId[i];

    tau_trk1_pt[i] = tau_trk_pt;
    tau_trk1_eta[i] = tau_trk_eta;
    tau_trk1_phi[i] = tau_trk_phi;
    tau_trk1_m[i] = tau_trk_m;
    tau_trk1_x[i] = tau_trk_x;
    tau_trk1_y[i] = tau_trk_y;
    tau_trk1_z[i] = tau_trk_z;
    tau_trk1_xd[i] = tau_trk_xd;
    tau_trk1_yd[i] = tau_trk_yd;
    tau_trk1_zd[i] = tau_trk_zd;
    tau_trk1_charge[i] = tau_trk_charge;
    tau_trk1_pdgId[i] = tau_trk_pdgId;
  }
  
  if(tau_trk2_pt[i]<tau_trk3_pt[i]) {
    Float_t tau_trk_pt = tau_trk3_pt[i];
    Float_t tau_trk_eta = tau_trk3_eta[i];
    Float_t tau_trk_phi = tau_trk3_phi[i];
    Float_t tau_trk_m = tau_trk3_m[i];
    Float_t tau_trk_x = tau_trk3_x[i];
    Float_t tau_trk_y = tau_trk3_y[i];
    Float_t tau_trk_z = tau_trk3_z[i];
    Float_t tau_trk_xd = tau_trk3_xd[i];
    Float_t tau_trk_yd = tau_trk3_yd[i];
    Float_t tau_trk_zd = tau_trk3_zd[i];
    Int_t   tau_trk_charge = tau_trk3_charge[i];
    Int_t   tau_trk_pdgId = tau_trk3_pdgId[i];
   
    tau_trk3_pt[i] = tau_trk2_pt[i];
    tau_trk3_eta[i] = tau_trk2_eta[i];
    tau_trk3_phi[i] = tau_trk2_phi[i];
    tau_trk3_m[i] = tau_trk2_m[i];
    tau_trk3_x[i] = tau_trk2_x[i];
    tau_trk3_y[i] = tau_trk2_y[i];
    tau_trk3_z[i] = tau_trk2_z[i];
    tau_trk3_xd[i] = tau_trk2_xd[i];
    tau_trk3_yd[i] = tau_trk2_yd[i];
    tau_trk3_zd[i] = tau_trk2_zd[i];
    tau_trk3_charge[i] = tau_trk2_charge[i];
    tau_trk3_pdgId[i] = tau_trk2_pdgId[i];

    tau_trk2_pt[i] = tau_trk_pt;
    tau_trk2_eta[i] = tau_trk_eta;
    tau_trk2_phi[i] = tau_trk_phi;
    tau_trk2_m[i] = tau_trk_m;
    tau_trk2_x[i] = tau_trk_x;
    tau_trk2_y[i] = tau_trk_y;
    tau_trk2_z[i] = tau_trk_z;
    tau_trk2_xd[i] = tau_trk_xd;
    tau_trk2_yd[i] = tau_trk_yd;
    tau_trk2_zd[i] = tau_trk_zd;
    tau_trk2_charge[i] = tau_trk_charge;
    tau_trk2_pdgId[i] = tau_trk_pdgId;
  }
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

Int_t GetMotherIdx(TClonesArray *branchGenParticle, GenParticle *particle, bool removeSelf=false) {
  if(particle->M1>=0 && particle->M2>=0) {
    for(Int_t i=particle->M1; i<=particle->M2; i++) {
      GenParticle *mon = (GenParticle*) branchGenParticle->At(i);
      if(mon->PID==particle->PID && removeSelf) return GetMotherIdx(branchGenParticle,mon,removeSelf);
    }
    return particle->M1;
  }
  else if(particle->M1>=0) {
    GenParticle *part = (GenParticle*) branchGenParticle->At(particle->M1);
    if(part->PID==particle->PID && removeSelf) return GetMotherIdx(branchGenParticle,part,removeSelf);
    return particle->M1;
  }
  return -1;
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

  TH1F* h1 = new TH1F("h1","",60,0.1,0.16);
  TH1F* h2 = new TH1F("h2","",50,0.9,1.1);
  TH1F* h3 = new TH1F("h3","",4,0,4);
  TH1F* h4 = new TH1F("h4","",100,10,11);

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
  TClonesArray *branchGenParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchNeutral = treeReader->UseBranch("EFlowPhoton");

  Float_t n_eve[2][8] = {
    {0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0}
  };

  TFile* fout = new TFile(out,"RECREATE"); //output file
  TTree* tree = new TTree("tau",str); //output tree
  SetBranch(tree);

  TRandom3 rd;

  // Loop over all events
  printf(" %d entries to be processed\n",numberOfEntries);
  for(Int_t entry = 0; entry < numberOfEntries; entry++)
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

    n_eve[0][0]++;
    n_eve[1][0] += mc_event_weight;

    el_n = 0; mu_n=0; trk_n=0;
    Int_t nK = 0; Int_t nP = 0;
    for(Int_t i=0; i<branchTrack->GetEntries(); i++) {
      Track *track = (Track*) branchTrack->At(i);

      if(abs(track->PID)==11) { // electron
	if(track->PT>0.2 && track->Eta>-1.317 && track->Eta<1.901) {
	  el_E[el_n] = track->PT*cosh(track->Eta)*GeV;
	  el_pt[el_n] = track->PT*GeV;
	  el_eta[el_n] = track->Eta;
	  el_phi[el_n] = track->Phi;
	  el_charge[el_n] = track->Charge;
	  el_x[el_n] = track->X;
	  el_y[el_n] = track->Y;
	  el_z[el_n] = track->Z;
	  el_xd[el_n] = track->Xd;
	  el_yd[el_n] = track->Yd;
	  el_zd[el_n] = track->Zd;
	  // brem recovery:
	  TLorentzVector clus0;
	  clus0.SetPtEtaPhiM(track->PT,track->Eta,track->Phi,0);
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
	  Float_t sigma = sqrt(pow(clus.E()*0.012,2) + pow(clus.E(),1.5)*pow(0.016,2) + pow(0.002,2));
	  Float_t energy = LogNormal(clus.E(), sigma);
	  if(energy>0.0) {
	    clus.SetPtEtaPhiE(energy/cosh(clus.Eta()),clus.Eta(),clus.Phi(),energy);
	    clus0 += clus;
	  }
	  el_E_cor[el_n] = clus0.E()*GeV;
	  el_pt_cor[el_n] = clus0.Pt()*GeV;
	  el_eta_cor[el_n] = clus0.Eta();
	  el_phi_cor[el_n] = clus0.Phi();
	  el_n++;
	}
      }
      else if(abs(track->PID)==13) { // muon
	if(track->PT>0.2 && track->Eta>-1.317 && track->Eta<1.901) {
	  mu_E[mu_n] = sqrt(pow(track->PT*cosh(track->Eta),2)+pow(0.105658,2))*GeV;
	  mu_pt[mu_n] = track->PT*GeV;
	  mu_eta[mu_n] = track->Eta;
	  mu_phi[mu_n] = track->Phi;
	  mu_charge[mu_n] = track->Charge;
	  mu_x[mu_n] = track->X;
	  mu_y[mu_n] = track->Y;
	  mu_z[mu_n] = track->Z;
	  mu_xd[mu_n] = track->Xd;
	  mu_yd[mu_n] = track->Yd;
	  mu_zd[mu_n] = track->Zd;
	  // brem recovery
	  TLorentzVector clus0;
	  clus0.SetPtEtaPhiM(track->PT,track->Eta,track->Phi,0);
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
	  Float_t sigma = sqrt(pow(clus.E()*0.012,2) + pow(clus.E(),1.5)*pow(0.016,2) + pow(0.002,2));
	  Float_t energy = LogNormal(clus.E(), sigma);
	  if(energy>0.0) {
	    clus.SetPtEtaPhiE(energy/cosh(clus.Eta()),clus.Eta(),clus.Phi(),energy);
	    clus0 += clus;
	  }
	  mu_E_cor[mu_n] = clus0.E()*GeV;
	  mu_pt_cor[mu_n] = clus0.Pt()*GeV;
	  mu_eta_cor[mu_n] = clus0.Eta();
	  mu_phi_cor[mu_n] = clus0.Phi();
	  mu_n++;
	}
      }
      else if(abs(track->PID)==211) { // pion
	if(track->PT>0.2 && track->Eta>-1.317 && track->Eta<1.901) {
	  trk_E[trk_n] = sqrt(pow(track->PT*cosh(track->Eta),2)+pow(0.139571,2))*GeV;
	  trk_pt[trk_n] = track->PT*GeV;
	  trk_eta[trk_n] = track->Eta;
	  trk_phi[trk_n] = track->Phi;
	  trk_charge[trk_n] = track->Charge;
	  trk_x[trk_n] = track->X;
	  trk_y[trk_n] = track->Y;
	  trk_z[trk_n] = track->Z;
	  trk_xd[trk_n] = track->Xd;
	  trk_yd[trk_n] = track->Yd;
	  trk_zd[trk_n] = track->Zd;
	  trk_n++;
	}
      }
      else if(abs(track->PID)==321) { // kaon
	if(track->PT>0.2 && track->Eta>-1.317 && track->Eta<1.901) nK++;
      }
      else if(abs(track->PID)==2212) { // proton
	if(track->PT>0.2 && track->Eta>-1.317 && track->Eta<1.901) nP++;
      }
      else { // unidentified tracks
	if(track->PT>0.2 && track->Eta>-1.317 && track->Eta<1.901) {
	  trk_E[trk_n] = sqrt(pow(track->PT*cosh(track->Eta),2)+pow(0.139571,2))*GeV;
	  trk_pt[trk_n] = track->PT*GeV;
	  trk_eta[trk_n] = track->Eta;
	  trk_phi[trk_n] = track->Phi;
	  trk_charge[trk_n] = track->Charge;
	  trk_x[trk_n] = track->X;
	  trk_y[trk_n] = track->Y;
	  trk_z[trk_n] = track->Z;
	  trk_xd[trk_n] = track->Xd;
	  trk_yd[trk_n] = track->Yd;
	  trk_zd[trk_n] = track->Zd;
	  trk_n++;
	}
      }
    }

    Int_t nKL = 0;
    pi0_n = 0; gam_n = 0;

    std::vector<bool> used_pi0;
    for(Int_t i=0; i<branchGenParticle->GetEntries(); i++) used_pi0.push_back(false);
    std::vector<bool> used_neut;
    for(Int_t i=0; i<branchNeutral->GetEntries(); i++) used_neut.push_back(false);
    for(Int_t i=0; i<branchGenParticle->GetEntries(); i++) {
      GenParticle *part = (GenParticle*) branchGenParticle->At(i);
      if(part->PID==111 && part->Status==2) {
	TLorentzVector gam1; TLorentzVector gam2;
	if(part->D1>=0 && part->D2>=0) {
	  for(Int_t j=part->D1; j<=part->D2; j++) {
	    GenParticle *dau = (GenParticle*) branchGenParticle->At(j);
	    if(dau->PID==22) {
	      if(gam1.E()==0) gam1.SetPtEtaPhiM(dau->PT,dau->Eta,dau->Phi,0);
 	      else if(gam2.E()==0) gam2.SetPtEtaPhiM(dau->PT,dau->Eta,dau->Phi,0);
	    }
	  }
	}
	if(gam1.Eta()>-1.317 && gam1.Eta()<1.901 && gam1.E()>0.5/GeV && gam2.Eta()>-1.317 && gam2.Eta()<1.901 && gam2.E()>0.5/GeV) {
	  TLorentzVector cand; Int_t nn = 0;
	  TLorentzVector neut1; TLorentzVector neut2;
	  std::vector<Int_t> used_idx;
	  for(Int_t j=0; j<branchNeutral->GetEntries(); j++) {
	    Tower *tower = (Tower*) branchNeutral->At(j);
	    if(tower->E<0.5/GeV) continue;
	    //h3->Fill(tower->E);
	    TLorentzVector neut;
	    neut.SetPtEtaPhiM(tower->ET, tower->Eta, tower->Phi, 0);
	    if(neut.DeltaR(gam1)<0.1) {
	      if(neut1.E()==0) neut1 = neut;
	      else if(neut2.E()==0) neut2 = neut;
	      cand += neut;
	      used_idx.push_back(j);
	      nn++;
	    }
	    else if(neut.DeltaR(gam2)<0.1) {
	      if(neut1.E()==0) neut1 = neut;
	      else if(neut2.E()==0) neut2 = neut;
	      cand += neut;
	      used_idx.push_back(j);
	      nn++;
	    }
	  }
	  if(nn==2) {
	    h1->Fill(cand.M());
	    //h2->Fill(neut1.E()/(gam1+gam2).E());
	    //h2->Fill(cand.DeltaR(cand-neut1));
	    h2->Fill(cand.E()/(gam1+gam2).E());
	  }
	  else if(nn==1) {
	    //h2->Fill(neut1.E()/(gam1+gam2).E());
	  }
	  //h2->Fill((gam1+gam2).Eta());
	  //h2->Fill(gam2.Eta());
	  h3->Fill(nn);
	  //h3->Fill(gam2.E());
	  //if(!(gam1.Eta()>-1.317 && gam1.Eta()<1.901 && gam2.Eta()>-1.317 && gam2.Eta()<1.901)) { h3->Fill(gam1.E()); }
	  if(cand.E()>0.2 && cand.Eta()>-1.317 && cand.Eta()<1.901) {
	    bool overlap = false;
	    for(Int_t k=0; k<el_n; k++) {
	      if(DR(el_eta[k],el_phi[k],cand.Eta(),cand.Phi())<0.2) {
		overlap = true; break;
	      }
	    }
	    if(!overlap) {
	      used_pi0[i] = true;
	      for(std::vector<Int_t>::size_type k=0; k<used_idx.size(); k++) {
		used_neut[used_idx[k]] = true;
	      }
	      bool isPi0 = false;
	      Double_t p = rd.Rndm();
	      if(cand.E()<1.5 && p<0.999)                      isPi0 = true;
	      else if(cand.E()>=1.5 && cand.E()<2.5 && p<0.95) isPi0 = true;
	      else if(cand.E()>=2.5 && cand.E()<3.5 && p<0.92) isPi0 = true;
	      else if(cand.E()>=3.5 && cand.E()<4.5 && p<0.88) isPi0 = true;
	      else if(cand.E()>=4.5 && p<0.79)                 isPi0 = true;
	      if(isPi0) {
		pi0_E[pi0_n] = cand.E()*GeV;
		pi0_pt[pi0_n] = sqrt(pow(cand.E(),2)-pow(0.134977,2))/cosh(cand.Eta())*GeV;
		pi0_eta[pi0_n] = cand.Eta();
		pi0_phi[pi0_n] = cand.Phi();
		pi0_n++;
	      }
	      else {
		if(neut1.E()>0.1) {
		  gam_E[gam_n] = neut1.E()*GeV;
		  gam_pt[gam_n] = neut1.E()/cosh(neut1.Eta())*GeV;
		  gam_eta[gam_n] = neut1.Eta();
		  gam_phi[gam_n] = neut1.Phi();
		  gam_n++;
		}
		if(neut2.E()>0.1) {
		  gam_E[gam_n] = neut2.E()*GeV;
		  gam_pt[gam_n] = neut2.E()/cosh(neut2.Eta())*GeV;
		  gam_eta[gam_n] = neut2.Eta();
		  gam_phi[gam_n] = neut2.Phi();
		  gam_n++;
		}
	      }
	    }
	  }
	}
      }
      if(part->PID==130) {
	Double_t p = rd.Rndm();
	if(p<0.7) nKL++;
      }
    }

    // gamma:
    for(Int_t i=0; i<branchGenParticle->GetEntries(); i++) {
      GenParticle *part = (GenParticle*) branchGenParticle->At(i);    
      if(part->PID==22 && (GetMotherID(branchGenParticle,part,true)!=111 || !used_pi0[GetMotherIdx(branchGenParticle,part,true)]) && part->E>0.05) {
	TLorentzVector truth_gam;
	truth_gam.SetPtEtaPhiM(part->PT, part->Eta, part->Phi, 0.0);
	TLorentzVector cand;
	std::vector<Int_t> used_idx;
	for(Int_t j=0; j<branchNeutral->GetEntries(); j++) {
	  Tower *tower = (Tower*) branchNeutral->At(j);
	  if(used_neut[j]) continue;
	  TLorentzVector neut;
	  neut.SetPtEtaPhiM(tower->ET, tower->Eta, tower->Phi, 0);
	  if(neut.DeltaR(truth_gam)<0.1) {
	    cand += neut;
	    used_idx.push_back(j);
	  }
	}
	if(cand.E()>0.1 && cand.Eta()>-1.317 && cand.Eta()<1.901) {
	  bool overlap = false;
	  for(Int_t k=0; k<el_n; k++) {
	    if(DR(el_eta[k],el_phi[k],cand.Eta(),cand.Phi())<0.2) {
	      overlap = true; break;
	    }
	  }
	  for(Int_t k=0; k<pi0_n; k++) {
	    if(DR(pi0_eta[k],pi0_phi[k],cand.Eta(),cand.Phi())<0.2) {
	      overlap = true; break;
	    }
	  }
	  for(Int_t k=0; k<gam_n; k++) {
	    if(DR(gam_eta[k],gam_phi[k],cand.Eta(),cand.Phi())<0.2) {
	      overlap = true; break;
	    }
	  }
	  if(!overlap) {
	    for(std::vector<Int_t>::size_type k=0; k<used_idx.size(); k++) {
	      used_neut[used_idx[k]] = true;
	    }
	    bool isPi0 = false;
	    Double_t p = rd.Rndm();
	    if(cand.E()<1.5 && p<0.01)                        isPi0 = true;
	    else if(cand.E()>=1.5 && cand.E()<2.5 && p<0.025) isPi0 = true;
	    else if(cand.E()>=2.5 && cand.E()<3.5 && p<0.030) isPi0 = true;
	    else if(cand.E()>=3.5 && p<0.10)                  isPi0 = true;
	    if(isPi0 && cand.E()>0.2) {
	      pi0_E[pi0_n] = cand.E()*GeV;
	      pi0_pt[pi0_n] = sqrt(pow(cand.E(),2)-pow(0.134977,2))/cosh(cand.Eta())*GeV;
	      pi0_eta[pi0_n] = cand.Eta();
	      pi0_phi[pi0_n] = cand.Phi();
	      pi0_n++;
	    }
	    else {
	      gam_E[gam_n] = cand.E()*GeV;
	      gam_pt[gam_n] = cand.E()/cosh(cand.Eta())*GeV;
	      gam_eta[gam_n] = cand.Eta();
	      gam_phi[gam_n] = cand.Phi();
	      gam_n++;
	    }
	  }
	}
      }
    }

    // other neutral:
    for(Int_t i=0; i<branchNeutral->GetEntries(); i++) {
      if(used_neut[i]) continue;
      Tower *tower = (Tower*) branchNeutral->At(i);
      TLorentzVector cand;
      cand.SetPtEtaPhiM(tower->ET, tower->Eta, tower->Phi, 0);
      if(cand.E()>0.1 && cand.Eta()>-1.317 && cand.Eta()<1.901) {
	bool overlap = false;
	for(Int_t k=0; k<el_n; k++) {
	  if(DR(el_eta[k],el_phi[k],cand.Eta(),cand.Phi())<0.2) {
	    overlap = true; break;
	  }
	}
	for(Int_t k=0; k<pi0_n; k++) {
	  if(DR(pi0_eta[k],pi0_phi[k],cand.Eta(),cand.Phi())<0.2) {
	    overlap = true; break;
	  }
	}
	for(Int_t k=0; k<gam_n; k++) {
	  if(DR(gam_eta[k],gam_phi[k],cand.Eta(),cand.Phi())<0.2) {
	    overlap = true; break;
	  }
	}
	if(!overlap) {
	  gam_E[gam_n] = cand.E()*GeV;
	  gam_pt[gam_n] = cand.E()/cosh(cand.Eta())*GeV;
	  gam_eta[gam_n] = cand.Eta();
	  gam_phi[gam_n] = cand.Phi();
	  gam_n++;
	}
      }
    }

    if(el_n!=0 || mu_n!=0 || nK!=0 || nP!=0 || nKL!=0) continue;    
    n_eve[0][1]++;
    n_eve[1][1] += mc_event_weight;

    // reconstruct the ditau topology
    std::vector<TLorentzVector> event;
    std::vector<int> idx;
    std::vector<int> cat;
    for(Int_t i=0; i<trk_n; i++) {
      TLorentzVector cand;
      cand.SetPtEtaPhiE(trk_pt[i]/GeV,trk_eta[i],trk_phi[i],trk_E[i]/GeV);
      event.push_back(cand);
      idx.push_back(i);
      cat.push_back(1);
    }
    for(Int_t i=0; i<pi0_n; i++) {
      TLorentzVector cand;
      cand.SetPtEtaPhiE(pi0_pt[i]/GeV,pi0_eta[i],pi0_phi[i],pi0_E[i]/GeV);
      event.push_back(cand);
      idx.push_back(i);
      cat.push_back(2);
    }
    for(Int_t i=0; i<gam_n; i++) {
      TLorentzVector cand;
      cand.SetPtEtaPhiM(gam_pt[i]/GeV,gam_eta[i],gam_phi[i],0);
      event.push_back(cand);
      idx.push_back(i);
      cat.push_back(0);
    }

    if(event.size()<2) continue;
    n_eve[0][2]++;
    n_eve[1][2] += mc_event_weight;

    TVector3 bv(0,0,-0.272727);
    thrust(event,bv.Px(),bv.Py(),bv.Pz());
    TVector3 thrst(eVec1.px(),eVec1.py(),eVec1.pz()); // thrust axis in CM

    Int_t npi[2] = {0,0};
    Int_t npi0[2] = {0,0};
    std::vector<Int_t> vtau1;
    std::vector<Int_t> vtau2;
    for(std::vector<TLorentzVector>::size_type i=0; i<event.size(); i++) {
      TLorentzVector cand = event[i];
      cand.Boost(bv);
      if(cand.Vect().Dot(thrst)>0) {
	if(cat[i]==1 || cat[i]==2) vtau1.push_back(i);
	if(cat[i]==1) npi[0]++;
	else if(cat[i]==2) npi0[0]++;
      }
      else {
	if(cat[i]==1 || cat[i]==2) vtau2.push_back(i);
	if(cat[i]==1) npi[1]++;
	else if(cat[i]==2) npi0[1]++;
      }
    }

    bool pass = false;
    if(npi[0]==3 && npi0[0]==0 && npi[1]==3 && npi0[1]==0) pass = true;
    if(npi[0]==3 && npi0[0]==0 && npi[1]==1 && npi0[1]==1) pass = true;
    if(npi[0]==3 && npi0[0]==0 && npi[1]==1 && npi0[1]==0) pass = true;
    if(npi[0]==1 && npi0[0]==1 && npi[1]==1 && npi0[1]==1) pass = true;
    if(npi[0]==1 && npi0[0]==1 && npi[1]==1 && npi0[1]==0) pass = true;
    if(npi[0]==1 && npi0[0]==0 && npi[1]==1 && npi0[1]==0) pass = true;
    if(!pass) continue;
    n_eve[0][3]++;
    n_eve[1][3] += mc_event_weight;

    // re-ordering
    if(vtau1.size()<vtau2.size()) {
      std::vector<Int_t> tmp = vtau1;
      vtau1 = vtau2;
      vtau2 = tmp;
    }

    // initialize
    for(Int_t i=0; i<2; i++) {
      tau_trk2_pt[i] = 0;
      tau_trk2_eta[i] = 0;
      tau_trk2_phi[i] = 0;
      tau_trk2_m[i] = 0;
      tau_trk2_x[i] = 0;
      tau_trk2_y[i] = 0;
      tau_trk2_z[i] = 0;
      tau_trk2_xd[i] = 0;
      tau_trk2_yd[i] = 0;
      tau_trk2_zd[i] = 0;
      tau_trk2_charge[i] = 0;
      tau_trk2_pdgId[i] = 0;
      tau_trk3_pt[i] = 0;
      tau_trk3_eta[i] = 0;
      tau_trk3_phi[i] = 0;
      tau_trk3_m[i] = 0;
      tau_trk3_x[i] = 0;
      tau_trk3_y[i] = 0;
      tau_trk3_z[i] = 0;
      tau_trk3_xd[i] = 0;
      tau_trk3_yd[i] = 0;
      tau_trk3_zd[i] = 0;
      tau_trk3_charge[i] = 0;
      tau_trk3_pdgId[i] = 0;
      tau_neut1_pt[i] = 0;
      tau_neut1_eta[i] = 0;
      tau_neut1_phi[i] = 0;
      tau_neut1_m[i] = 0;
    }

    tau_n = 0;
    tau_charge[tau_n] = 0;
    tau_numTrack[tau_n] = 0;
    tau_numNeut[tau_n] = 0;
    TLorentzVector v1;
    for(std::vector<Int_t>::size_type i=0; i<vtau1.size(); i++) {
      v1 += event[vtau1[i]];
      Int_t j = idx[vtau1[i]];
      if(cat[vtau1[i]]==1) {
	if(tau_numTrack[tau_n]==0) {
	  tau_trk1_pt[tau_n] = trk_pt[j];
	  tau_trk1_eta[tau_n] = trk_eta[j];
	  tau_trk1_phi[tau_n] = trk_phi[j];
	  tau_trk1_m[tau_n] = event[vtau1[i]].M()*GeV;
	  tau_trk1_x[tau_n] = trk_x[j];
	  tau_trk1_y[tau_n] = trk_y[j];
	  tau_trk1_z[tau_n] = trk_z[j];
	  tau_trk1_xd[tau_n] = trk_xd[j];
	  tau_trk1_yd[tau_n] = trk_yd[j];
	  tau_trk1_zd[tau_n] = trk_zd[j];
	  tau_trk1_charge[tau_n] = trk_charge[j]>0 ? 1 : -1;
	  tau_trk1_pdgId[tau_n] = trk_charge[j]>0 ? 211 : -211;
	}
	else if(tau_numTrack[tau_n]==1) {
	  tau_trk2_pt[tau_n] = trk_pt[j];
	  tau_trk2_eta[tau_n] = trk_eta[j];
	  tau_trk2_phi[tau_n] = trk_phi[j];
	  tau_trk2_m[tau_n] = event[vtau1[i]].M()*GeV;
	  tau_trk2_x[tau_n] = trk_x[j];
	  tau_trk2_y[tau_n] = trk_y[j];
	  tau_trk2_z[tau_n] = trk_z[j];
	  tau_trk2_xd[tau_n] = trk_xd[j];
	  tau_trk2_yd[tau_n] = trk_yd[j];
	  tau_trk2_zd[tau_n] = trk_zd[j];
	  tau_trk2_charge[tau_n] = trk_charge[j]>0 ? 1 : -1;
	  tau_trk2_pdgId[tau_n] = trk_charge[j]>0 ? 211 : -211;
	}
	else if(tau_numTrack[tau_n]==2) {
	  tau_trk3_pt[tau_n] = trk_pt[j];
	  tau_trk3_eta[tau_n] = trk_eta[j];
	  tau_trk3_phi[tau_n] = trk_phi[j];
	  tau_trk3_m[tau_n] = event[vtau1[i]].M()*GeV;
	  tau_trk3_x[tau_n] = trk_x[j];
	  tau_trk3_y[tau_n] = trk_y[j];
	  tau_trk3_z[tau_n] = trk_z[j];
	  tau_trk3_xd[tau_n] = trk_xd[j];
	  tau_trk3_yd[tau_n] = trk_yd[j];
	  tau_trk3_zd[tau_n] = trk_zd[j];
	  tau_trk3_charge[tau_n] = trk_charge[j]>0 ? 1 : -1;
	  tau_trk3_pdgId[tau_n] = trk_charge[j]>0 ? 211 : -211;
	}
	tau_numTrack[tau_n]++;
	tau_charge[tau_n] += trk_charge[j]>0 ? 1 : -1;
      }
      else if(cat[vtau1[i]]==2) {
	if(tau_numNeut[tau_n]==0) {
	  tau_neut1_pt[tau_n] = pi0_pt[j];
	  tau_neut1_eta[tau_n] = pi0_eta[j];
	  tau_neut1_phi[tau_n] = pi0_phi[j];
	  tau_neut1_m[tau_n] = event[vtau1[i]].M()*GeV;
	}
	tau_numNeut[tau_n]++;
      }
    }
    tau_vis_E[tau_n] = v1.E()*GeV;
    tau_vis_pt[tau_n] = v1.Pt()*GeV;
    tau_vis_eta[tau_n] = v1.Eta();
    tau_vis_phi[tau_n] = v1.Phi();
    tau_n++;

    tau_charge[tau_n] = 0;
    tau_numTrack[tau_n] = 0;
    tau_numNeut[tau_n] = 0;
    TLorentzVector v2;
    for(std::vector<Int_t>::size_type i=0; i<vtau2.size(); i++) {
      v2 += event[vtau2[i]];
      Int_t j = idx[vtau2[i]];
      if(cat[vtau2[i]]==1) {
	if(tau_numTrack[tau_n]==0) {
	  tau_trk1_pt[tau_n] = trk_pt[j];
	  tau_trk1_eta[tau_n] = trk_eta[j];
	  tau_trk1_phi[tau_n] = trk_phi[j];
	  tau_trk1_m[tau_n] = event[vtau2[i]].M()*GeV;
	  tau_trk1_x[tau_n] = trk_x[j];
	  tau_trk1_y[tau_n] = trk_y[j];
	  tau_trk1_z[tau_n] = trk_z[j];
	  tau_trk1_xd[tau_n] = trk_xd[j];
	  tau_trk1_yd[tau_n] = trk_yd[j];
	  tau_trk1_zd[tau_n] = trk_zd[j];
	  tau_trk1_charge[tau_n] = trk_charge[j]>0 ? 1 : -1;
	  tau_trk1_pdgId[tau_n] = trk_charge[j]>0 ? 211 : -211;
	}
	else if(tau_numTrack[tau_n]==1) {
	  tau_trk2_pt[tau_n] = trk_pt[j];
	  tau_trk2_eta[tau_n] = trk_eta[j];
	  tau_trk2_phi[tau_n] = trk_phi[j];
	  tau_trk2_m[tau_n] = event[vtau2[i]].M()*GeV;
	  tau_trk2_x[tau_n] = trk_x[j];
	  tau_trk2_y[tau_n] = trk_y[j];
	  tau_trk2_z[tau_n] = trk_z[j];
	  tau_trk2_xd[tau_n] = trk_xd[j];
	  tau_trk2_yd[tau_n] = trk_yd[j];
	  tau_trk2_zd[tau_n] = trk_zd[j];
	  tau_trk2_charge[tau_n] = trk_charge[j]>0 ? 1 : -1;
	  tau_trk2_pdgId[tau_n] = trk_charge[j]>0 ? 211 : -211;
	}
	else if(tau_numTrack[tau_n]==2) {
	  tau_trk3_pt[tau_n] = trk_pt[j];
	  tau_trk3_eta[tau_n] = trk_eta[j];
	  tau_trk3_phi[tau_n] = trk_phi[j];
	  tau_trk3_m[tau_n] = event[vtau2[i]].M()*GeV;
	  tau_trk3_x[tau_n] = trk_x[j];
	  tau_trk3_y[tau_n] = trk_y[j];
	  tau_trk3_z[tau_n] = trk_z[j];
	  tau_trk3_xd[tau_n] = trk_xd[j];
	  tau_trk3_yd[tau_n] = trk_yd[j];
	  tau_trk3_zd[tau_n] = trk_zd[j];
	  tau_trk3_charge[tau_n] = trk_charge[j]>0 ? 1 : -1;
	  tau_trk3_pdgId[tau_n] = trk_charge[j]>0 ? 211 : -211;
	}
	tau_numTrack[tau_n]++;
	tau_charge[tau_n] += trk_charge[j]>0 ? 1 : -1;
      }
      else if(cat[vtau2[i]]==2) {
	if(tau_numNeut[tau_n]==0) {
	  tau_neut1_pt[tau_n] = pi0_pt[j];
	  tau_neut1_eta[tau_n] = pi0_eta[j];
	  tau_neut1_phi[tau_n] = pi0_phi[j];
	  tau_neut1_m[tau_n] = event[vtau2[i]].M()*GeV;
	}
	tau_numNeut[tau_n]++;
      }
    }
    tau_vis_E[tau_n] = v2.E()*GeV;
    tau_vis_pt[tau_n] = v2.Pt()*GeV;
    tau_vis_eta[tau_n] = v2.Eta();
    tau_vis_phi[tau_n] = v2.Phi();
    tau_n++;

    for(Int_t i=0; i<tau_n; i++) {
      if(tau_numTrack[i]==3) SortTauTrack(i);
    }

    if(abs(tau_charge[0])!=1 || abs(tau_charge[1])!=1 || tau_charge[0]*tau_charge[1]>0) continue;
    n_eve[0][4]++;
    n_eve[1][4] += mc_event_weight;

    truthTau_n = 0;
    Int_t it[2] = {-1,-1};
    Int_t chrg[2] = {0,0};
    for(Int_t i=0; i<branchGenParticle->GetEntries(); i++) {
      GenParticle *part = (GenParticle*) branchGenParticle->At(i);
      if(part && abs(part->PID)==15 && isNotSelf(branchGenParticle,part)) {
	if(it[0]<0) { it[0] = i; chrg[0] = part->Charge; }
	else if(it[1]<0) { it[1] = i; chrg[1] = part->Charge; }
      }
    }
    if(it[0]>=0 && it[1]>=0) {
      if(tau_charge[0]*chrg[0]<0) {
	Int_t tmp = it[0]; it[0] = it[1]; it[1] = tmp;
	tmp = chrg[0]; chrg[0] = chrg[1]; chrg[1] = tmp;
      }
    }

    for(Int_t i=0; i<2; i++) {
      if(it[i]<0) continue;
      GenParticle *dau = (GenParticle*) branchGenParticle->At(it[i]);

      TLorentzVector vis;
      TLorentzVector invis;
      TLorentzVector nu3;
      int id;
      std::vector<TLorentzVector> trk;
      std::vector<TVector3> trk_pos;
      std::vector<TVector3> trk_d;
      std::vector<int> trk_chrg;
      std::vector<int> trk_id;
      std::vector<TLorentzVector> neut;
      std::vector<TVector3> neut_pos;
      std::vector<int> neut_id;
      AnalyzeTau(branchGenParticle,dau,vis,invis,nu3,id,trk,trk_pos,trk_d,trk_chrg,trk_id,neut,neut_pos,neut_id);
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
	truthTau_trk1_charge[truthTau_n] = trk_chrg[itrk[0]];
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
	truthTau_trk2_charge[truthTau_n] = trk_chrg[itrk[1]];
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
	truthTau_trk3_charge[truthTau_n] = trk_chrg[itrk[2]];
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
    }

    // truth plot
    if(truthTau_n>=2) {
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
      h4->Fill((tau1+tau2).M()/GeV);
    }

    tree->Fill();
  }
  printf(" Tot %1.0f, Sel %1.0f %1.0f %1.0f %1.0f\n",n_eve[0][0],n_eve[0][1],n_eve[0][2],n_eve[0][3],n_eve[0][4]);
  printf(" Tot %g, Sel %g %g %g %g\n",n_eve[1][0],n_eve[1][1],n_eve[1][2],n_eve[1][3],n_eve[1][4]);

  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  h1->SetXTitle("m_{#pi0} [GeV]");
  h1->SetMinimum(0);
  h1->Draw();

  TCanvas* c2 = new TCanvas("c2","c2",600,600);
  h2->SetXTitle("#DeltaR(#gamma#gamma)");
  h2->SetMinimum(0);
  h2->Draw();

  TCanvas* c3 = new TCanvas("c3","c3",600,600);
  h3->SetXTitle("N_{#gamma}");
  h3->SetMinimum(0);
  h3->Draw();

  TCanvas* c4 = new TCanvas("c4","c4",600,600);
  h4->SetXTitle("m_{#tau#tau} [GeV]");
  h4->SetMinimum(0);
  h4->Draw();

  fout->cd();
  tree->AutoSave();
}
