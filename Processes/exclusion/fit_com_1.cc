#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TRandom1.h"
#include "TLegend.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <vector>
#include <sstream> //istringstream
#include "../Elizabeth.h"

using namespace std;

Double_t f(Double_t x, Double_t* s, Double_t* b, Double_t* n) {
  return n[0]*s[0]/(s[0]*x+b[0]) + n[1]*s[1]/(s[1]*x+b[1]) + n[2]*s[2]/(s[2]*x+b[2]) - (s[0]+s[1]+s[2]);
}

Double_t fprime(Double_t x, Double_t* s, Double_t* b, Double_t* n) {
  return -n[0]*pow(s[0]/(s[0]*x+b[0]),2)-n[1]*pow(s[1]/(s[1]*x+b[1]),2)-n[2]*pow(s[2]/(s[2]*x+b[2]),2);
}
// get the solution using Newton's method
Double_t getMuBest(Double_t* s, Double_t* b, Double_t* n) {
  Double_t lo = TMath::Min(-b[0]/s[0],-b[1]/s[1]);
  lo = TMath::Min(lo,-b[2]/s[2])+1e-4;
  Double_t x1 = 1;
  Double_t x2 = 0;
  while(fabs(x2-x1)>1e-4) {
    x1 = x2;
    x2 = x1 - f(x1,s,b,n)/fprime(x1,s,b,n);
    //if(x2<lo) printf("Error: solution out of range in Newton's iteration\n");
  }
  return x2;
}

Double_t getSigmaExcl(Double_t* s, Double_t* b, Double_t* n) {
  Double_t mu_best = getMuBest(s,b,n);
  Double_t lnL_best = n[0]*log(s[0]*mu_best+b[0])-(s[0]*mu_best+b[0]) +
    n[1]*log(s[1]*mu_best+b[1])-(s[1]*mu_best+b[1]) +
    n[2]*log(s[2]*mu_best+b[2])-(s[2]*mu_best+b[2]);
  Double_t lnL_1 = n[0]*log(s[0]+b[0])-(s[0]+b[0]) +
    n[1]*log(s[1]+b[1])-(s[1]+b[1]) +
    n[2]*log(s[2]+b[2])-(s[2]+b[2]);
  return sqrt(2*(lnL_best-lnL_1))*(mu_best>1. ? -1. : 1.);
}

void Fill(TH1F* h, Double_t x) {
  Double_t tmp = x;
  if(x<h->GetXaxis()->GetXmin()) tmp = h->GetXaxis()->GetXmin();
  if(x>=h->GetXaxis()->GetXmax()) tmp = h->GetXaxis()->GetXmax()*0.999999;
  h->Fill(tmp);
}

string to_string(int val) {
    char buf[200];
    sprintf(buf, "%d", val);
    return string(buf);
}
 

void plot() {

ofstream fout("alpha_1.txt");
fout<<"file"<<" "<<""<<"mass"<<" "<<"fw"<<" "<<"fww"<<" "<<"alpha"<<" "<<"sigma1"<<" "<<"sigma2"<<endl;

   ifstream fr("../../../param_2lep600.txt");
  if(!fr.is_open()){
  cout<<"unable to open file"<<endl;
 }

 vector<string> vec;
 Double_t cs;
 Double_t fw;
 Double_t fww;
 Double_t mass;
 Int_t num;
 Int_t i=0;

 while(!fr.eof())     //while 循环
 {
  string temp;
  getline(fr, temp, '\n');
  vec.push_back(temp);
  i ++;
}
 int Row=0; //行

//  //for(auto it=vec.begin(); it != vec.end(); it++)
 for(std::vector<string>::iterator it=vec.begin(); it != vec.end(); it++)
  {  
    // cout << *it << endl;
    istringstream is(*it);
    string s;
    Row++;
    num=0;
    int column=0; //列

    //if(  (Row>1 && Row<46) ){
    if(  (Row>1 && Row<212) ){

    while(is>>s)
    {
      
      if(column == 0 )
      {
        num=atof(s.c_str());
      }
      if(column == 1 )
      {
        mass=atof(s.c_str());
      }

      if(column == 2 )
      {
        fw=atof(s.c_str());
      }
      if(column == 3 )
      {
        fww=atof(s.c_str());
      }

      if(column == 4 )
      {
        cs=atof(s.c_str());
      }

      column++;
    }

    string s_fw = to_string(int(fabs(fw)));
    string s_fww = to_string(int(fabs(fww)));
    string in_fw = "";
    string in_fww = "";
    if(fw<0.)  in_fw = "m";
    if(fww<0.)  in_fww = "m";
    TString hist_sig_name = "hsig_fw_"+in_fw+s_fw+"_fww_"+in_fww+s_fww;

    cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<:"<<hist_sig_name<<endl;

  
  char str[200];
  TFile *file_2_bkg = TFile::Open("../../ntuple/mHH_600_hist_bkg_2lep.root");
  if (!file_2_bkg) return;
  TFile *file_2_sig = TFile::Open("../../ntuple/mHH_600_sig_2lep_005.root");//30,40,50,60,70
  if (!file_2_sig) return;

  TH1F* hsig_1 = (TH1F*) file_2_sig->Get(hist_sig_name);
  TH1F* hbkg_1 = (TH1F*) file_2_bkg->Get("hbkg1");
  hsig_1->SetName("hsig_1");
  hbkg_1->SetName("hbkg_1");
  hsig_1->Scale(140./300.);
  hbkg_1->Scale(140./300.);
  hsig_1->SetDirectory(0);
  hbkg_1->SetDirectory(0);
  Double_t nsig_exp_1 = hsig_1->Integral();
  Double_t nbkg_exp_1 = hbkg_1->Integral();

  TFile *file_3_bkg = TFile::Open("../../ntuple/mHH_600_hist_bkg_3lep.root");
  if (!file_3_bkg) return;
  TFile *file_3_sig = TFile::Open("../../ntuple/mHH_600_sig_3lep_005.root");//30,40,50,60
  if (!file_3_sig) return;

  TH1F* hsig_2 = (TH1F*) file_3_sig->Get(hist_sig_name);
  TH1F* hbkg_2 = (TH1F*) file_3_bkg->Get("hbkg1");
  TH1F* hbkg2_2 = (TH1F*) file_3_bkg->Get("hbkg2");
  hsig_2->SetName("hsig_2");
  hbkg_2->SetName("hbkg_2");
  hbkg2_2->SetName("hbkg2_2");
  hsig_2->Scale(140./300.);
  hbkg_2->Scale(140./300.);
  hbkg2_2->Scale(140./300.);
  hsig_2->SetDirectory(0);
  hbkg_2->SetDirectory(0);
  hbkg2_2->SetDirectory(0);
  hbkg_2->Add(hbkg2_2);
  Double_t nsig_exp_2 = hsig_2->Integral();
  Double_t nbkg_exp_2 = hbkg_2->Integral();
  ////////////////////////// www ////////////
  TFile *file_4_bkg = TFile::Open("../../ntuple/mHH_600_hist_bkg_www.root");
  if (!file_4_bkg) return;
  TFile *file_4_sig = TFile::Open("../../ntuple/mHH_600_sig_www_005.root");//30,40,50,60
  if (!file_4_sig) return;

  TH1F* hsig_3 = (TH1F*) file_4_sig->Get(hist_sig_name);
  TH1F* hbkg_3 = (TH1F*) file_4_bkg->Get("hbkg1");
  TH1F* hbkg2_3 = (TH1F*) file_4_bkg->Get("hbkg2");
  TH1F* hbkg3_3 = (TH1F*) file_4_bkg->Get("hbkg3");
  hsig_3->SetName("hsig_3");
  hbkg_3->SetName("hbkg_3");
  hbkg2_3->SetName("hbkg2_3");
  hbkg3_3->SetName("hbkg3_3");
  hsig_3->Scale(140./300.);
  hbkg_3->Scale(140./300.);
  hbkg2_3->Scale(140./300.);
  hbkg3_3->Scale(140./300.);
  hsig_3->SetDirectory(0);
  hbkg_3->SetDirectory(0);
  hbkg2_3->SetDirectory(0);
  hbkg3_3->SetDirectory(0);
  hbkg_3->Add(hbkg2_3);
  hbkg_3->Add(hbkg3_3);
  Double_t nsig_exp_3 = hsig_3->Integral();
  Double_t nbkg_exp_3 = hbkg_3->Integral();
  //////////////////////////////////////////
  Double_t s[3] = {nsig_exp_1, nsig_exp_2, nsig_exp_3};
  Double_t b[3] = {nbkg_exp_1, nbkg_exp_2, nbkg_exp_3};
  Double_t n[3];
  TH1F* h1 = new TH1F("h1","",1000,-3,5); // b-only
  TH1F* h2 = new TH1F("h2","",1000,-3,5); // s+b

  TRandom1 rdm;

  for(Int_t i=0; i<10000; i++) {
    n[0] = rdm.Poisson(b[0]);
    n[1] = rdm.Poisson(b[1]);
    n[2] = rdm.Poisson(b[2]);
    Fill(h1,getSigmaExcl(s,b,n));
  }

  for(Int_t i=0; i<10000; i++) {
    n[0] = rdm.Poisson(b[0]+s[0]);
    n[1] = rdm.Poisson(b[1]+s[1]);
    n[2] = rdm.Poisson(b[2]+s[2]);
    Fill(h2,getSigmaExcl(s,b,n));
  }
  
  h2->SetLineColor(kBlue);
  TLegend* lg = new TLegend(0.77,0.74,0.90,0.89,"");
  lg->AddEntry(h1,"B-only","L");
  lg->AddEntry(h2,"S+B","L");
  lg->SetBorderSize(0);
  lg->SetMargin(0.25);
  lg->SetFillColor(kWhite);
  
  TCanvas* c1 = new TCanvas("c1","c1",1000,500);
  h1->Draw("hist");
  h2->Draw("hist same");
  lg->Draw("same");

  // Double_t F[1] = {0.5};
  //   // Double_t q[1];
  //     // h1->GetQuantiles(1,q,F);
  Int_t iqS1=0;
  Int_t iqS2=0;
  Int_t iq=0;
  Double_t minS1 = 999;
  Double_t minS2 = 999;
  Double_t min = 999;
  Double_t ptmp = SUSYStat_Pval(1);
  Double_t xq[3] = {ptmp,0.5,1-ptmp};

  for(Int_t i=1; i<h1->GetNbinsX(); i++) {
    if(fabs(h1->Integral(1,i)/h1->Integral()-xq[0])<minS1) {
        minS1 = fabs(h1->Integral(1,i)/h1->Integral()-xq[0]);
        iqS1 = i;
    }
  }
  for(Int_t i=1; i<h1->GetNbinsX(); i++) {
    if(fabs(h1->Integral(1,i)/h1->Integral()-xq[2])<minS2) {
        minS2 = fabs(h1->Integral(1,i)/h1->Integral()-xq[2]);
        iqS2 = i;
    }
  }

  for(Int_t i=1; i<h1->GetNbinsX(); i++) {
    if(fabs(h1->Integral(1,i)/h1->Integral()-0.5)<min) {
        min = fabs(h1->Integral(1,i)/h1->Integral()-0.5);
        iq = i;
    }
  }
  Double_t q[3] = {h1->GetBinCenter(iqS1),h1->GetBinCenter(iq),h1->GetBinCenter(iqS2)};
  printf("q  = %f\n",q[1]);
  TLine* line2lo = new TLine(q[0],0,q[0],h2->GetMaximum());
  line2lo->SetLineWidth(1);
  line2lo->SetLineColor(kRed);
  line2lo->SetLineStyle(2);
  line2lo->Draw();
  TLine* line2 = new TLine(q[1],0,q[1],h2->GetMaximum());
  line2->SetLineWidth(2);
  line2->SetLineColor(kRed);
  line2->Draw();
  TLine* line2hi = new TLine(q[2],0,q[2],h2->GetMaximum());
  line2hi->SetLineWidth(1);
  line2hi->SetLineColor(kRed);
  line2hi->SetLineStyle(2);
  line2hi->Draw();
                                          
  printf("Excl. CL = %%%.1f\n",h2->Integral(1,h2->FindBin(q[1]))/h2->Integral()*100);
  Double_t alpha = 1.-h2->Integral(1,h2->FindBin(q[1]))/h2->Integral();
  Double_t sigma1 = 1.-h2->Integral(1,h2->FindBin(q[0]))/h2->Integral();
  Double_t sigma2 = 1.-h2->Integral(1,h2->FindBin(q[2]))/h2->Integral();
  fout<<num<<" "<<mass<<" "<<fw<<" "<<fww<<" "<<alpha<<" "<<sigma1<<" "<<sigma2<<endl; 
  c1->SaveAs(hist_sig_name+"_h12.eps");
}
}
}
