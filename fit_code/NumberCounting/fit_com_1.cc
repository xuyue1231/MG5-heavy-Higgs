#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TRandom.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <vector>
#include <sstream> //istringstream
#include "../Elizabeth.h"

using namespace std;

string to_string(int val) {
    char buf[200];
    sprintf(buf, "%d", val);
    return string(buf);
}
 

void plot() {

ofstream fout("alpha_1.txt");
fout<<"file"<<" "<<""<<"mass"<<" "<<"fw"<<" "<<"fww"<<" "<<"alpha"<<" "<<"sigma1"<<" "<<"sigma2"<<endl;

   ifstream fr("../../parameter/param_2lep900.txt");
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
    if(  (Row>1 && Row<150) ){

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
  TFile *file_2_bkg = TFile::Open("../../ntuple/mHH_900_hist_bkg_2lep_up.root");
  if (!file_2_bkg) return;
  TFile *file_2_sig = TFile::Open("../../ntuple/mHH_900_sig_2lep_005.root");//30,40,50,60,70
  if (!file_2_sig) return;

  TH1F* hsig_1 = (TH1F*) file_2_sig->Get(hist_sig_name);
  TH1F* hbkg_1 = (TH1F*) file_2_bkg->Get("hbkg1");
  hsig_1->Scale(140./300.);
  hbkg_1->Scale(140./300.);
  hsig_1->SetName("hsig_1");
  hbkg_1->SetName("hbkg_1");
  hsig_1->SetDirectory(0);
  hbkg_1->SetDirectory(0);
  Double_t nsig_exp_1 = hsig_1->Integral();
  Double_t nbkg_exp_1 = hbkg_1->Integral();

  TFile *file_3_bkg = TFile::Open("../../ntuple/mHH_900_hist_bkg_3lep.root");
  if (!file_3_bkg) return;
  TFile *file_3_sig = TFile::Open("../../ntuple/mHH_900_sig_3lep_005.root");//30,40,50,60
  if (!file_3_sig) return;

  TH1F* hsig_2 = (TH1F*) file_3_sig->Get(hist_sig_name);
  TH1F* hbkg_2 = (TH1F*) file_3_bkg->Get("hbkg1");
  TH1F* hbkg2_2 = (TH1F*) file_3_bkg->Get("hbkg2");
  hsig_2->Scale(140./300.);
  hbkg_2->Scale(140./300.);
  hbkg2_2->Scale(140./300.);
  hsig_2->SetName("hsig_2");
  hbkg_2->SetName("hbkg_2");
  hbkg2_2->SetName("hbkg2_2");
  hsig_2->SetDirectory(0);
  hbkg_2->SetDirectory(0);
  hbkg2_2->SetDirectory(0);
  hbkg_2->Add(hbkg2_2);
  Double_t nsig_exp_2 = hsig_2->Integral();
  Double_t nbkg_exp_2 = hbkg_2->Integral();
  ////////////////////////// www ////////////
  TFile *file_4_bkg = TFile::Open("../../ntuple/mHH_900_hist_bkg_www.root");
  if (!file_4_bkg) return;
  TFile *file_4_sig = TFile::Open("../../ntuple/mHH_900_sig_www_005.root");//30,40,50,60
  if (!file_4_sig) return;

  TH1F* hsig_3 = (TH1F*) file_4_sig->Get(hist_sig_name);
  TH1F* hbkg_3 = (TH1F*) file_4_bkg->Get("hbkg1");
  TH1F* hbkg2_3 = (TH1F*) file_4_bkg->Get("hbkg2");
  TH1F* hbkg3_3 = (TH1F*) file_4_bkg->Get("hbkg3");
  hsig_3->Scale(140./300.);
  hbkg_3->Scale(140./300.);
  hbkg2_3->Scale(140./300.);
  hbkg3_3->Scale(140./300.);
  hsig_3->SetName("hsig_3");
  hbkg_3->SetName("hbkg_3");
  hbkg2_3->SetName("hbkg2_3");
  hbkg3_3->SetName("hbkg3_3");
  hsig_3->SetDirectory(0);
  hbkg_3->SetDirectory(0);
  hbkg2_3->SetDirectory(0);
  hbkg3_3->SetDirectory(0);
  hbkg_3->Add(hbkg2_3);
  hbkg_3->Add(hbkg3_3);
  Double_t nsig_exp_3 = hsig_3->Integral();
  Double_t nbkg_exp_3 = hbkg_3->Integral();
  //////////////////////////////////////////
  Int_t ntoy = 30000;
  TH1F* h1 = new TH1F("h1","",600,0,1000);
  TH1F* h2 = new TH1F("h2","",600,0,1000);
  Float_t num_s = nsig_exp_1 + nsig_exp_2 + nsig_exp_3;
  Float_t num_bkg = nbkg_exp_1 + nbkg_exp_2 + nbkg_exp_3;
  Float_t num_sig_bkg = num_s+num_bkg;
  for(int i=0; i<ntoy; i++){
  	h1->Fill(gRandom->PoissonD(num_bkg));
  	h2->Fill(gRandom->PoissonD(num_sig_bkg));
  }
  Double_t ptmp = SUSYStat_Pval(1);
  Double_t xq[3] = {ptmp,0.5,1-ptmp};
  Double_t yq[3];
  h1->GetQuantiles(3, yq, xq);
  printf("yq[1] = %f\n",yq[1]);

  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  c1->SetLogy();
  h1->Draw();
  TLine* line1lo = new TLine(yq[0],0,yq[0],h1->GetMaximum());
  line1lo->SetLineWidth(1);
  line1lo->SetLineColor(kRed);
  line1lo->SetLineStyle(2);
  line1lo->Draw();
  TLine* line1 = new TLine(yq[1],0,yq[1],h1->GetMaximum());
  line1->SetLineWidth(2);
  line1->SetLineColor(kRed);
  line1->Draw();
  TLine* line1hi = new TLine(yq[2],0,yq[2],h1->GetMaximum());
  line1hi->SetLineWidth(1);
  line1hi->SetLineColor(kRed);
  line1hi->SetLineStyle(2);
  line1hi->Draw();
  //c1->SaveAs(hist_sig_name+"_h1_.eps");

  TCanvas* c2 = new TCanvas("c2","c2",600,600);
  c2->SetLogy();
  h2->Draw();
  TLine* line2lo = new TLine(yq[0],0,yq[0],h2->GetMaximum());
  line2lo->SetLineWidth(1);
  line2lo->SetLineColor(kRed);
  line2lo->SetLineStyle(2);
  line2lo->Draw();
  TLine* line2 = new TLine(yq[1],0,yq[1],h2->GetMaximum());
  line2->SetLineWidth(2);
  line2->SetLineColor(kRed);
  line2->Draw();
  TLine* line2hi = new TLine(yq[2],0,yq[2],h2->GetMaximum());
  line2hi->SetLineWidth(1);
  line2hi->SetLineColor(kRed);
  line2hi->SetLineStyle(2);
  line2hi->Draw();
 //c2->SaveAs(hist_sig_name+"_h2_.eps");
  printf("%f %f %f\n",h2->Integral(1,h2->FindBin(yq[0]))/ntoy,h2->Integral(1,h2->FindBin(yq[1]))/ntoy,h2->Integral(1,h2->FindBin(yq[2]))/ntoy);
  Double_t alpha = h2->Integral(1,h2->FindBin(yq[1]))/ntoy;
  Double_t sigma1 = h2->Integral(1,h2->FindBin(yq[0]))/ntoy;
  Double_t sigma2 = h2->Integral(1,h2->FindBin(yq[2]))/ntoy;
  fout<<num<<" "<<mass<<" "<<fw<<" "<<fww<<" "<<alpha<<" "<<sigma1<<" "<<sigma2<<endl;
}
}
}
