#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLine.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooProduct.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooAddPdf.h"
#include "RooHistPdf.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooPlot.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <vector>
#include <sstream> //istringstream
#include "../Elizabeth.h"

using namespace std;
using namespace RooFit;

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
    if(  (Row>1 && Row<7) ){

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

  RooRealVar mVV("mVV","m_{VV} [GeV]",150,2000);
  mVV.setBins(hsig_1->GetNbinsX());

  RooRealVar weight("weight","",0,200);

  RooCategory cat("cat","");;
  cat.defineType("2lep");
  cat.defineType("3lep");
  cat.defineType("www");
  
  RooDataHist dsig_1("dsig_1","",RooArgList(mVV),hsig_1);
  RooDataHist dbkg_1("dbkg_1","",RooArgList(mVV),hbkg_1);
  RooHistPdf psig_1("psig_1","",RooArgSet(mVV),dsig_1);
  RooHistPdf pbkg_1("pbkg_1","",RooArgSet(mVV),dbkg_1);
  RooRealVar nsig_1_("nsig_1_","",nsig_exp_1);
  RooRealVar nbkg_1("nbkg_1","",nbkg_exp_1);

  RooDataHist dsig_2("dsig_2","",RooArgList(mVV),hsig_2);
  RooDataHist dbkg_2("dbkg_2","",RooArgList(mVV),hbkg_2);
  RooHistPdf psig_2("psig_2","",RooArgSet(mVV),dsig_2);
  RooHistPdf pbkg_2("pbkg_2","",RooArgSet(mVV),dbkg_2);
  RooRealVar nsig_2_("nsig_2_","",nsig_exp_2);
  RooRealVar nbkg_2("nbkg_2","",nbkg_exp_2);

  RooDataHist dsig_3("dsig_3","",RooArgList(mVV),hsig_3);
  RooDataHist dbkg_3("dbkg_3","",RooArgList(mVV),hbkg_3);
  RooHistPdf psig_3("psig_3","",RooArgSet(mVV),dsig_3);
  RooHistPdf pbkg_3("pbkg_3","",RooArgSet(mVV),dbkg_3);
  RooRealVar nsig_3_("nsig_3_","",nsig_exp_3);
  RooRealVar nbkg_3("nbkg_3","",nbkg_exp_3);

  RooRealVar mu("mu","",-3,3);
  RooProduct nsig_1("nsig_1","",RooArgList(nsig_1_,mu));
  RooProduct nsig_2("nsig_2","",RooArgList(nsig_2_,mu));
  RooProduct nsig_3("nsig_3","",RooArgList(nsig_3_,mu));

  RooAddPdf* p1 = new RooAddPdf("p1","",RooArgList(pbkg_1,psig_1),RooArgList(nbkg_1,nsig_1));
  RooAddPdf* p2 = new RooAddPdf("p2","",RooArgList(pbkg_2,psig_2),RooArgList(nbkg_2,nsig_2));
  RooAddPdf* p3 = new RooAddPdf("p3","",RooArgList(pbkg_3,psig_3),RooArgList(nbkg_3,nsig_3));

  RooSimultaneous pall("pall","",cat) ;//facilitates simultaneous fitting of multiple PDFs to subsets of a given dataset.
  pall.addPdf(*p1,"2lep") ;
  pall.addPdf(*p2,"3lep") ;
  pall.addPdf(*p3,"www") ;

  TH1F* h1 = new TH1F("h1","",600,-3,3);
  TH1F* h2 = new TH1F("h2","",600,-3,3);
  
  mu.setVal(1);
  nbkg_1.setVal(nbkg_exp_1);
  nbkg_2.setVal(nbkg_exp_2);
  nbkg_3.setVal(nbkg_exp_3);
  RooDataSet* d1;
  RooDataSet* d2;
  RooDataSet* d3;
  RooDataSet* dall;
  RooArgSet* pars1 = (RooArgSet*) p1->getParameters(RooArgSet(mVV))->snapshot();
  RooArgSet* pars2 = (RooArgSet*) p2->getParameters(RooArgSet(mVV))->snapshot();
  RooArgSet* pars3 = (RooArgSet*) p3->getParameters(RooArgSet(mVV))->snapshot();
  ///对2lep和3lep作为两个subset，用RooCategory区分两个子集，fit的时候是两个子集分别fit
  Int_t ntoy = 10000;
  for(Int_t i=0; i<ntoy; i++) {
    (*p1->getParameters(RooArgSet(mVV))) = (*pars1);
    (*p2->getParameters(RooArgSet(mVV))) = (*pars2);
    (*p3->getParameters(RooArgSet(mVV))) = (*pars3);
    mu.setVal(0);
    d1 = p1->generate(RooArgSet(mVV),AutoBinned(true),Extended());
    d2 = p2->generate(RooArgSet(mVV),AutoBinned(true),Extended());
    d3 = p3->generate(RooArgSet(mVV),AutoBinned(true),Extended());
    dall = new RooDataSet("dall","",RooArgSet(mVV,cat,weight),WeightVar("weight"));
    for(Int_t i=0; i<d1->numEntries(); i++) {
      const RooArgSet* row = d1->get(i);
      cat.setLabel("2lep");
      RooRealVar* var = (RooRealVar*) row->find("mVV");
      mVV.setVal(var->getVal());
      dall->add(RooArgSet(mVV,cat),d1->weight());
    } 
    for(Int_t i=0; i<d2->numEntries(); i++) {
      const RooArgSet* row = d2->get(i);
      cat.setLabel("3lep");
      RooRealVar* var = (RooRealVar*) row->find("mVV");
      mVV.setVal(var->getVal());
      dall->add(RooArgSet(mVV,cat),d2->weight());
    }
    for(Int_t i=0; i<d3->numEntries(); i++) {
      const RooArgSet* row = d3->get(i);
      cat.setLabel("www");
      RooRealVar* var = (RooRealVar*) row->find("mVV");
      mVV.setVal(var->getVal());
      dall->add(RooArgSet(mVV,cat),d3->weight());
    }
    pall.fitTo(*dall,Extended(true),SumW2Error(true));
    if(mu.getVal()<h1->GetXaxis()->GetXmin()) h1->Fill(h1->GetXaxis()->GetXmin());
    else if(mu.getVal()>=h1->GetXaxis()->GetXmax()) h1->Fill(h1->GetXaxis()->GetXmax()*0.999);
    else h1->Fill(mu.getVal());
    delete d1;
    delete d2;
    delete d3;
    delete dall;

    (*p1->getParameters(RooArgSet(mVV))) = (*pars1);
    (*p2->getParameters(RooArgSet(mVV))) = (*pars2);
    (*p3->getParameters(RooArgSet(mVV))) = (*pars3);
    mu.setVal(1);
    d1 = p1->generate(RooArgSet(mVV),AutoBinned(true),Extended());
    d2 = p2->generate(RooArgSet(mVV),AutoBinned(true),Extended());
    d3 = p3->generate(RooArgSet(mVV),AutoBinned(true),Extended());
    dall = new RooDataSet("dall","",RooArgSet(mVV,cat,weight),WeightVar("weight"));
    for(Int_t i=0; i<d1->numEntries(); i++) {
      const RooArgSet* row = d1->get(i);
      cat.setLabel("2lep");
      RooRealVar* var = (RooRealVar*) row->find("mVV");
      mVV.setVal(var->getVal());
      dall->add(RooArgSet(mVV,cat),d1->weight());
    }
    for(Int_t i=0; i<d2->numEntries(); i++) {
      const RooArgSet* row = d2->get(i);
      cat.setLabel("3lep");
      RooRealVar* var = (RooRealVar*) row->find("mVV");
      mVV.setVal(var->getVal());
      dall->add(RooArgSet(mVV,cat),d2->weight());
    }
    for(Int_t i=0; i<d3->numEntries(); i++) {
      const RooArgSet* row = d3->get(i);
      cat.setLabel("www");
      RooRealVar* var = (RooRealVar*) row->find("mVV");
      mVV.setVal(var->getVal());
      dall->add(RooArgSet(mVV,cat),d3->weight());
    }
    pall.fitTo(*dall,Extended(true),SumW2Error(true));
    if(mu.getVal()<h2->GetXaxis()->GetXmin()) h2->Fill(h2->GetXaxis()->GetXmin());
    else if(mu.getVal()>=h2->GetXaxis()->GetXmax()) h2->Fill(h2->GetXaxis()->GetXmax()*0.999);
    else h2->Fill(mu.getVal());
    delete d1;
    delete d2;
    delete d3;
    delete dall;
  }

  Double_t ptmp = SUSYStat_Pval(1);
  Double_t xq[3] = {ptmp,0.5,1-ptmp};
  Double_t yq[3];
  h1->GetQuantiles(3,yq,xq);
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
//  c1->SaveAs(hist_sig_name+"_h1_.eps");

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
 // c2->SaveAs(hist_sig_name+"_h2_.eps");
  printf("%f %f %f\n",h2->Integral(1,h2->FindBin(yq[0]))/ntoy,h2->Integral(1,h2->FindBin(yq[1]))/ntoy,h2->Integral(1,h2->FindBin(yq[2]))/ntoy);
  Double_t alpha = h2->Integral(1,h2->FindBin(yq[1]))/ntoy;
  Double_t sigma1 = h2->Integral(1,h2->FindBin(yq[0]))/ntoy;
  Double_t sigma2 = h2->Integral(1,h2->FindBin(yq[2]))/ntoy;
  fout<<num<<" "<<mass<<" "<<fw<<" "<<fww<<" "<<alpha<<" "<<sigma1<<" "<<sigma2<<endl;
}
}
}
