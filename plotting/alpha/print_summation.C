#include <iostream>
#include <vector>
#include <algorithm>
#include "TCanvas.h"
#include "TError.h"
#include "TPad.h"
#include "TString.h"
#include "TRandom1.h"
#include "TLorentzVector.h"
#include "TH1F.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TCut.h"
#include "TNtuple.h"
#include "TLine.h"  

void drawText(const char *text, float xp, float yp, int size=22){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(43);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw("same");
}

void print_summation(){
 TH1D::SetDefaultSumw2();
 int domp=0;
 int gensigbkgreco=3;
 int cent_min=60;
 int cent_max=200;
 double etadijet=0.5; 
 int doMC=0; 
 int docorr=1;
 if(gensigbkgreco!=3) docorr=0; 
 int jetpt1=120;
 int jetpt2=30;
 
 TString smpt[]={"mpt","mp","mpt_boosted"};
 TString skind[]={"","_s","_b","_tracks","tracks_uncorr"};
 TString spart[]={"gen. part.","gen. sig.","gen. part.","corr reco. part.","uncorr. reco."};


 int npt=6;
 double ptmin[]={0.5,1,2,4,8  ,0.5};
 double ptmax[]={1  ,2,4,8,100,100};
 TFile *f[npt];
 TTree *t[npt];
 TTree *t_jet[npt];
 int nalpha=30;
 double mpt[npt][nalpha];
 double mpterr[npt][nalpha];
 double mptsum[nalpha];
 double mptsumpos[nalpha];
 double mptsumneg[nalpha];
 TH1D* h[npt][nalpha];
 TH1D *hsum[nalpha];
 TH1D * hmpt_vs_alpha[npt];
 // int col[]={kAzure-9,kYellow-9,kOrange+1,kGreen+3,kPink-2,kRed+1,1};
 int col[]={kAzure-9,kYellow-9,kOrange+1,kGreen+3,kRed+1,1};

 for(int ialpha=0;ialpha<nalpha;ialpha++){
  mptsum[ialpha]=0;
  mptsumpos[ialpha]=0;
  mptsumneg[ialpha]=0;
 }
   

 TLegend *leg,*leg2;
 if(domp==0)leg= new TLegend(0.28,0.45,0.56,0.65);
 if(domp>0)leg= new TLegend(0.64,0.65,0.98,0.8);
 leg->SetBorderSize(0); 
 leg->SetFillStyle(0);
 if(domp==0)leg2 = new TLegend(0.64,0.45,0.98,0.65);
 if(domp>0)leg2 = new TLegend(0.64,0.8,1.,0.97);
 leg2->SetBorderSize(0); 
 leg2->SetFillStyle(0); 
 
  double frac[]={0,0.01,0.02,0.03,0.04,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.925,0.95,0.96,0.97,0.98,0.99,1.01};
  
 for(int ipt=0;ipt<npt;ipt++){  
  if(!doMC){ 
   f[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_data_20140310/full_ntuple_HiForest_PbPb_Jet80_Track8_Jet17_GR_R_53_LV6_merged_forest_0_pt%d_%d_akPu3Calo.root",(int)ptmin[ipt],(int)ptmax[ipt])); 
   // if(!docorr)f[ipt] = TFile::Open(Form("root://eoscms//eos/cms/store/caf/user/dgulhan/missingPt/ntuples_data/full_ntuple_HiForest_PbPb_Jet80_Track8_Jet17_GR_R_53_LV6_merged_forest_0_pt%d_%d_akPu3Calo.root",(int)ptmin[ipt],(int)ptmax[ipt]));   
  }
  if(doMC){
   // if(!docorr)f[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_MC_uncorr/full_ntuple_HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_0_pt%d_%d_akPu3Calo.root",(int)ptmin[ipt],(int)ptmax[ipt])); 
   if(docorr)f[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_MC_20140310/full_ntuple_HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_0_pt%d_%d_akPu3Calo.root",(int)ptmin[ipt],(int)ptmax[ipt])); 
  }
  t[ipt]=(TTree*)f[ipt]->Get(Form("nt_%s%s",smpt[domp].Data(),skind[gensigbkgreco].Data()));
  t_jet[ipt]=(TTree*)f[ipt]->Get("nt_jet");
  t[ipt]->AddFriend(t_jet[ipt]);
  // hmpt_vs_alpha[ipt]= new TH1D(Form("hmpt_%d",ipt),"",nalpha,0.025,1.025);
  // hmpt_vs_alpha[ipt]= new TH1D(Form("hmpt_%d",ipt),"",nalpha/2,0.025,0.525);
  hmpt_vs_alpha[ipt]= new TH1D(Form("hmpt_%d",ipt),"",nalpha/2,frac);
  // for(int ialpha=0;ialpha<nalpha;ialpha++ ){
  for(int ialpha=0;ialpha<nalpha/2;ialpha++ ){
   h[ipt][ialpha]=new TH1D(Form("h_%d_%d",ialpha,ipt),"",200000,-100000,100000);
   // h[ipt][ialpha]=new TH1D(Form("h_%d_%d",ialpha,ipt),"",10000,-5000,5000);
   // t[ipt]->Draw(Form("alpha_%d>>h_%d_%d",ialpha,ialpha,ipt),Form("pt1>%d && pt2>%d && abs((eta1+eta2)/2)<%.1f && dphi>(7*TMath::Pi()/8)",jetpt1,jetpt2,etadijet),"");
   if(domp<2)t[ipt]->Draw(Form("-(alpha_%d+alpha_%d)>>h_%d_%d",ialpha,nalpha-ialpha-1,ialpha,ipt),Form("pt1>%d && pt2>%d && abs(eta1)<%.2f && abs(eta2)<%.2f && dphi>(7*TMath::Pi()/8) && pt3<30 && cent<%d && cent>=%d",jetpt1,jetpt2,etadijet,etadijet,cent_max,cent_min),"");
   if(domp==2)t[ipt]->Draw(Form("-(alpha_%d+alpha_%d)>>h_%d_%d",ialpha,nalpha-ialpha-1,ialpha,ipt),Form("pt1>%d && pt2>%d && (abs(eta1+eta2)/2)<%.2f && dphi>(7*TMath::Pi()/8) && cent<%d && cent>=%d",jetpt1,jetpt2,etadijet,cent_max,cent_min),"");
   // t[ipt]->Draw(Form("mpt>>h_%d_%d",ialpha,ipt),"");
   mpt[ipt][ialpha]=h[ipt][ialpha]->GetMean();
   mpterr[ipt][ialpha]=h[ipt][ialpha]->GetMeanError();
   if(ipt<npt-1){
    mptsum[ialpha]+=mpt[ipt][ialpha];
    if(mpt[ipt][ialpha]>=0){
     mptsumpos[ialpha]+=mpt[ipt][ialpha];
     hmpt_vs_alpha[ipt]->SetBinContent(ialpha+1,mptsumpos[ialpha]);
     hmpt_vs_alpha[ipt]->SetBinError(ialpha+1,mpterr[ipt][ialpha]);
    }
    if(mpt[ipt][ialpha]<0){
     // cout<<ialpha<<" "<<mpt[ipt][ialpha]<<endl;
     mptsumneg[ialpha]+=mpt[ipt][ialpha];
     hmpt_vs_alpha[ipt]->SetBinContent(ialpha+1,mptsumneg[ialpha]);
     hmpt_vs_alpha[ipt]->SetBinError(ialpha+1,mpterr[ipt][ialpha]);
    }
	  cout<<ptmin[ipt]<<" "<<ptmax[ ipt]<<" "<<mpt[ipt][ialpha]<<endl;
    hmpt_vs_alpha[ipt]->SetLineColor(1);
    hmpt_vs_alpha[ipt]->SetMarkerColor(1);
    hmpt_vs_alpha[ipt]->SetMarkerSize(0);
    hmpt_vs_alpha[ipt]->SetFillColor(col[ipt]);
   }else{
     hmpt_vs_alpha[ipt]->SetBinContent(ialpha+1,mpt[ipt][ialpha]);
     hmpt_vs_alpha[ipt]->SetBinError(ialpha+1,mpterr[ipt][ialpha]);
   }
  }    
 }
   
 // cout<<"all events"<<endl; 
 for(int ialpha=0;ialpha<nalpha/2;ialpha++){
  cout<<(mptsum[ialpha])<<" "<<mpt[npt-1][ialpha]<<endl;
 }
 
 TCanvas * c1 = new TCanvas("c1","",600,600);
 // TH1D * empty=new TH1D("empty",";#alpha/#pi;<#slash{p}_{T}^{#parallel}>",100,0.025,1.025);
 TH1D * empty=new TH1D("empty",";#alpha^{#pm}/#pi;<#slash{p}_{T}^{#parallel}> (GeV)",nalpha/2,frac);
 empty->Fill(0.5,1000);
 empty->SetMaximum(35); 
 empty->SetMinimum(-60); 
 if(domp==2){
 empty->SetMaximum(85); 
 empty->SetMinimum(-100);  
 }
 
 empty->Draw();
 
 for(int ipt=npt-2;ipt>=0;ipt--){ 
  if(ipt<2)leg->AddEntry(hmpt_vs_alpha[ipt],Form("%.1f - %.1f GeV",ptmin[ipt],ptmax[ipt]),"f l");
  if(ipt>=2)leg2->AddEntry(hmpt_vs_alpha[ipt],Form("%.1f - %.1f GeV",ptmin[ipt],ptmax[ipt]),"f l");
  hmpt_vs_alpha[ipt]->Draw("same");
  hmpt_vs_alpha[ipt]->Draw("same hist");
 }
 leg->AddEntry(hmpt_vs_alpha[npt-1],">0.5 GeV","p");
 hmpt_vs_alpha[npt-1]->Draw("same");
 leg->Draw("same");
 leg2->Draw("same");
 if(!doMC)drawText(Form("PbPb #sqrt{s_{NN}}=2.76 TeV"),0.28,0.23);
 if(doMC) drawText(Form("PYTHIA+HYDJET"),0.28,0.23);
 if(domp<2)drawText(Form("|#eta_{1}|,|#eta_{2}|<%.2f",etadijet),0.28,0.3);
 if(domp==2)drawText(Form("|#eta_{dijet}|<%.2f",etadijet),0.28,0.3);
 drawText(Form("p_{T,1}>%d, p_{T,2}>%d",jetpt1,jetpt2),0.28,0.37);
 drawText(Form("|#eta|<2, #Delta#phi_{1,2}>7#pi/8",jetpt1,jetpt2),0.64,0.23);
 drawText(Form("akPu3Calo"),0.64,0.3);
 drawText(Form("%s",spart[gensigbkgreco].Data()),0.64,0.37);
 drawText(Form("%d-%d %%",(int)(0.5*cent_min),(int)(0.5*cent_max)),0.22,0.93);
 if(!doMC){
  c1->SaveAs(Form("mpt_vs_alpha_%s%s_dijet%d_cent%d_%d.png",smpt[domp].Data(),skind[gensigbkgreco].Data(),(int)etadijet,cent_min,cent_max));
  c1->SaveAs(Form("mpt_vs_alpha_%s%s_dijet%d_cent%d_%d.pdf",smpt[domp].Data(),skind[gensigbkgreco].Data(),(int)etadijet,cent_min,cent_max));
 }
 if(doMC){
  c1->SaveAs(Form("mpt_vs_alpha_%s%s_dijet%d_cent%d_%d_MC.png",smpt[domp].Data(),skind[gensigbkgreco].Data(),(int)etadijet,cent_min,cent_max));
  c1->SaveAs(Form("mpt_vs_alpha_%s%s_dijet%d_cent%d_%d_MC.pdf",smpt[domp].Data(),skind[gensigbkgreco].Data(),(int)etadijet,cent_min,cent_max));
 }
}