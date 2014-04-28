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

void incone(){
 TH1D::SetDefaultSumw2();
 int domp=0;
 int gensigbkgreco=1;
 int cent_min=0;
 int cent_max=60;
 double etadijet=1.6; 
 int doMC=1; 
 int docorr=1;
 if(gensigbkgreco!=3) docorr=0; 
 int jetpt1=120;
 int jetpt2=50;
 int doGenJet=1;
 TString smpt[]={"mpt","mp","mpt_boosted"};
 TString skind[]={"","_s","_b","_tracks","tracks_uncorr"};
 TString spart[]={"gen. part.","gen. sig.","gen. part.","corr reco. part.","uncorr. reco."};
 TString weight_MC;
 // if(doMC) weight_MC="cent_weight";
 if(doMC) weight_MC="1";
 else weight_MC="1";
 int npt=6;
 double ptmin[]={  8,4,2,1,0.5,0.5};
 double ptmax[]={300,8,4,2,  1,300};
 // double ptmin[]={0.5,1,2,4,8  ,0.5};
 // double ptmax[]={1  ,2,4,8,100,100};
 TFile *f[npt];
 TTree *t[npt];
 TTree *t_jet[npt];
 int nAj=4;
 double mpt[npt][nAj];
 double mpterr[npt][nAj];
 double mptsum[nAj];
 double mptsumpos[nAj];
 double mptsumneg[nAj];
 TH1D* h[npt][nAj];
 TH1D *hsum[nAj];
 TH1D * hmpt_vs_Aj[npt];
 int col[]={kRed+1,kGreen+3,kOrange+1,kYellow-9,kAzure-9,1};

 for(int iAj=0;iAj<nAj;iAj++){
  mptsum[iAj]=0;
  mptsumpos[iAj]=0;
  mptsumneg[iAj]=0;
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
 
  double Aj[]={0,0.11,0.22,0.33,0.5};
  
 for(int ipt=0;ipt<npt;ipt++){  
  if(!doMC){  
   f[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_PbPb_mptonly_20140325_v5/full_ntuple_hiForest_Jet80or95_GR_R_53_LV6_03Mar2014_1600CET_CMSSW_5_3_16_merged_pt%d_%d_akVs3Calo.root",(int)ptmin[ipt],(int)ptmax[ipt])); 
   // f_ref[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_pp_mptonly_20140325_v4/full_ntuple_PP2013_HiForest_PromptReco_JsonPP_Jet80_PPReco_forestv82_pt%d_%d_ak3Calo.root",(int)ptmin[ipt],(int)ptmax[ipt])); 
  }
  
  if(doMC){
   f[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_MC_mptonly_20140325_pthat120/full_ntuple_HydjetDrum_Pyquen_Dijet120_FOREST_Track8_Jet24_FixedPtHat_v0_pt%d_%d_akVs3Calo_doGenJet%d.root",(int)ptmin[ipt],(int)ptmax[ipt],doGenJet)); 
   // f[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_PYTHIA_HIRECO_20140327_v5/full_ntuple_hiForest_QCDpT120_STARTHI53_LV1_Track8_Jet22_1GeVcut_pt%d_%d_ak3Calo_doGenJet%d.root",(int)ptmin[ipt],(int)ptmax[ipt],doGenJet)); 
   // f[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_PYTHIA_mptonly_20140325_v5/full_ntuple_pt120_pp2013_P01_prod22_v81_merged_forest_0_pt%d_%d_ak3Calo_doGenJet%d.root",(int)ptmin[ipt],(int)ptmax[ipt],doGenJet)); 
  }
  t[ipt]=(TTree*)f[ipt]->Get(Form("nt_%s%s_dR_ave",smpt[domp].Data(),skind[gensigbkgreco].Data()));
  t_jet[ipt]=(TTree*)f[ipt]->Get("nt_jet");
  t[ipt]->AddFriend(t_jet[ipt]);
  hmpt_vs_Aj[ipt]= new TH1D(Form("hmpt_%d",ipt),"",nAj,Aj);
  for(int iAj=0;iAj<nAj;iAj++ ){
   h[ipt][iAj]=new TH1D(Form("h_%d_%d",iAj,ipt),"",200000,-10000,10000);
   // if(iAj<nAj-1)t[ipt]->Draw(Form("-(dR_8)>>h_%d_%d",iAj,ipt),Form("%s*(pt1>%d && pt2>%d && abs(eta1)<%.2f && abs(eta2)<%.2f && dphi>(5*TMath::Pi()/6) && cent<%d && cent>=%d && ((pt1-pt2)/(pt1+pt2))>=%.2f && ((pt1-pt2)/(pt1+pt2))<%.2f)",weight_MC.Data(),jetpt1,jetpt2,etadijet,etadijet,cent_max,cent_min,Aj[iAj],Aj[iAj+1]),"");
   // else t[ipt]->Draw(Form("-(dR_8)>>h_%d_%d",iAj,ipt),Form("%s*(pt1>%d && pt2>%d && abs(eta1)<%.2f && abs(eta2)<%.2f && dphi>(5*TMath::Pi()/6) && cent<%d && cent>=%d && ((pt1-pt2)/(pt1+pt2))>=%.2f)",weight_MC.Data(),jetpt1,jetpt2,etadijet,etadijet,cent_max,cent_min,Aj[iAj]),"");
   if(iAj<nAj-1)t[ipt]->Draw(Form("-(dR_7)>>h_%d_%d",iAj,ipt),Form("%s*(pt1>%d && pt2>%d && abs(eta1)<%.2f && abs(eta2)<%.2f && dphi>(5*TMath::Pi()/6) && ((pt1-pt2)/(pt1+pt2))>=%.2f && ((pt1-pt2)/(pt1+pt2))<%.2f && abs(vz)<15)",weight_MC.Data(),jetpt1,jetpt2,etadijet,etadijet,Aj[iAj],Aj[iAj+1]),"");
   else t[ipt]->Draw(Form("-(dR_7)>>h_%d_%d",iAj,ipt),Form("%s*(pt1>%d && pt2>%d && abs(eta1)<%.2f && abs(eta2)<%.2f && dphi>(5*TMath::Pi()/6)&& ((pt1-pt2)/(pt1+pt2))>=%.2f && abs(vz)<15)",weight_MC.Data(),jetpt1,jetpt2,etadijet,etadijet,Aj[iAj]),"");
   mpt[ipt][iAj]=h[ipt][iAj]->GetMean();
   mpterr[ipt][iAj]=h[ipt][iAj]->GetMeanError();
   if(ipt<npt-1){
    mptsum[iAj]+=mpt[ipt][iAj];
    if(mpt[ipt][iAj]>=0){
     mptsumpos[iAj]+=mpt[ipt][iAj];
     hmpt_vs_Aj[ipt]->SetBinContent(iAj+1,mptsumpos[iAj]);
     hmpt_vs_Aj[ipt]->SetBinError(iAj+1,mpterr[ipt][iAj]);
    }
    if(mpt[ipt][iAj]<0){
     mptsumneg[iAj]+=mpt[ipt][iAj];
     hmpt_vs_Aj[ipt]->SetBinContent(iAj+1,mptsumneg[iAj]);
     hmpt_vs_Aj[ipt]->SetBinError(iAj+1,mpterr[ipt][iAj]);
    }
	  cout<<ptmin[ipt]<<" "<<ptmax[ ipt]<<" "<<mpt[ipt][iAj]<<endl;
    hmpt_vs_Aj[ipt]->SetLineColor(1);
    hmpt_vs_Aj[ipt]->SetMarkerColor(1);
    hmpt_vs_Aj[ipt]->SetMarkerSize(0);
    hmpt_vs_Aj[ipt]->SetFillColor(col[ipt]);
   }else{
     hmpt_vs_Aj[ipt]->SetBinContent(iAj+1,mpt[ipt][iAj]);
     hmpt_vs_Aj[ipt]->SetBinError(iAj+1,mpterr[ipt][iAj]);
   }
  }    
 }
   
 // cout<<"all events"<<endl; 
 for(int iAj=0;iAj<nAj;iAj++){
  cout<<(mptsum[iAj])<<" "<<mpt[npt-1][iAj]<<endl;
 }
 
 TCanvas * c1 = new TCanvas("c1","",600,600);
 TH1D * empty=new TH1D("empty",";A_{J};<#slash{p}_{T}^{#parallel}> (GeV)",nAj,Aj);
 empty->Fill(0.5,1000);
 empty->SetMaximum(60); 
 empty->SetMinimum(-60); 
 if(domp==2){
 empty->SetMaximum(85); 
 empty->SetMinimum(-100);  
 }
 
 empty->Draw();
 
 for(int ipt=npt-2;ipt>=0;ipt--){ 
  if(ipt<2)leg->AddEntry(hmpt_vs_Aj[ipt],Form("%.1f - %.1f GeV",ptmin[ipt],ptmax[ipt]),"f l");
  if(ipt>=2)leg2->AddEntry(hmpt_vs_Aj[ipt],Form("%.1f - %.1f GeV",ptmin[ipt],ptmax[ipt]),"f l");
  hmpt_vs_Aj[ipt]->Draw("same");
  hmpt_vs_Aj[ipt]->Draw("same hist");
 }
 leg->AddEntry(hmpt_vs_Aj[npt-1],">0.5 GeV","p");
 hmpt_vs_Aj[npt-1]->Draw("same");
 // leg->Draw("same");
 // leg2->Draw("same");
 if(!doMC)drawText(Form("PbPb #sqrt{s_{NN}}=2.76 TeV"),0.28,0.23);
 // if(doMC) drawText(Form("PYTHIA+HYDJET"),0.28,0.23);
 if(doMC) drawText(Form("PYTHIA"),0.28,0.23);
 // if(domp<2)drawText(Form("|#eta_{1}|,|#eta_{2}|<%.2f",etadijet),0.28,0.3);
 // if(domp==2)drawText(Form("|#eta_{dijet}|<%.2f",etadijet),0.28,0.3);
 drawText(Form("p_{T,1}>%d, p_{T,2}>%d",jetpt1,jetpt2),0.28,0.37);
 drawText(Form("|#eta|<2, #Delta#phi_{1,2}>5#pi/6",jetpt1,jetpt2),0.64,0.23);
 // drawText(Form("akVs3Calo"),0.64,0.3);
 drawText(Form("gen jet"),0.64,0.3);
 drawText(Form("%s",spart[gensigbkgreco].Data()),0.28,0.3);
 drawText(Form("%d-%d %%",(int)(0.5*cent_min),(int)(0.5*cent_max)),0.22,0.93);
 if(!doMC){
  c1->SaveAs(Form("incone_mpt_vs_Aj_%s%s_dijet%d_cent%d_%d_gen_gen_ph.png",smpt[domp].Data(),skind[gensigbkgreco].Data(),(int)etadijet,cent_min,cent_max));
  c1->SaveAs(Form("incone_mpt_vs_Aj_%s%s_dijet%d_cent%d_%d_gen_gen_ph.pdf",smpt[domp].Data(),skind[gensigbkgreco].Data(),(int)etadijet,cent_min,cent_max));
 }
 if(doMC){
  c1->SaveAs(Form("incone_mpt_vs_Aj_%s%s_dijet%d_cent%d_%d_MC_gen_gen_ph.png",smpt[domp].Data(),skind[gensigbkgreco].Data(),(int)etadijet,cent_min,cent_max));
  c1->SaveAs(Form("incone_mpt_vs_Aj_%s%s_dijet%d_cent%d_%d_MC_gen_gen_ph.pdf",smpt[domp].Data(),skind[gensigbkgreco].Data(),(int)etadijet,cent_min,cent_max));
 }
}