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

void drawText(const char *text, float xp, float yp, int size=18){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(43);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw("same");
}

void integrated(int cent_min=0,int cent_max=20,int gensigbkgreco=0,double etadijet=2){
 TH1D::SetDefaultSumw2();
 int domp=0;
 bool doGenJet=true;
 TString smpt[]={"mpt","mp","mpt_boosted"};
 TString skind[]={"","_s","_b","_tracks","_tracks_uncorr","_parts"};
 TString spart[]={"gen. part.","gen. sig.","gen. part.","corr reco. part.","uncorr. reco.","gen. part."};
 int doMC=1;
 int doAve=1;
 TString strave;
 if(doAve) strave="_ave";
 else strave="";
 int docorr=0;
 if(gensigbkgreco==3)docorr=1;
 int jetpt1=120; 
 int jetpt2=50;
 double Aj=0.2;
 int npt=6;
 double ptmin[]={  8,4,2,1,0.5,0.5};
 double ptmax[]={300,8,4,2,  1,300};
 TFile *f[npt];
 TTree *t[npt];
 TTree *t_jet[npt]; 
 int nalpha=33;
 double mpt[npt][nalpha+1];
 double mpterr[npt][nalpha+1]; 
 double mptsum[nalpha+1];
 double mptsumpos[nalpha+1];
 double mptsumneg[nalpha+1];
 TH1D* h[npt][nalpha+1];
 TH1D *hsum[nalpha+1];
 TH1D * hmpt_vs_alpha[npt];
 
 TFile *f_ref[npt];
 TTree *t_ref[npt];
 TTree *t_jet_ref[npt]; 
 double mpt_ref[npt][nalpha+1];
 double mpterr_ref[npt][nalpha+1];
 double mptsum_ref[nalpha+1];
 double mptsumpos_ref[nalpha+1];
 double mptsumneg_ref[nalpha+1];
 TH1D* h_ref[npt][nalpha+1];
 TH1D *hsum_ref[nalpha+1];
 TH1D * hmpt_vs_alpha_ref[npt];
 
 double mpt_diff[npt][nalpha+1];
 double mpterr_diff[npt][nalpha+1];
 double mptsum_diff[nalpha+1];
 double mptsumpos_diff[nalpha+1];
 double mptsumneg_diff[nalpha+1];
 TH1D *hsum_diff[nalpha+1];
 TH1D * hmpt_vs_alpha_diff[npt];
 
 // int col[]={kAzure-9,kYellow-9,kOrange+1,kGreen+3,kPink-2,kRed+1,1};
 int col[]={kRed+1,kGreen+3,kOrange+1,kYellow-9,kAzure-9,1};
 TString sintegrate[nalpha+1][npt];
 TString s[npt];
 for(int ipt=0; ipt<npt;ipt++){
  s[ipt]="";
 }
 for(int ialpha=0;ialpha<nalpha+1;ialpha++){
  mptsum[ialpha]=0;
  mptsumpos[ialpha]=0;
  mptsumneg[ialpha]=0;
  mptsum_ref[ialpha]=0;
  mptsumpos_ref[ialpha]=0;
  mptsumneg_ref[ialpha]=0;
  mptsum_diff[ialpha]=0;
  mptsumpos_diff[ialpha]=0;
  mptsumneg_diff[ialpha]=0;
 }
   
 TLegend *leg,*leg2;
 leg= new TLegend(0.35,0.65,0.56,0.95);
 leg->SetBorderSize(0); 
 leg->SetFillStyle(0);
 leg2 = new TLegend(0.7,0.75,0.98,0.95);
 leg2->SetBorderSize(0); 
 leg2->SetFillStyle(0); 
   
 double frac[]={0,0.01,0.02,0.03,0.04,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,1,1.1,1.25,1.5,2,2.5,3,5,6};
 
 for(int ipt=0;ipt<npt;ipt++){ 
  cout<<1<<endl;
  if(!doMC){  
    cout<<2<<endl;

   f[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_PbPb_mptonly_20140322/full_ntuple_hiForest_Jet80or95_GR_R_53_LV6_12Mar2014_0000CET_Track8_Jet21_0_pt%d_%d_akVs3Calo.root",(int)ptmin[ipt],(int)ptmax[ipt])); 
   f_ref[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_pp_mptonly_20140322/full_ntuple_pt80_pp2013_P01_prod22_v81_merged_forest_0_pt%d_%d_ak3Calo.root",(int)ptmin[ipt],(int)ptmax[ipt])); 
  }
  
  if(doMC){
   f[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_MC_mptonly_20140322/full_ntuple_HiForest_Pythia_Hydjet_Jet80_Track8_Jet21_STARTHI53_LV1_merged_forest_0_pt%d_%d_akVs3Calo_doGenJet%d.root",(int)ptmin[ipt],(int)ptmax[ipt],doGenJet)); 
   f_ref[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_PYTHIA_mptonly_20140321/full_ntuple_pt80_pp2013_P01_prod22_v81_merged_forest_0_pt%d_%d_ak3Calo_doGenJet%d.root",(int)ptmin[ipt],(int)ptmax[ipt],doGenJet)); 
  }
  cout<<3<<endl;

  t[ipt]=(TTree*)f[ipt]->Get(Form("nt_%s%s_dR%s",smpt[domp].Data(),skind[gensigbkgreco].Data(),strave.Data()));
      cout<<3.1<<endl;
  t_jet[ipt]=(TTree*)f[ipt]->Get("nt_jet");
  
  
      cout<<3.2<<endl;

  t[ipt]->AddFriend(t_jet[ipt]);
  
      cout<<3.3<<endl;

  hmpt_vs_alpha[ipt]= new TH1D(Form("hmpt_%d",ipt),"",nalpha+1,frac);
    cout<<4<<endl;

  t_ref[ipt]=(TTree*)f_ref[ipt]->Get(Form("nt_%s%s_dR%s",smpt[domp].Data(),skind[gensigbkgreco].Data(),strave.Data()));
  t_jet_ref[ipt]=(TTree*)f_ref[ipt]->Get("nt_jet");
  t_ref[ipt]->AddFriend(t_jet_ref[ipt]);
  hmpt_vs_alpha_ref[ipt]= new TH1D(Form("hmpt_ref_%d",ipt),"",nalpha+1,frac);
  hmpt_vs_alpha_diff[ipt]= new TH1D(Form("hmpt_diff_%d",ipt),"",nalpha+1,frac);
      cout<<5<<endl;

  for(int ialpha=0;ialpha<nalpha+1;ialpha++ ){
   if(ialpha<nalpha)s[ipt]+=Form("-dR_%d",ialpha);
   // if(ialpha<nalpha)s[ipt]+=Form("-alpha_%d",ialpha);
   cout<<s[ipt].Data()<<endl;
   if(ialpha<nalpha)sintegrate[ialpha][ipt]=s[ipt];
   else  sintegrate[ialpha][ipt]=Form("-%s%s%s",smpt[domp].Data(),skind[gensigbkgreco].Data(),strave.Data());
   h[ipt][ialpha]=new TH1D(Form("h_%d_%d",ialpha,ipt),"",100000,-100000,100000);
   t[ipt]->Draw(Form("%s>>h_%d_%d",sintegrate[ialpha][ipt].Data(),ialpha,ipt),Form("pt1>%d && pt2>%d && abs(eta1)<%.1f && abs(eta2)<%.1f && dphi>(9*TMath::Pi()/10) && cent>=%d && cent<%d && ((pt1-pt2)/(pt2+pt1))>%.2f",jetpt1,jetpt2,etadijet,etadijet, cent_min,cent_max,Aj),"");
   mpt[ipt][ialpha]=h[ipt][ialpha]->GetMean(); 
   mpterr[ipt][ialpha]=h[ipt][ialpha]->GetMeanError();
   
   h_ref[ipt][ialpha]=new TH1D(Form("h_ref_%d_%d",ialpha,ipt),"",100000,-100000,100000);
   t_ref[ipt]->Draw(Form("%s>>h_ref_%d_%d",sintegrate[ialpha][ipt].Data(),ialpha,ipt),Form("pt1>%d && pt2>%d && abs(eta1)<%.1f && abs(eta2)<%.1f && dphi>(9*TMath::Pi()/10)&& ((pt1-pt2)/(pt2+pt1))>%.2f",jetpt1,jetpt2,etadijet,etadijet,Aj),"");
   mpt_ref[ipt][ialpha]=h_ref[ipt][ialpha]->GetMean();
   mpterr_ref[ipt][ialpha]=h_ref[ipt][ialpha]->GetMeanError();
   
   mpt_diff[ipt][ialpha]=mpt[ipt][ialpha]-mpt_ref[ipt][ialpha];
   mpterr_diff[ipt][ialpha]=sqrt(pow(mpterr[ipt][ialpha],2)+pow(mpterr_ref[ipt][ialpha],2));
   // mpterr_diff[ipt][ialpha]=0;
   
   if(ipt<npt-1){
    cout<<mptsum[ialpha]<<endl;
    mptsum[ialpha]+=mpt[ipt][ialpha];
    if(mpt[ipt][ialpha]>=0){
     mptsumpos[ialpha]+=mpt[ipt][ialpha];
     hmpt_vs_alpha[ipt]->SetBinContent(ialpha+1,mptsumpos[ialpha]);
     hmpt_vs_alpha[ipt]->SetBinError(ialpha+1,mpterr[ipt][ialpha]);
    }
    if(mpt[ipt][ialpha]<0){
     mptsumneg[ialpha]+=mpt[ipt][ialpha];
     hmpt_vs_alpha[ipt]->SetBinContent(ialpha+1,mptsumneg[ialpha]);
     hmpt_vs_alpha[ipt]->SetBinError(ialpha+1,mpterr[ipt][ialpha]);
    } 
    hmpt_vs_alpha[ipt]->SetLineColor(1);
    hmpt_vs_alpha[ipt]->SetMarkerColor(1);
    hmpt_vs_alpha[ipt]->SetMarkerSize(0);
    hmpt_vs_alpha[ipt]->SetFillColor(col[ipt]); 
   }else{
     hmpt_vs_alpha[ipt]->SetBinContent(ialpha+1,mpt[ipt][ialpha]);
     hmpt_vs_alpha[ipt]->SetBinError(ialpha+1,mpterr[ipt][ialpha]);
   }   
   
   if(ipt<npt-1){
    mptsum_ref[ialpha]+=mpt_ref[ipt][ialpha];
    if(mpt_ref[ipt][ialpha]>=0){
     mptsumpos_ref[ialpha]+=mpt_ref[ipt][ialpha];
     hmpt_vs_alpha_ref[ipt]->SetBinContent(ialpha+1,mptsumpos_ref[ialpha]);
     hmpt_vs_alpha_ref[ipt]->SetBinError(ialpha+1,mpterr_ref[ipt][ialpha]);
    }
    if(mpt_ref[ipt][ialpha]<0){  
     cout<<ialpha<<" "<<mpt_ref[ipt][ialpha]<<endl;
     mptsumneg_ref[ialpha]+=mpt_ref[ipt][ialpha];
     hmpt_vs_alpha_ref[ipt]->SetBinContent(ialpha+1,mptsumneg_ref[ialpha]);
     hmpt_vs_alpha_ref[ipt]->SetBinError(ialpha+1,mpterr_ref[ipt][ialpha]);
    }
    hmpt_vs_alpha_ref[ipt]->SetLineColor(1);
    hmpt_vs_alpha_ref[ipt]->SetMarkerColor(1);
    hmpt_vs_alpha_ref[ipt]->SetMarkerSize(0);
    hmpt_vs_alpha_ref[ipt]->SetFillColor(col[ipt]);
   }else{ 
     hmpt_vs_alpha_ref[ipt]->SetBinContent(ialpha+1,mpt_ref[ipt][ialpha]);
     hmpt_vs_alpha_ref[ipt]->SetBinError(ialpha+1,mpterr_ref[ipt][ialpha]);
   }
   
   if(ipt<npt-1){
    mptsum_diff[ialpha]+=mpt_diff[ipt][ialpha];
    if(mpt_diff[ipt][ialpha]>=0){
     mptsumpos_diff[ialpha]+=mpt_diff[ipt][ialpha];
     hmpt_vs_alpha_diff[ipt]->SetBinContent(ialpha+1,mptsumpos_diff[ialpha]);
     hmpt_vs_alpha_diff[ipt]->SetBinError(ialpha+1,mpterr_diff[ipt][ialpha]);
    }
    if(mpt_diff[ipt][ialpha]<0){  
     mptsumneg_diff[ialpha]+=mpt_diff[ipt][ialpha];
     hmpt_vs_alpha_diff[ipt]->SetBinContent(ialpha+1,mptsumneg_diff[ialpha]);
     hmpt_vs_alpha_diff[ipt]->SetBinError(ialpha+1,mpterr_diff[ipt][ialpha]);
    }
    hmpt_vs_alpha_diff[ipt]->SetLineColor(1);
    hmpt_vs_alpha_diff[ipt]->SetMarkerColor(1);
    hmpt_vs_alpha_diff[ipt]->SetMarkerSize(0);
    hmpt_vs_alpha_diff[ipt]->SetFillColor(col[ipt]);
   }else{ 
     hmpt_vs_alpha_diff[ipt]->SetBinContent(ialpha+1,mpt_diff[ipt][ialpha]);
     hmpt_vs_alpha_diff[ipt]->SetBinError(ialpha+1,mpterr_diff[ipt][ialpha]);
   }
  }
 }
 
 
 // cout<<"all events"<<endl;  
 for(int ialpha=0;ialpha<nalpha+1;ialpha++){
  cout<<(mptsum[ialpha])<<" "<<mpt[npt-1][ialpha]<<endl;
  cout<<(mptsum_ref[ialpha])<<" "<<mpt_ref[npt-1][ialpha]<<endl;
 }

 TCanvas * c1 = new TCanvas("c1","",600,600);
 TH1D * empty=new TH1D("empty",";#DeltaR;<#slash{p}_{T}^{#parallel}> (GeV)",nalpha+1,frac);
 empty->Fill(0.5,1000); 
 empty->SetMaximum(40); 
 empty->SetMinimum(-30); 
 empty->Draw();
  
 for(int ipt=npt-2;ipt>=0;ipt--){ 
  leg2->AddEntry(hmpt_vs_alpha[ipt],Form("%.1f - %.1f GeV",ptmin[ipt],ptmax[ipt]),"f l");
  hmpt_vs_alpha[ipt]->Draw("same");
  hmpt_vs_alpha[ipt]->Draw("same hist");
 }
 leg2->AddEntry(hmpt_vs_alpha[npt-1],">0.5 GeV","p");
 hmpt_vs_alpha[npt-1]->Draw("same");
 TLine l(0.025,0,0.525,0);
 l.SetLineWidth(2);
 // l.Draw("same");
 // leg->Draw("same");
 leg2->Draw("same"); 
 if(!doMC)drawText(Form("PbPb"),0.22,0.87);
 if(doMC) drawText(Form("(PYT.+HYD.)"),0.22,0.87);
 if(domp<2)drawText(Form("|#eta_{1}|,|#eta_{2}|<%.2f",etadijet),0.22,0.81);
 if(domp==2)drawText(Form("|#eta_{dijet}|<%.2f",etadijet),0.22,0.81);
 drawText(Form("p_{T,1}>%d, p_{T,2}>%d",jetpt1,jetpt2),0.22,0.75);
 drawText(Form("|#eta|<2, #Delta#phi_{1,2}>9#pi/10",jetpt1,jetpt2),0.44,0.93);
 drawText(Form("akVs3Calo, ak3Calo"),0.44,0.87);
 drawText(Form("%s",spart[gensigbkgreco].Data()),0.44,0.81);
 drawText(Form("%d-%d %%",(int)(0.5*cent_min),(int)(0.5*cent_max)),0.22,0.93);
 if(!doMC){
  c1->SaveAs(Form("plots_integrated_data/integrated_mpt_vs_alpha_%s%s_dijet%d_cent%d_%d_doAve%d_dR.png",smpt[domp].Data(),skind[gensigbkgreco].Data(),(int)etadijet,cent_min,cent_max,doAve));
  c1->SaveAs(Form("plots_integrated_data/integrated_mpt_vs_alpha_%s%s_dijet%d_cent%d_%d_doAve%d_dR.pdf",smpt[domp].Data(),skind[gensigbkgreco].Data(),(int)etadijet,cent_min,cent_max,doAve));
 }
 if(doMC){
  c1->SaveAs(Form("plots_integrated_mc%d/integrated_mpt_vs_alpha_%s%s_dijet%d_cent%d_%d_doAve%d_dR_MC.png",doGenJet,smpt[domp].Data(),skind[gensigbkgreco].Data(),(int)etadijet,cent_min,cent_max,doAve));
  c1->SaveAs(Form("plots_integrated_mc%d/integrated_mpt_vs_alpha_%s%s_dijet%d_cent%d_%d_doAve%d_dR_MC.pdf",doGenJet,smpt[domp].Data(),skind[gensigbkgreco].Data(),(int)etadijet,cent_min,cent_max,doAve));
 } 
 
 TCanvas * c2 = new TCanvas("c2","",600,600);
 empty->Draw();
 
 for(int ipt=npt-2;ipt>=0;ipt--){ 
  hmpt_vs_alpha_ref[ipt]->Draw("same");
  hmpt_vs_alpha_ref[ipt]->Draw("same hist");
 }
 hmpt_vs_alpha_ref[npt-1]->Draw("same");
 l.SetLineWidth(2);
 leg2->Draw("same");
 leg2->Draw("same"); 
 if(!doMC)drawText(Form("pp"),0.22,0.87);
 if(doMC) drawText(Form("PYT."),0.22,0.87);
 if(domp<2)drawText(Form("|#eta_{1}|,|#eta_{2}|<%.2f",etadijet),0.22,0.81);
 if(domp==2)drawText(Form("|#eta_{dijet}|<%.2f",etadijet),0.22,0.81);
 drawText(Form("p_{T,1}>%d, p_{T,2}>%d",jetpt1,jetpt2),0.22,0.75);
 drawText(Form("|#eta|<2, #Delta#phi_{1,2}>9#pi/10",jetpt1,jetpt2),0.44,0.93);
 drawText(Form("akVs3Calo, ak3Calo"),0.44,0.87);
 drawText(Form("%s",spart[gensigbkgreco].Data()),0.44,0.81);
 drawText(Form("%d-%d %%",(int)(0.5*cent_min),(int)(0.5*cent_max)),0.22,0.93);
 // drawText(Form("%d-%d %%",(int)(0.5*cent_min),(int)(0.5*cent_max)),0.22,0.93);
 if(!doMC){
  c2->SaveAs(Form("plots_integrated_data/integrated_mpt_vs_alpha_%s%s_dijet%d_cent%d_%d_doAve%d_dR_ref.png",smpt[domp].Data(),skind[gensigbkgreco].Data(),(int)etadijet,cent_min,cent_max,doAve));
  c2->SaveAs(Form("plots_integrated_data/integrated_mpt_vs_alpha_%s%s_dijet%d_cent%d_%d_doAve%d_dR_ref.pdf",smpt[domp].Data(),skind[gensigbkgreco].Data(),(int)etadijet,cent_min,cent_max,doAve));
 }
 if(doMC){
  c2->SaveAs(Form("plots_integrated_mc%d/integrated_mpt_vs_alpha_%s%s_dijet%d_cent%d_%d_doAve%d_dR_MC_ref.png",doGenJet,smpt[domp].Data(),skind[gensigbkgreco].Data(),(int)etadijet,cent_min,cent_max,doAve));
  c2->SaveAs(Form("plots_integrated_mc%d/integrated_mpt_vs_alpha_%s%s_dijet%d_cent%d_%d_doAve%d_dR_MC_ref.pdf",doGenJet,smpt[domp].Data(),skind[gensigbkgreco].Data(),(int)etadijet,cent_min,cent_max,doAve));
 }
 
 TCanvas * c3 = new TCanvas("c3","",600,600);
 empty->Draw();
 
 for(int ipt=npt-2;ipt>=0;ipt--){ 
  hmpt_vs_alpha_diff[ipt]->Draw("same");
  hmpt_vs_alpha_diff[ipt]->Draw("same hist");
 }
 hmpt_vs_alpha_diff[npt-1]->Draw("same"); 
 l.SetLineWidth(2);
 leg2->Draw("same"); 
 leg2->Draw("same");  
 if(!doMC)drawText(Form("PbPb-pp"),0.22,0.87);
 if(doMC) drawText(Form("(PYT.+HYD.)-PYT."),0.22,0.87);
 if(domp<2)drawText(Form("|#eta_{1}|,|#eta_{2}|<%.2f",etadijet),0.22,0.81);
 if(domp==2)drawText(Form("|#eta_{dijet}|<%.2f",etadijet),0.22,0.81);
 drawText(Form("p_{T,1}>%d, p_{T,2}>%d",jetpt1,jetpt2),0.22,0.75);
 drawText(Form("|#eta|<2, #Delta#phi_{1,2}>9#pi/10",jetpt1,jetpt2),0.44,0.93);
 drawText(Form("akVs3Calo, ak3Calo"),0.44,0.87);
 drawText(Form("%s",spart[gensigbkgreco].Data()),0.44,0.81);
 drawText(Form("%d-%d %%",(int)(0.5*cent_min),(int)(0.5*cent_max)),0.22,0.93);
 if(!doMC){
  c3->SaveAs(Form("plots_integrated_data/integrated_mpt_vs_alpha_%s%s_dijet%d_cent%d_%d_doAve%d__dR_diff.png",smpt[domp].Data(),skind[gensigbkgreco].Data(),(int)etadijet,cent_min,cent_max,doAve));
  c3->SaveAs(Form("plots_integrated_data/integrated_mpt_vs_alpha_%s%s_dijet%d_cent%d_%d_doAve%d_dR_diff.pdf",smpt[domp].Data(),skind[gensigbkgreco].Data(),(int)etadijet,cent_min,cent_max,doAve));
 }
 if(doMC){
  c3->SaveAs(Form("plots_integrated_mc%d/integrated_mpt_vs_alpha_%s%s_dijet%d_cent%d_%d_doAve%d_dR_MC_diff.png",doGenJet,smpt[domp].Data(),skind[gensigbkgreco].Data(),(int)etadijet,cent_min,cent_max,doAve));
  c3->SaveAs(Form("plots_integrated_mc%d/integrated_mpt_vs_alpha_%s%s_dijet%d_cent%d_%d_doAve%d_dR_MC_diff.pdf",doGenJet,smpt[domp].Data(),skind[gensigbkgreco].Data(),(int)etadijet,cent_min,cent_max,doAve));
 }
}