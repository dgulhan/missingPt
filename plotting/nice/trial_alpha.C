#include <iostream>
#include "TCanvas.h"
#include "TError.h"
#include "TPad.h"
#include "TString.h"
#include "TRandom.h"
#include "TH1F.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TString.h"
#include "TCut.h"
#include "TNtuple.h"
#include "TLine.h"



void makeMultiPanelCanvas(TCanvas*& canv,
                          const Int_t columns,
                          const Int_t rows,
                          const Float_t leftOffset,
                          const Float_t bottomOffset,
                          const Float_t leftMargin,
                          const Float_t bottomMargin,
                          const Float_t edge) {
   if (canv==0) {
      Error("makeMultiPanelCanvas","Got null canvas.");
      return;
   }
   canv->Clear();
   
   TPad* pad[columns][rows];

   Float_t Xlow[columns];
   Float_t Xup[columns];
   Float_t Ylow[rows];
   Float_t Yup[rows];
   Float_t PadWidth = 
   (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
   (1.0/(1.0-edge))+(Float_t)columns-2.0);
   Float_t PadHeight =
   (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
   (1.0/(1.0-edge))+(Float_t)rows-2.0);
   Xlow[0] = leftOffset;
   Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
   Xup[columns-1] = 1;
   Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);

   Yup[0] = 1;
   Ylow[0] = 1.0-PadHeight/(1.0-edge);
   Ylow[rows-1] = bottomOffset;
   Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);

   for(Int_t i=1;i<columns-1;i++) {
      Xlow[i] = Xup[0] + (i-1)*PadWidth;
      Xup[i] = Xup[0] + (i)*PadWidth;
   }
   Int_t ct = 0;
   for(Int_t i=rows-2;i>0;i--) {
      Ylow[i] = Yup[rows-1] + ct*PadHeight;
      Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
      ct++;
   }

   TString padName;
   for(Int_t i=0;i<columns;i++) {
      for(Int_t j=0;j<rows;j++) {
         canv->cd();
         padName = Form("p_%d_%d",i,j);
         pad[i][j] = new TPad(padName.Data(),padName.Data(),
            Xlow[i],Ylow[j],Xup[i],Yup[j]);
         if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
         else pad[i][j]->SetLeftMargin(0);

         if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
         else pad[i][j]->SetRightMargin(0);

         if(j==0) pad[i][j]->SetTopMargin(edge);
         else pad[i][j]->SetTopMargin(0);

         if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
         else pad[i][j]->SetBottomMargin(0);

         pad[i][j]->Draw();
         pad[i][j]->cd();
         pad[i][j]->SetNumber(columns*j+i+1);
      }
   }
}

void drawText(const char *text, float xp, float yp, int size=18){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(43);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw("same");
}

void trial(int gensigbkgreco=3){
TH1D::SetDefaultSumw2();
  TH1D::SetDefaultSumw2();
 int domp=0;
 bool doGenJet=true;
 TString smpt[]={"mpt","mp","mpt_boosted"};
 TString skind[]={"","_s","_b","_tracks","_tracks_uncorr","_parts"};
 TString spart[]={"gen. part.","gen. sig.","gen. part.","corr reco. part.","uncorr. reco.","gen. part."};
 int doMC=1;
 int doAve=0;
 TString strave;
 if(doAve) strave="_ave";
 else strave="";
 int docorr=0;
 if(gensigbkgreco==3)docorr=1;
 int jetpt1=120; 
 int jetpt2=50;
 double Aj=0.0;
 int npt=6;
 double ptmin[]={  8,4,2,1,0.5,0.5};
 double ptmax[]={300,8,4,2,  1,300};
 TFile *f[npt];
 TTree *t[npt];
 TTree *t_jet[npt]; 
 int nalpha=30;
 double etadijet=0.5;
 TFile *f_ref[npt];
 TTree *t_ref[npt];
 TTree *t_jet_ref[npt]; 
 double mpt_ref[npt][nalpha/2];
 double mpterr_ref[npt][nalpha/2];
 double mptsum_ref[nalpha/2];
 double mptsumpos_ref[nalpha/2];
 double mptsumneg_ref[nalpha/2];
 TH1D* h_ref[npt][nalpha/2];
 TH1D *hsum_ref[nalpha/2];
 TH1D * hmpt_vs_alpha_ref[npt];
 
 int col[]={kRed+1,kGreen+3,kOrange+1,kYellow-9,kAzure-9,1};
 TString sintegrate[nalpha/2][npt];
 TString s[npt];
 for(int ipt=0; ipt<npt;ipt++){
  s[ipt]="";
 }
 for(int ialpha=0;ialpha<nalpha/2;ialpha++){
  mptsum_ref[ialpha]=0;
  mptsumpos_ref[ialpha]=0;
  mptsumneg_ref[ialpha]=0;
 }
   
 TLegend *leg,*leg2;
 leg= new TLegend(0.35,0.65,0.56,0.95);
 leg->SetBorderSize(0); 
 leg->SetFillStyle(0);
 leg2 = new TLegend(0.7,0.75,0.98,0.95);
 leg2->SetBorderSize(0); 
 leg2->SetFillStyle(0); 
   
 double frac[]={0,0.01,0.02,0.03,0.04,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.925,0.95,0.96,0.97,0.98,0.99,1.01};
  

  if(!doMC){  
    f_ref[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_pp_mptonly_20140323_v3/full_ntuple_pt80_pp2013_P01_prod22_v81_merged_forest_0_pt%d_%d_ak3Calo.root",(int)ptmin[ipt],(int)ptmax[ipt])); 
  }
  
  if(doMC){

  f_ref[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_PYTHIA_mptonly_20140323_v3/full_ntuple_pt80_pp2013_P01_prod22_v81_merged_forest_0_pt%d_%d_ak3Calo_doGenJet%d.root",(int)ptmin[ipt],(int)ptmax[ipt],doGenJet)); 
  }
   
   
   cout<<Form("nt_%s%s_dphi%s",smpt[domp].Data(),skind[gensigbkgreco].Data(),strave.Data())<<endl;
  t_ref[ipt]=(TTree*)f_ref[ipt]->Get(Form("nt_%s%s%s",smpt[domp].Data(),skind[gensigbkgreco].Data(),strave.Data()));
  t_jet_ref[ipt]=(TTree*)f_ref[ipt]->Get("nt_jet");
  t_ref[ipt]->AddFriend(t_jet_ref[ipt]);
  hmpt_vs_alpha_ref[ipt]= new TH1D(Form("hmpt_ref_%d",ipt),"",((int)(nalpha/2)),frac);
      cout<<5<<endl;

  for(int ialpha=0;ialpha<nalpha/2;ialpha++ ){
   if(ialpha<nalpha)s[ipt]+=Form("-alpha_%d-alpha_%d",ialpha,nalpha-ialpha-1);
   sintegrate[ialpha][ipt]=s[ipt];
   h_ref[ipt][ialpha]=new TH1D(Form("h_ref_%d_%d",ialpha,ipt),"",100000,-100000,100000);
   t_ref[ipt]->Draw(Form("%s>>h_ref_%d_%d",sintegrate[ialpha][ipt].Data(),ialpha,ipt),Form("pt1>%d && pt2>%d && abs(eta1)<%.1f && abs(eta2)<%.1f && dphi>(5*TMath::Pi()/6)&& ((pt1-pt2)/(pt2+pt1))>%.2f",jetpt1,jetpt2,etadijet,etadijet,Aj),"");
   mpt_ref[ipt][ialpha]=h_ref[ipt][ialpha]->GetMean();
   mpterr_ref[ipt][ialpha]=h_ref[ipt][ialpha]->GetMeanError();
   
   if(ipt<npt-1){
    mptsum_ref[ialpha]+=mpt_ref[ipt][ialpha];
    if(ipt==4)cout<<mptsum_ref[ialpha]<<endl;
    if(mpt_ref[ipt][ialpha]>=0){
	 cout<<Form("-alpha_%d-alpha_%d",ialpha,nalpha-ialpha-1)<<" positive "<<ipt<<" " <<mpt_ref[ipt][ialpha]<<endl;

     mptsumpos_ref[ialpha]+=mpt_ref[ipt][ialpha];
     hmpt_vs_alpha_ref[ipt]->SetBinContent(ialpha,mptsumpos_ref[ialpha]);
     hmpt_vs_alpha_ref[ipt]->SetBinError(ialpha,mpterr_ref[ipt][ialpha]);
    }
    if(mpt_ref[ipt][ialpha]<0){  
	 cout<<Form("-alpha_%d-alpha_%d",ialpha,nalpha-ialpha-1)<<" negative "<<ipt<<" " <<mpt_ref[ipt][ialpha]<<endl;
     mptsumneg_ref[ialpha]+=mpt_ref[ipt][ialpha];
     hmpt_vs_alpha_ref[ipt]->SetBinContent(ialpha,mptsumneg_ref[ialpha]);
     hmpt_vs_alpha_ref[ipt]->SetBinError(ialpha,mpterr_ref[ipt][ialpha]);
    }
    hmpt_vs_alpha_ref[ipt]->SetLineColor(1);
    hmpt_vs_alpha_ref[ipt]->SetMarkerColor(1);
    hmpt_vs_alpha_ref[ipt]->SetMarkerSize(0);
    hmpt_vs_alpha_ref[ipt]->SetFillColor(col[ipt]);
   }else{ 
     hmpt_vs_alpha_ref[ipt]->SetBinContent(ialpha,mpt_ref[ipt][ialpha]);
     hmpt_vs_alpha_ref[ipt]->SetBinError(ialpha,mpterr_ref[ipt][ialpha]);
   }
   
  }
 }
 
 
 TCanvas * c1 = new TCanvas("c1","",600,600);
 // TH1D * empty=new TH1D("empty",";#DeltaR;<#slash{p}_{T}^{#parallel}> (GeV)",nalpha/2,frac);
 TH1D * empty=new TH1D("empty",";#Delta#eta;<#slash{p}_{T}^{#parallel}> (GeV)",((int)(nalpha/2)),frac);
 empty->Fill(0.5,1000); 
 // empty->SetMaximum(30); 
 // empty->SetMinimum(-30); 
 empty->SetMaximum(100); 
 empty->SetMinimum(-100); 
 empty->Draw();
    
 for(int ipt=npt-2;ipt>=0;ipt--){ 
  hmpt_vs_alpha_ref[ipt]->Draw("same");
  hmpt_vs_alpha_ref[ipt]->Draw("same hist");
 }
 hmpt_vs_alpha_ref[npt-1]->Draw("same"); 
 // leg2->Draw("same"); 
 // leg2->Draw("same");  
 // if(!doMC)drawText(Form("PbPb-pp"),0.22,0.87);
 // if(doMC) drawText(Form("(PYT.+HYD.)-PYT."),0.22,0.87);
 drawText(Form("PYTHIA",etadijet),0.22,0.93);
 drawText(Form("p_{T,1}>%d, p_{T,2}>%d",jetpt1,jetpt2),0.22,0.87);
  if(domp<2)drawText(Form("|#eta_{1}|,|#eta_{2}|<%.2f",etadijet),0.22,0.81);
 drawText(Form("|#eta|<2, #Delta#phi_{1,2}>5#pi/6",jetpt1,jetpt2),0.44,0.93);
 if(!doGenJet)drawText(Form("akVs3Calo, ak3Calo"),0.44,0.87);
 else drawText(Form("gen jet"),0.44,0.87);
 drawText(Form("%s",spart[gensigbkgreco].Data()),0.44,0.81);
 drawText(Form("%s",spart[gensigbkgreco].Data()),0.44,0.81);
 c1->SaveAs("trial_alpha_midrapidity_diff.png");
}