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

void multiIntegrated(int gensigbkgreco=3, int doMC=0,double etadijet=0.5, int index_var=2, bool doIntegrate=true, int doAve=1){
 TH1D::SetDefaultSumw2(); 
 int domp=0;
 int ncent=2;  
 bool doGenJet=false;
 TString smpt[]={"mpt","mp","mpt_boosted"};
 TString skind[]={"","_s","_b","_tracks","_tracks_uncorr","_parts"};
 TString spart[]={"gen. part.","gen. sig.","gen. part.","corr reco. part.","uncorr. reco.","gen. part."};
 TString svariable[]={"alpha","deta","dphi","dR"};
 TString axistitle[]={"#alpha","#Delta#eta","#Delta#phi","#DeltaR"};
 TString strave;
 if(doAve) strave="_ave";
 else strave=""; 
 int docorr=0;
 if(gensigbkgreco==3)docorr=1;
 int jetpt1=120; 
 int jetpt2=50;
 double Aj=0.0;
 int npt=6;
 // int  cent_min[]={0,20,40,60,100};
 // int cent_max[]={20,40,60,100,200};
 int  cent_min[]={0,60};
 int cent_max[]={60,200};
 double ptmin[]={  8,4,2,1,0.5,0.5};
 double ptmax[]={300,8,4,2,  1,300};
 
 TFile *f[npt];
 TTree *t[npt];
 TTree *t_jet[npt];  
 int nalpha;
 if(index_var==0 ||index_var==2) nalpha=16;
 else nalpha=16;
 double mpt;
 double mpterr; 
 double mptsum[nalpha+1][ncent];
 double mptsumpos[nalpha+1][ncent];
 double mptsumneg[nalpha+1][ncent];
 TH1D* h[npt][nalpha+1][ncent];
 TH1D *hsum[nalpha+1][ncent];
 TH1D * hmpt_vs_alpha[npt][ncent];
 
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
 
 double mpt_diff;
 double mpterr_diff;
 double mptsum_diff[nalpha+1][ncent];
 double mptsumpos_diff[nalpha+1][ncent];
 double mptsumneg_diff[nalpha+1][ncent];
 TH1D *hsum_diff[nalpha+1][ncent];
 TH1D * hmpt_vs_alpha_diff[npt][ncent];
 
 // int col[]={kAzure-9,kYellow-9,kOrange+1,kGreen+3,kPink-2,kRed+1,1};
 int col[]={kRed+1,kGreen+3,kOrange+1,kYellow-9,kAzure-9,1};
 TString sintegrate[nalpha+1][npt];
 TString s[npt];
 for(int ipt=0; ipt<npt;ipt++){
  s[ipt]="";
 }
 for(int ialpha=0;ialpha<nalpha+1;ialpha++){
 
   mptsum_ref[ialpha]=0;
   mptsumpos_ref[ialpha]=0;
   mptsumneg_ref[ialpha]=0;
  for(int icent=0;icent<ncent;icent++){
   mptsum[ialpha][icent]=0;
   mptsumpos[ialpha][icent]=0;
   mptsumneg[ialpha][icent]=0;
   mptsum_diff[ialpha][icent]=0;
   mptsumpos_diff[ialpha][icent]=0;
   mptsumneg_diff[ialpha][icent]=0;
  }
 }
   
 TLegend *leg,*leg2;
 leg= new TLegend(0.52,0.2,0.99,0.5);
 leg->SetBorderSize(0); 
 leg->SetFillStyle(0);
 leg2 = new TLegend(0.22,0.2,0.6,0.5);
 leg2->SetBorderSize(0); 
 leg2->SetFillStyle(0); 
   
 double frac[]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,3,3.5,4,4.4,5};
 
 for(int ipt=0;ipt<npt;ipt++){ 
  if(!doMC){  
   f[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_PbPb_mptonly_20140323_v3/full_ntuple_hiForest_Jet80or95_GR_R_53_LV6_12Mar2014_0000CET_Track8_Jet21_0_pt%d_%d_akVs3Calo.root",(int)ptmin[ipt],(int)ptmax[ipt])); 
   f_ref[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_pp_mptonly_20140323_v3/full_ntuple_pt80_pp2013_P01_prod22_v81_merged_forest_0_pt%d_%d_ak3Calo.root",(int)ptmin[ipt],(int)ptmax[ipt])); 
  }
  
  if(doMC){
   // f[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_MC_mptonly_20140322_v3/full_ntuple_HiForest_Pythia_Hydjet_Jet80_Track8_Jet21_STARTHI53_LV1_merged_forest_0_pt%d_%d_akVs3Calo_doGenJet%d.root",(int)ptmin[ipt],(int)ptmax[ipt],doGenJet)); 
   f[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_MC_mptonly_20140322_v3/full_ntuple_HiForest_Pythia_Hydjet_Jet80_Track8_Jet21_STARTHI53_LV1_merged_forest_0_pt%d_%d_akVs3Calo_doGenJet%d.root",(int)ptmin[ipt],(int)ptmax[ipt],doGenJet)); 
   // f[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_MC_mptonly_20140325_pthat120/full_ntuple_HydjetDrum_Pyquen_Dijet120_FOREST_Track8_Jet24_FixedPtHat_v0_pt%d_%d_akVs3Calo_doGenJet%d.root",(int)ptmin[ipt],(int)ptmax[ipt],doGenJet)); 
   f_ref[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_PYTHIA_mptonly_20140323_v3/full_ntuple_pt80_pp2013_P01_prod22_v81_merged_forest_0_pt%d_%d_ak3Calo_doGenJet%d.root",(int)ptmin[ipt],(int)ptmax[ipt],doGenJet)); 
  }

  cout<<Form("nt_%s%s_%s%s",smpt[domp].Data(),skind[gensigbkgreco].Data(),svariable[index_var].Data(),strave.Data())<<endl; 

  f[ipt]->cd();
  t_jet[ipt]=(TTree*)f[ipt]->Get("nt_jet");
  f[ipt]->cd();
  if(index_var==0) t[ipt]=(TTree*)f[ipt]->Get(Form("nt_%s%s%s",smpt[domp].Data(),skind[gensigbkgreco].Data(),strave.Data()));
  else t[ipt]=(TTree*)f[ipt]->Get(Form("nt_%s%s_%s%s",smpt[domp].Data(),skind[gensigbkgreco].Data(),svariable[index_var].Data(),strave.Data()));
  cout<<0<<endl;
  cout<<0.1<<endl;
  t[ipt]->AddFriend(t_jet[ipt]);
  cout<<0.2<<endl;

 cout<<Form("nt_%s%s_%s%s",smpt[domp].Data(),skind[gensigbkgreco].Data(),svariable[index_var].Data(),strave.Data())<<endl; 

  cout<<Form("nt_%s%s%s",smpt[domp].Data(),skind[gensigbkgreco].Data(),strave.Data())<<endl;
  if(index_var==0) t_ref[ipt]=(TTree*)f_ref[ipt]->Get(Form("nt_%s%s%s",smpt[domp].Data(),skind[gensigbkgreco].Data(),strave.Data()));
  else t_ref[ipt]=(TTree*)f_ref[ipt]->Get(Form("nt_%s%s_%s%s",smpt[domp].Data(),skind[gensigbkgreco].Data(),svariable[index_var].Data(),strave.Data()));
  t_jet_ref[ipt]=(TTree*)f_ref[ipt]->Get("nt_jet");
  t_ref[ipt]->AddFriend(t_jet_ref[ipt]);
  
  cout<<0.3<<endl;
  hmpt_vs_alpha_ref[ipt]= new TH1D(Form("hmpt_ref_%d",ipt),"",nalpha+1,frac);
   
  for(int ialpha=0;ialpha<nalpha+1;ialpha++){
  cout<<1<<endl;
   // if(ialpha<nalpha)s[ipt]=Form("-%s_%d",variable.Data(),ialpha);
   // if(ialpha<nalpha)sintegrate[ialpha][ipt]=s[ipt];
   // else  sintegrate[ialpha][ipt]=Form("-%s%s%s",smpt[domp].Data(),skind[gensigbkgreco].Data(),strave.Data());
    if(doIntegrate){
    if(ialpha<nalpha)s[ipt]=Form("-%s_%d",svariable[index_var].Data(),ialpha);
    if(ialpha<nalpha)sintegrate[ialpha][ipt]=s[ipt];
    else  sintegrate[ialpha][ipt]=Form("-%s%s%s",smpt[domp].Data(),skind[gensigbkgreco].Data(),strave.Data());
   }else{
    if(ialpha<nalpha && ialpha>0) sintegrate[ialpha][ipt]=Form("-%s_%d+%s_%d",svariable[index_var].Data(),ialpha,svariable[index_var].Data(),ialpha-1);
    else if(ialpha==0)sintegrate[ialpha][ipt]=Form("-%s_%d",svariable[index_var].Data(),ialpha);
    else  sintegrate[ialpha][ipt]=Form("-%s%s%s+%s_%d",smpt[domp].Data(),skind[gensigbkgreco].Data(),strave.Data(),svariable[index_var].Data(),nalpha-1);
   }
   h_ref[ipt][ialpha]=new TH1D(Form("h_ref_%d_%d",ialpha,ipt),"",100000,-100000,100000);
   t_ref[ipt]->Draw(Form("%s>>h_ref_%d_%d",sintegrate[ialpha][ipt].Data(),ialpha,ipt),Form("pt1>%d && pt2>%d && abs(eta1)<%.1f && abs(eta2)<%.1f && dphi>(9*TMath::Pi()/10)&& ((pt1-pt2)/(pt2+pt1))>%.2f",jetpt1,jetpt2,etadijet,etadijet,Aj),"");
   mpt_ref[ipt][ialpha]=h_ref[ipt][ialpha]->GetMean();
   mpterr_ref[ipt][ialpha]=h_ref[ipt][ialpha]->GetMeanError();
   
   if(ipt<npt-1){
    mptsum_ref[ialpha]+=mpt_ref[ipt][ialpha];
    if(mpt_ref[ipt][ialpha]>=0){
     mptsumpos_ref[ialpha]+=mpt_ref[ipt][ialpha];
     hmpt_vs_alpha_ref[ipt]->SetBinContent(ialpha+1,mptsumpos_ref[ialpha]);
     hmpt_vs_alpha_ref[ipt]->SetBinError(ialpha+1,mpterr_ref[ipt][ialpha]);
    }
    if(mpt_ref[ipt][ialpha]<0){  
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
  }
  
  for(int icent=0;icent<ncent;icent++){
  hmpt_vs_alpha[ipt][icent]= new TH1D(Form("hmpt_%d_%d",ipt,icent),"",nalpha+1,frac);
  hmpt_vs_alpha_diff[ipt][icent]= new TH1D(Form("hmpt_diff_%d_%d",ipt,icent),"",nalpha+1,frac);

  for(int ialpha=0;ialpha<nalpha+1;ialpha++ ){
   h[ipt][ialpha][icent]=new TH1D(Form("h_%d_%d_%d",ialpha,ipt,icent),"",100000,-100000,100000);
   t[ipt]->Draw(Form("%s>>h_%d_%d_%d",sintegrate[ialpha][ipt].Data(),ialpha,ipt,icent),Form("pt1>%d && pt2>%d && abs(eta1)<%.1f && abs(eta2)<%.1f && dphi>(9*TMath::Pi()/10) && cent>=%d && cent<%d && ((pt1-pt2)/(pt2+pt1))>%.2f",jetpt1,jetpt2,etadijet,etadijet, cent_min[icent],cent_max[icent],Aj),"");
   mpt=h[ipt][ialpha][icent]->GetMean(); 
   mpterr=h[ipt][ialpha][icent]->GetMeanError();
   
   mpt_diff=mpt-mpt_ref[ipt][ialpha];
   mpterr_diff=sqrt(pow(mpterr,2)+pow(mpterr_ref[ipt][ialpha],2));
   
   if(ipt<npt-1){
    mptsum[ialpha][icent]+=mpt;
    if(mpt>=0){
     mptsumpos[ialpha][icent]+=mpt;
     hmpt_vs_alpha[ipt][icent]->SetBinContent(ialpha+1,mptsumpos[ialpha][icent]);
     hmpt_vs_alpha[ipt][icent]->SetBinError(ialpha+1,mpterr);
    }
    if(mpt<0){
     mptsumneg[ialpha][icent]+=mpt;
     hmpt_vs_alpha[ipt][icent]->SetBinContent(ialpha+1,mptsumneg[ialpha][icent]);
     hmpt_vs_alpha[ipt][icent]->SetBinError(ialpha+1,mpterr);
    } 
    hmpt_vs_alpha[ipt][icent]->SetLineColor(1);
    hmpt_vs_alpha[ipt][icent]->SetMarkerColor(1);
    hmpt_vs_alpha[ipt][icent]->SetMarkerSize(0);
    hmpt_vs_alpha[ipt][icent]->SetFillColor(col[ipt]); 
   }else{
     hmpt_vs_alpha[ipt][icent]->SetBinContent(ialpha+1,mpt);
     hmpt_vs_alpha[ipt][icent]->SetBinError(ialpha+1,mpterr);
   }   
   
   if(ipt<npt-1){
    mptsum_diff[ialpha][icent]+=mpt_diff;
    if(mpt_diff>=0){
     mptsumpos_diff[ialpha][icent]+=mpt_diff;
     hmpt_vs_alpha_diff[ipt][icent]->SetBinContent(ialpha+1,mptsumpos_diff[ialpha][icent]);
     hmpt_vs_alpha_diff[ipt][icent]->SetBinError(ialpha+1,mpterr_diff);
    }
    if(mpt_diff<0){  
     mptsumneg_diff[ialpha][icent]+=mpt_diff;
     hmpt_vs_alpha_diff[ipt][icent]->SetBinContent(ialpha+1,mptsumneg_diff[ialpha][icent]);
     hmpt_vs_alpha_diff[ipt][icent]->SetBinError(ialpha+1,mpterr_diff);
    }
    hmpt_vs_alpha_diff[ipt][icent]->SetLineColor(1);
    hmpt_vs_alpha_diff[ipt][icent]->SetMarkerColor(1);
    hmpt_vs_alpha_diff[ipt][icent]->SetMarkerSize(0);
    hmpt_vs_alpha_diff[ipt][icent]->SetFillColor(col[ipt]);
   }else{ 
     hmpt_vs_alpha_diff[ipt][icent]->SetBinContent(ialpha+1,mpt_diff);
     hmpt_vs_alpha_diff[ipt][icent]->SetBinError(ialpha+1,mpterr_diff);
   }
  }
 }
 }
 
 // cout<<"all events"<<endl;  

  TCanvas *c1 = new TCanvas("c1","",(ncent+1)*300,700);
  makeMultiPanelCanvas(c1,ncent+1,2,0.0,0.0,0.2,0.2,0.02);
 TH1D * empty=new TH1D("empty",Form(";%s;<#slash{p}_{T}^{#parallel}> (GeV)",axistitle[index_var].Data()),nalpha+1,frac);
 empty->Fill(0.5,1000); 
 if(doIntegrate){
  if(index_var==0){
   empty->SetMaximum(30); 
   empty->SetMinimum(-60); 
  
  }else{
   empty->SetMaximum(30); 
   empty->SetMinimum(-40); 
  }
 }else{
  empty->SetMaximum(5); 
  empty->SetMinimum(-10); 
 
 }
    empty->GetXaxis()->SetTitleSize(28);
   empty->GetXaxis()->SetTitleFont(43); 
   empty->GetXaxis()->SetTitleOffset(2.2);
   empty->GetXaxis()->SetLabelSize(22);
   empty->GetXaxis()->SetLabelFont(43);
   empty->GetYaxis()->SetTitleSize(28);
   empty->GetYaxis()->SetTitleFont(43); 
   empty->GetYaxis()->SetTitleOffset(2.2);
   empty->GetYaxis()->SetLabelSize(22);
   empty->GetYaxis()->SetLabelFont(43);
   
   
 c1->cd(ncent+2); 
 leg2->Draw("same");
 leg->Draw("same");

 drawText(Form("|#eta_{1}|,|#eta_{2}|<%.2f, A_{J} > %.1f",etadijet,Aj),0.22,0.54);
 drawText(Form("p_{T,1}>%d, p_{T,2}>%d",jetpt1,jetpt2),0.22,0.9);
 drawText(Form("|#eta|<2, #Delta#phi_{1,2}>9#pi/10",jetpt1,jetpt2),0.22,0.81);
 if(!doGenJet)drawText(Form("akVs3Calo, ak3Calo"),0.22,0.72);
 else drawText(Form("gen jet"),0.22,0.72);
 drawText(Form("%s",spart[gensigbkgreco].Data()),0.22,0.63);
 
 c1->cd(1);
 empty->Draw();

 for(int ipt=npt-2;ipt>=0;ipt--){
  hmpt_vs_alpha_ref[ipt]->Draw("same");
  hmpt_vs_alpha_ref[ipt]->Draw("same hist");
 }
 hmpt_vs_alpha_ref[npt-1]->Draw("same");
 c1->cd(1)->RedrawAxis();
 if(!doMC)drawText("pp",0.22,0.9);
 else drawText("PYTHIA",0.22,0.9);

 for(int icent=0;icent<ncent;icent++){
 c1->cd(ncent+1-icent);
 empty->Draw();
  for(int ipt=npt-2;ipt>=0;ipt--){ 
   if(icent==0){
    if(ipt>2) leg2->AddEntry(hmpt_vs_alpha[ipt][icent],Form("%.1f - %.1f GeV",ptmin[ipt],ptmax[ipt]),"f l");
    else leg->AddEntry(hmpt_vs_alpha[ipt][icent],Form("%.1f - %.1f GeV",ptmin[ipt],ptmax[ipt]),"f l");
   }
   hmpt_vs_alpha[ipt][icent]->Draw("same");
   hmpt_vs_alpha[ipt][icent]->Draw("same hist");
  }
  if(icent==0)leg2->AddEntry(hmpt_vs_alpha[npt-1][icent],">0.5 GeV","p");
  hmpt_vs_alpha[npt-1][icent]->Draw("same");
  drawText(Form("%d-%d %%",(int)(0.5*cent_min[icent]),(int)(0.5*cent_max[icent])),0.22,0.85);
   c1->cd(ncent+1-icent)->RedrawAxis();
   if(!doMC)drawText("PbPb",0.22,0.9);
   else drawText("PYtHIA+HYDJET",0.22,0.9);

 }
 for(int icent=0;icent<ncent;icent++){
 c1->cd((2*(ncent+1))-icent);
  empty->Draw();
  for(int ipt=npt-2;ipt>=0;ipt--){ 
   hmpt_vs_alpha_diff[ipt][icent]->Draw("same");
   hmpt_vs_alpha_diff[ipt][icent]->Draw("same hist");
  }
  hmpt_vs_alpha_diff[npt-1][icent]->Draw("same");
  drawText(Form("%d-%d %%",(int)(0.5*cent_min[icent]),(int)(0.5*cent_max[icent])),0.22,0.85);
   c1->cd((2*(ncent+1))-icent)->RedrawAxis();
   if(!doMC)drawText("PbPb-pp",0.22,0.9);
   else drawText("(PYT.+HYD.)-PYT.",0.22,0.9);
 }
 
 c1->SaveAs(Form("multipanel_eta%d_doMC%d_%s_doIntegrate%d_doGenJet%d_kind%d_Aj%d.png",(int)(etadijet*10),doMC,svariable[index_var].Data(),doIntegrate,doGenJet,gensigbkgreco,(int)(Aj*10)));
 c1->SaveAs(Form("multipanel_eta%d_doMC%d_%s_doIntegrate%d_doGenJet%d_kind%d_Aj%d.pdf",(int)(etadijet*10),doMC,svariable[index_var].Data(),doIntegrate,doGenJet,gensigbkgreco,(int)(Aj*10)));
 
}