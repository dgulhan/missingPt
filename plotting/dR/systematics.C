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
#include "TGraphErrors.h"
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

void drawPatch(float x1, float y1, float x2, float y2){
  TLegend *t1=new TLegend(x1,y1,x2,y2);
  t1->SetFillColor(kWhite);
  t1->SetBorderSize(0);
  t1->SetFillStyle(1001);
  t1->Draw("");
}


void cumulative(int gensigbkgreco=3, int doMC=0, int index_var=2, bool doIntegrate=true,  int iAj=0, int doAve=1){
 TH1D::SetDefaultSumw2(); 
 bool doCompareGenReco= false;
 int domp=0;
 int ncent=2;  
 bool doGenJet=false;
 double etadijet=0.5;
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
 int npt=6;
 // int  cent_min[]={0,20,40,60,100};
 // int cent_max[]={20,40,60,100,200};
 int  cent_min[]={0,60};
 int cent_max[]={60,200};
 double ptmin[]={  8,4,2,1,0.5,0.5};
 double ptmax[]={300,8,4,2,  1,300};
 double systematics[] = { 3, 3,2.5};
 double Ajmin[]={0,0.22,0};
 double Ajmax[]={0.22,1,1};
 TFile *f[npt];
 TTree *t[npt];
 TTree *t_jet[npt];  
 int nalpha;
 if(index_var==0 ||index_var==2) nalpha=16;
 else nalpha=18;
 double mpt;
 double mpterr;  
 double mptsum[nalpha+1][ncent];
 double mptsumpos[nalpha+1][ncent];
 double mptsumneg[nalpha+1][ncent];
 TH1D* h[npt][nalpha+1][ncent];
 TH1D* h_stat[npt][nalpha+1][ncent];
 TH1D *hsum[nalpha+1][ncent];
 TH1D * hmpt_vs_alpha[npt][ncent];
 
 TFile *f_ref[npt];
 TTree *t_ref[npt];
 TTree *t_jet_ref[npt]; 
 double mpt_ref[npt][nalpha+1][ncent];
 double mpterr_ref[npt][nalpha+1][ncent];
 double mptsum_ref[nalpha+1][ncent];
 double mptsumpos_ref[nalpha+1][ncent];
 double mptsumneg_ref[nalpha+1][ncent];
 TH1D* h_ref[npt][nalpha+1][ncent];
 TH1D* h_ref_stat[npt][nalpha+1][ncent];
 TH1D *hsum_ref[nalpha+1][ncent];
 TH1D * hmpt_vs_alpha_ref[npt][ncent];
 
 double mpt_diff;
 double mpterr_diff;
 double mptsum_diff[nalpha+1][ncent];
 double mptsumpos_diff[nalpha+1][ncent];
 double mptsumneg_diff[nalpha+1][ncent];
 TH1D *hsum_diff[nalpha+1][ncent];
 TH1D * hmpt_vs_alpha_diff[npt][ncent];
 
 // int col[]={kAzure-9,kYellow-9,kOrange+1,kGreen+3,kPink-2,kRed+1,1};
 int col[]={kRed+1,kGreen+3,kOrange+1,kYellow-9,kBlue-9,1};
 TString sintegrate[nalpha+1][npt];
 TString sintegrate_stat[nalpha+1][npt];
 TString sintegrate_ref[nalpha+1][npt];
 TString sintegrate_stat_ref[nalpha+1][npt];
 TString s[npt];

 for(int ipt=0; ipt<npt;ipt++){
  s[ipt]="";
 }
 for(int ialpha=0;ialpha<nalpha+1;ialpha++){
 
  for(int icent=0;icent<ncent;icent++){
   mptsum_ref[ialpha][icent]=0;
   mptsumpos_ref[ialpha][icent]=0;
   mptsumneg_ref[ialpha][icent]=0;
   mptsum[ialpha][icent]=0;
   mptsumpos[ialpha][icent]=0;
   mptsumneg[ialpha][icent]=0;
   mptsum_diff[ialpha][icent]=0;
   mptsumpos_diff[ialpha][icent]=0;
   mptsumneg_diff[ialpha][icent]=0;
  }
 }
   
 TLegend *leg,*leg2, *leg3,*leg4;
 leg= new TLegend(0.45,0.1,0.99,0.4);
 leg->SetBorderSize(0); 
 leg->SetFillStyle(0);
 leg->SetTextFont(43);
 leg->SetTextSize(18);
 leg3 = new TLegend(0.35,0.1,0.95,0.4);
 leg3->SetBorderSize(0); 
 leg3->SetFillStyle(0); 
 leg3->SetTextFont(43);
 leg3->SetTextSize(18);
 leg4 = new TLegend(0.35,0.3,0.95,0.25);
 leg4->SetBorderSize(0); 
 leg4->SetFillStyle(0); 
 leg4->SetTextFont(43);
 leg4->SetTextSize(18);
 leg2 = new TLegend(0.05,0.1,0.6,0.4); 
 leg2->SetBorderSize(0); 
 leg2->SetFillStyle(0); 
 leg2->SetTextFont(43);
 leg2->SetTextSize(18);
 double frac[]={0,0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2,2.2,2.4,3,4,5};
 
 for(int ipt=0;ipt<npt;ipt++){ 
  if(!doMC){  
   f[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_PbPb_mptonly_20140325_v5/full_ntuple_hiForest_Jet80or95_GR_R_53_LV6_03Mar2014_1600CET_CMSSW_5_3_16_merged_pt%d_%d_akVs3Calo.root",(int)ptmin[ipt],(int)ptmax[ipt])); 
   f_ref[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_pp_mptonly_20140325_v4/full_ntuple_PP2013_HiForest_PromptReco_JsonPP_Jet80_PPReco_forestv82_pt%d_%d_ak3Calo.root",(int)ptmin[ipt],(int)ptmax[ipt])); 
  }
  
  if(doMC){
   f[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_MC_20140419//full_ntuple_HydjetDrum_Pyquen_Dijet_FOREST_Track8_Jet24_FixedPtHatJES_v0_pt%d_%d_akVs3Calo_doGenJet0.root",(int)ptmin[ipt],(int)ptmax[ipt],doGenJet)); 
   f_ref[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_MC_20140419//full_ntuple_HydjetDrum_Pyquen_Dijet_FOREST_Track8_Jet24_FixedPtHatJES_v0_pt%d_%d_akVs3Calo_doGenJet1.root",(int)ptmin[ipt],(int)ptmax[ipt],doGenJet)); 
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

 cout<<Form("nt_%s%s_%s%s",smpt[domp].Data(),skind[0].Data(),svariable[index_var].Data(),strave.Data())<<endl; 

  cout<<Form("nt_%s%s%s",smpt[domp].Data(),skind[0].Data(),strave.Data())<<endl;
  if(index_var==0) t_ref[ipt]=(TTree*)f_ref[ipt]->Get(Form("nt_%s%s%s",smpt[domp].Data(),skind[0].Data(),strave.Data()));
  else t_ref[ipt]=(TTree*)f_ref[ipt]->Get(Form("nt_%s%s_%s%s",smpt[domp].Data(),skind[0].Data(),svariable[index_var].Data(),strave.Data()));
  t_jet_ref[ipt]=(TTree*)f_ref[ipt]->Get("nt_jet");
  t_ref[ipt]->AddFriend(t_jet_ref[ipt]);
  
  cout<<0.3<<endl;
 for(int icent=0;icent<ncent;icent++){

  hmpt_vs_alpha_ref[ipt][icent]= new TH1D(Form("hmpt_ref_%d_%d",ipt,icent),"",nalpha/2+1,frac);
   
  for(int ialpha=0;ialpha<nalpha/2+1;ialpha++){
  cout<<1<<endl;
   if(icent==0){
   if(ipt==npt-1){
    // if(ialpha<nalpha/2)s[ipt]=Form("-%s_%d",svariable[index_var].Data(),2*ialpha+1);
    // if(ialpha<nalpha/2)sintegrate[ialpha][ipt]=s[ipt];
    // else  sintegrate[ialpha][ipt]=Form("-%s%s%s",smpt[domp].Data(),skind[gensigbkgreco].Data(),strave.Data());
	
    if(ialpha<nalpha/2 && ialpha>0) sintegrate[ialpha][ipt]=Form("-%s_%d+%s_%d",svariable[index_var].Data(),2*ialpha+1,svariable[index_var].Data(),2*(ialpha-1)+1);
    else if(ialpha==0)sintegrate[ialpha][ipt]=Form("-%s_%d",svariable[index_var].Data(),ialpha+1);
    else  sintegrate[ialpha][ipt]=Form("-%s%s%s+%s_%d",smpt[domp].Data(),skind[gensigbkgreco].Data(),strave.Data(),svariable[index_var].Data(),nalpha-1);
    if(ialpha<nalpha/2 && ialpha>0) sintegrate_stat[ialpha][ipt]=Form("-%s_%d+%s_%d",svariable[index_var].Data(),2*ialpha+1,svariable[index_var].Data(),2*(ialpha-1)+1);
    else if(ialpha==0)sintegrate_stat[ialpha][ipt]=Form("-%s_%d",svariable[index_var].Data(),ialpha+1);
    else  sintegrate_stat[ialpha][ipt]=Form("-%s%s%s+%s_%d",smpt[domp].Data(),skind[gensigbkgreco].Data(),strave.Data(),svariable[index_var].Data(),nalpha-1);
   }else{
    if(ialpha<nalpha/2 && ialpha>0) sintegrate[ialpha][ipt]=Form("-%s_%d+%s_%d",svariable[index_var].Data(),2*ialpha+1,svariable[index_var].Data(),2*(ialpha-1)+1);
    else if(ialpha==0)sintegrate[ialpha][ipt]=Form("-%s_%d",svariable[index_var].Data(),ialpha+1);
    else  sintegrate[ialpha][ipt]=Form("-%s%s%s+%s_%d",smpt[domp].Data(),skind[gensigbkgreco].Data(),strave.Data(),svariable[index_var].Data(),nalpha-1);
   }
   if(ipt==npt-1){
    if(ialpha<nalpha/2 && ialpha>0) sintegrate_ref[ialpha][ipt]=Form("-%s_%d+%s_%d",svariable[index_var].Data(),2*ialpha+1,svariable[index_var].Data(),2*(ialpha-1)+1);
    else if(ialpha==0)sintegrate_ref[ialpha][ipt]=Form("-%s_%d",svariable[index_var].Data(),ialpha+1);
    else  sintegrate_ref[ialpha][ipt]=Form("-%s%s%s+%s_%d",smpt[domp].Data(),skind[0].Data(),strave.Data(),svariable[index_var].Data(),nalpha-1);
    if(ialpha<nalpha/2 && ialpha>0) sintegrate_stat_ref[ialpha][ipt]=Form("-%s_%d+%s_%d",svariable[index_var].Data(),2*ialpha+1,svariable[index_var].Data(),2*(ialpha-1)+1);
    else if(ialpha==0)sintegrate_stat_ref[ialpha][ipt]=Form("-%s_%d",svariable[index_var].Data(),ialpha+1);
    else  sintegrate_stat_ref[ialpha][ipt]=Form("-%s%s%s+%s_%d",smpt[domp].Data(),skind[0].Data(),strave.Data(),svariable[index_var].Data(),nalpha-1);
   }else{
    if(ialpha<nalpha/2 && ialpha>0) sintegrate_ref[ialpha][ipt]=Form("-%s_%d+%s_%d",svariable[index_var].Data(),2*ialpha+1,svariable[index_var].Data(),2*(ialpha-1)+1);
    else if(ialpha==0)sintegrate_ref[ialpha][ipt]=Form("-%s_%d",svariable[index_var].Data(),ialpha+1);
    else  sintegrate_ref[ialpha][ipt]=Form("-%s%s%s+%s_%d",smpt[domp].Data(),skind[0].Data(),strave.Data(),svariable[index_var].Data(),nalpha-1);
   }
   }
   h_ref[ipt][ialpha][icent]=new TH1D(Form("h_ref_%d_%d_%d",ialpha,ipt,icent),"",100000,-100000,100000);
   t_ref[ipt]->Draw(Form("%s>>h_ref_%d_%d_%d",sintegrate_ref[ialpha][ipt].Data(),ialpha,ipt,icent),Form("weight*cent_weight*(pthat>80 && pt1>%d && pt2>%d && abs(eta1)<%.1f && abs(eta2)<%.1f && dphi>(5*TMath::Pi()/6)&& ((pt1-pt2)/(pt2+pt1))>=%.2f && ((pt1-pt2)/(pt2+pt1))<%.2f && abs(vz)<15 && cent>=%d && cent<%d)",jetpt1,jetpt2,etadijet,etadijet,Ajmin[iAj],Ajmax[iAj],cent_min[icent],cent_max[icent]),"");
   
   mpt_ref[ipt][ialpha][icent]=h_ref[ipt][ialpha][icent]->GetMean();
   mpterr_ref[ipt][ialpha][icent]=h_ref[ipt][ialpha][icent]->GetMeanError();
   if(ipt==npt-1){
    h_ref_stat[ipt][ialpha][icent]=new TH1D(Form("h_ref_stat_%d_%d_%d",ialpha,ipt,icent),"",100000,-100000,100000);
    t_ref[ipt]->Draw(Form("%s>>h_ref_stat_%d_%d_%d",sintegrate_stat_ref[ialpha][ipt].Data(),ialpha,ipt,icent),Form("weight*cent_weight*(pthat>80 && pt1>%d && pt2>%d && abs(eta1)<%.1f && abs(eta2)<%.1f && dphi>(5*TMath::Pi()/6)&& ((pt1-pt2)/(pt2+pt1))>=%.2f && ((pt1-pt2)/(pt2+pt1))<%.2f && abs(vz)<15 && cent>=%d && cent<%d)",jetpt1,jetpt2,etadijet,etadijet,Ajmin[iAj],Ajmax[iAj],cent_min[icent],cent_max[icent]),"");
    mpterr_ref[ipt][ialpha][icent]=h_ref_stat[ipt][ialpha][icent]->GetMeanError();
   }
   
   if(ipt<npt-1){
    mptsum_ref[ialpha][icent]+=mpt_ref[ipt][ialpha][icent];
    if(mpt_ref[ipt][ialpha][icent]>=0){
     mptsumpos_ref[ialpha][icent]+=mpt_ref[ipt][ialpha][icent];
     hmpt_vs_alpha_ref[ipt][icent]->SetBinContent(ialpha+1,mptsumpos_ref[ialpha][icent]);
     hmpt_vs_alpha_ref[ipt][icent]->SetBinError(ialpha+1,mpterr_ref[ipt][ialpha][icent]);
    }
    if(mpt_ref[ipt][ialpha][icent]<0){  
     mptsumneg_ref[ialpha][icent]+=mpt_ref[ipt][ialpha][icent];
     hmpt_vs_alpha_ref[ipt][icent]->SetBinContent(ialpha+1,mptsumneg_ref[ialpha][icent]);
     hmpt_vs_alpha_ref[ipt][icent]->SetBinError(ialpha+1,mpterr_ref[ipt][ialpha][icent]);
    }
    hmpt_vs_alpha_ref[ipt][icent]->SetLineColor(1);
    hmpt_vs_alpha_ref[ipt][icent]->SetMarkerColor(1);
    hmpt_vs_alpha_ref[ipt][icent]->SetMarkerSize(0);
    hmpt_vs_alpha_ref[ipt][icent]->SetFillColor(col[ipt]);
	
   }else{ 
     hmpt_vs_alpha_ref[ipt][icent]->SetBinContent(ialpha+1,mpt_ref[ipt][ialpha][icent]);
     hmpt_vs_alpha_ref[ipt][icent]->SetBinError(ialpha+1,mpterr_ref[ipt][ialpha][icent]);
     hmpt_vs_alpha_ref[ipt][icent]->SetMarkerStyle(25);
   }
  }
  }
  for(int icent=0;icent<ncent;icent++){
  hmpt_vs_alpha[ipt][icent]= new TH1D(Form("hmpt_%d_%d",ipt,icent),"",nalpha/2+1,frac);
  hmpt_vs_alpha_diff[ipt][icent]= new TH1D(Form("hmpt_diff_%d_%d",ipt,icent),"",nalpha/2+1,frac);

  for(int ialpha=0;ialpha<nalpha/2+1;ialpha++ ){
   h[ipt][ialpha][icent]=new TH1D(Form("h_%d_%d_%d",ialpha,ipt,icent),"",100000,-100000,100000);
   h_stat[ipt][ialpha][icent]=new TH1D(Form("h_stat_%d_%d_%d",ialpha,ipt,icent),"",100000,-100000,100000);
   if(!doMC)t[ipt]->Draw(Form("%s>>h_%d_%d_%d",sintegrate[ialpha][ipt].Data(),ialpha,ipt,icent),Form("pt1>%d && pt2>%d && abs(eta1)<%.1f && abs(eta2)<%.1f && dphi>(5*TMath::Pi()/6) && cent>=%d && cent<%d && ((pt1-pt2)/(pt2+pt1))>=%.2f && ((pt1-pt2)/(pt2+pt1))<%.2f && abs(vz)<15",jetpt1,jetpt2,etadijet,etadijet, cent_min[icent],cent_max[icent],Ajmin[iAj],Ajmax[iAj]),"");
   else t[ipt]->Draw(Form("%s>>h_%d_%d_%d",sintegrate[ialpha][ipt].Data(),ialpha,ipt,icent),Form("weight*cent_weight*(pthat>80 && pt1>%d && pt2>%d && abs(eta1)<%.1f && abs(eta2)<%.1f && dphi>(5*TMath::Pi()/6) && cent>=%d && cent<%d && ((pt1-pt2)/(pt2+pt1))>=%.2f && ((pt1-pt2)/(pt2+pt1))<%.2f && abs(vz)<15)",jetpt1,jetpt2,etadijet,etadijet, cent_min[icent],cent_max[icent],Ajmin[iAj],Ajmax[iAj]),"");
    mpt=h[ipt][ialpha][icent]->GetMean(); 
    mpterr=h[ipt][ialpha][icent]->GetMeanError();
   
   mpt=h[ipt][ialpha][icent]->GetMean(); 
   mpterr=h[ipt][ialpha][icent]->GetMeanError();
   
   if(ipt==npt-1){
    if(!doMC)t[ipt]->Draw(Form("%s>>h_stat_%d_%d_%d",sintegrate_stat[ialpha][ipt].Data(),ialpha,ipt,icent),Form("pt1>%d && pt2>%d && abs(eta1)<%.1f && abs(eta2)<%.1f && dphi>(5*TMath::Pi()/6) && cent>=%d && cent<%d && ((pt1-pt2)/(pt2+pt1))>=%.2f && ((pt1-pt2)/(pt2+pt1))<%.2f && abs(vz)<15",jetpt1,jetpt2,etadijet,etadijet, cent_min[icent],cent_max[icent],Ajmin[iAj],Ajmax[iAj]),"");
    else t[ipt]->Draw(Form("%s>>h_stat_%d_%d_%d",sintegrate_stat[ialpha][ipt].Data(),ialpha,ipt,icent),Form("weight*cent_weight*(pthat>80 &&pt1>%d && pt2>%d && abs(eta1)<%.1f && abs(eta2)<%.1f && dphi>(5*TMath::Pi()/6) && cent>=%d && cent<%d && ((pt1-pt2)/(pt2+pt1))>=%.2f && ((pt1-pt2)/(pt2+pt1))<%.2f && abs(vz)<15)",jetpt1,jetpt2,etadijet,etadijet, cent_min[icent],cent_max[icent],Ajmin[iAj],Ajmax[iAj]),""); 
    mpterr=h_stat[ipt][ialpha][icent]->GetMeanError();
   }
   mpt_diff=mpt-mpt_ref[ipt][ialpha][icent];
   mpterr_diff=sqrt(pow(mpterr,2)+pow(mpterr_ref[ipt][ialpha][icent],2));
   
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
   hmpt_vs_alpha_diff[ipt][icent]->SetMarkerStyle(24);

  }
 }
 }
 
 // cout<<"all events"<<endl;  
  TLine* zeroLine_p = new TLine(0., 0., 2., 0.);
  zeroLine_p->SetLineColor(1);
  zeroLine_p->SetLineStyle(2);
  zeroLine_p->Draw();
  TGraphErrors * g_hmpt_diff[ncent];
  double graphyerr[ncent][hmpt_vs_alpha_diff[npt-1][0]->GetNbinsX()];
  double graphy[ncent][hmpt_vs_alpha_diff[npt-1][0]->GetNbinsX()];
  double graphx[hmpt_vs_alpha_diff[npt-1][0]->GetNbinsX()];
  double graphxerr[hmpt_vs_alpha_diff[npt-1][0]->GetNbinsX()];
  
  TCanvas *c1 = new TCanvas("c1","",(ncent)*300,1050);
  makeMultiPanelCanvas(c1,ncent,3,0.0,0.0,0.2,0.2,0.02);
 TH1D * empty=new TH1D("empty",Form(";%s;<#slash{p}_{T}^{#parallel}> (GeV)",axistitle[index_var].Data()),nalpha/2+1,frac);
 TH1D * empty2=new TH1D("empty2",Form(";%s;<#slash{p}_{T}^{#parallel}> (GeV)",axistitle[index_var].Data()),nalpha/2+1,frac);
 empty->Fill(0.5,1000); 
 empty2->Fill(0.5,1000); 
 if(doIntegrate){
  if(index_var==0){
   empty->SetMaximum(20); 
   empty->SetMinimum(-55); 
  
  }else{
   empty->SetMaximum(20); 
   empty->SetMinimum(-55); 
  }
 }else{
  empty->SetMaximum(20); 
  empty->SetMinimum(-55); 
  empty2->SetMaximum(10); 
  empty2->SetMinimum(-10); 
 
 }
    empty->GetXaxis()->CenterTitle();
    empty->GetYaxis()->CenterTitle();
    empty2->GetXaxis()->CenterTitle();
    empty2->GetYaxis()->CenterTitle();
    empty->GetXaxis()->SetNdivisions(505);
    empty->GetXaxis()->SetNdivisions(505);
    empty2->GetXaxis()->SetNdivisions(505);
    empty2->GetYaxis()->SetNdivisions(505);
    empty->GetYaxis()->SetNdivisions(505);
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
    empty2->GetXaxis()->SetTitleSize(28);
   empty2->GetXaxis()->SetTitleFont(43); 
   empty2->GetXaxis()->SetTitleOffset(2.2);
   empty2->GetXaxis()->SetLabelSize(22);
   empty2->GetXaxis()->SetLabelFont(43);
   empty2->GetYaxis()->SetTitleSize(28);
   empty2->GetYaxis()->SetTitleFont(43); 
   empty2->GetYaxis()->SetTitleOffset(1.8);
   empty2->GetYaxis()->SetLabelSize(22);
   empty2->GetYaxis()->SetLabelFont(43);
   
   
 // c1->cd(ncent+2); 
 // leg2->Draw("same");
 // leg->Draw("same");

 // if(Ajmin[iAj]==0 && Ajmax[iAj]==1) drawText(Form("|#eta_{1}|,|#eta_{2}|<%.2f",etadijet),0.05,0.42);
 // else if(Ajmin[iAj]==0) drawText(Form("|#eta_{1}|,|#eta_{2}|<%.2f, A_{J} < %.2f",etadijet,Ajmax[iAj]),0.05,0.42);
 // else if(Ajmax[iAj]==1)drawText(Form("|#eta_{1}|,|#eta_{2}|<%.2f,  A_{J} > %.2f",etadijet,Ajmin[iAj]),0.05,0.42);
 // else drawText(Form("|#eta_{1}|,|#eta_{2}|<%.2f, %.2f < A_{J} < %.0f",etadijet,Ajmin[iAj],Ajmax[iAj]),0.05,0.42);
 // drawText(Form("p_{T,1}>%d, p_{T,2}>%d",jetpt1,jetpt2),0.05,0.9);
 // drawText(Form("|#eta|<2, #Delta#phi_{1,2}>5#pi/6"),0.05,0.78);
 // if(!doMC){
   // drawText(Form("PbPb #sqrt{s_{NN}}=2.76 TeV 150/#mub"),0.05,0.66);
   // drawText(Form("pp #sqrt{s_{NN}}=2.76 TeV 5.3/pb"),0.05,0.54);
 // }
 // if(doMC){
  // if(!doGenJet)drawText(Form("akVs3Calo, ak3Calo"),0.05,0.66);
  // else drawText(Form("gen jet"),0.05,0.66);
  // drawText(Form("%s",spart[gensigbkgreco].Data()),0.05,0.54);
 // }
 
 
 for(int icent=0;icent<ncent;icent++){ 
 c1->cd((3*(ncent))-icent);
 
 if(doIntegrate) empty->Draw();
   else empty2->Draw();
  for(int ipt=npt-2;ipt>=0;ipt--){ 
   hmpt_vs_alpha_diff[ipt][icent]->Draw("same");
   hmpt_vs_alpha_diff[ipt][icent]->Draw("same hist");
  }
  for(int ialpha = 0 ;ialpha<hmpt_vs_alpha_diff[npt-1][icent]->GetNbinsX();ialpha++){
   double  y = hmpt_vs_alpha_diff[npt-1][icent]->GetBinContent(ialpha+1);
   double yerr = systematics[iAj];
   double x=hmpt_vs_alpha_diff[npt-1][icent]->GetBinCenter(ialpha+1);
  }
  

  hmpt_vs_alpha_diff[npt-1][icent]->Draw("same");
  zeroLine_p->Draw("same");
  if(icent==1)drawText(Form("%d-%d %%",(int)(0.5*cent_min[icent]),(int)(0.5*cent_max[icent])),0.3,0.85);
  else drawText(Form("%d-%d %%",(int)(0.5*cent_min[icent]),(int)(0.5*cent_max[icent])),0.22,0.85);
   c1->cd((3*(ncent))-icent)->RedrawAxis();
   if(!doMC)drawText("PbPb-pp",0.22,0.9);
   else drawText("reco jet-gen jet w. gen part.",0.22,0.9);
   // else drawText("reco-gen",0.22,0.9);
   drawText("PYTHIA+HYDJET",0.25,0.45);
 if(icent==0){ 
  if(Ajmin[iAj]==0 && Ajmax[iAj]==1) drawText(Form("|#eta_{1}|,|#eta_{2}|<%.2f",etadijet),0.25,0.25);
  else if(Ajmin[iAj]==0) drawText(Form("|#eta_{1}|,|#eta_{2}|<%.2f, A_{J} < %.2f",etadijet,Ajmax[iAj]),0.25,0.25);
  else if(Ajmax[iAj]==1)drawText(Form("|#eta_{1}|,|#eta_{2}|<%.2f,  A_{J} > %.2f",etadijet,Ajmin[iAj]),0.25,0.25);
  else drawText(Form("|#eta_{1}|,|#eta_{2}|<%.2f, %.2f < A_{J} < %.0f",etadijet,Ajmin[iAj],Ajmax[iAj]),0.25,0.25);
 }else{
 
  drawText(Form("p_{T,1}>%d, p_{T,2}>%d",jetpt1,jetpt2),0.25,0.25);
  drawText(Form("|#eta|<2, #Delta#phi_{1,2}>5#pi/6"),0.25,0.35);

 }
 }

 for(int icent=0;icent<ncent;icent++){
 c1->cd(2*ncent-icent);
 empty->Draw();
  for(int ipt=npt-2;ipt>=0;ipt--){ 
   if(icent==0){
    if(ipt>2) leg2->AddEntry(hmpt_vs_alpha[ipt][icent],Form("%.1f - %.1f",ptmin[ipt],ptmax[ipt]),"f");
    else leg->AddEntry(hmpt_vs_alpha[ipt][icent],Form("%.1f - %.1f",ptmin[ipt],ptmax[ipt]),"f");
   }
   hmpt_vs_alpha[ipt][icent]->Draw("same");
   hmpt_vs_alpha[ipt][icent]->Draw("same hist");
  }
  
  if(icent==0)leg2->AddEntry(hmpt_vs_alpha[npt-1][icent],">0.5","p");
  hmpt_vs_alpha[npt-1][icent]->Draw("same");
  drawText(Form("%d-%d %%",(int)(0.5*cent_min[icent]),(int)(0.5*cent_max[icent])),0.22,0.85);
  drawText(Form("%d-%d %%",(int)(0.5*cent_min[icent]),(int)(0.5*cent_max[icent])),0.22,0.85);
   zeroLine_p->Draw("same");
   c1->cd(2*ncent-icent)->RedrawAxis();
   if(!doMC){
     if(icent==1) drawText("PbPb",0.3,0.9);
     else drawText("PbPb",0.22,0.9);
   }else drawText("PYTHIA+HYDJET,reco",0.22,0.9);

 } 
  
 for(int icent=0;icent<ncent;icent++){
 c1->cd(ncent-icent);
 empty->Draw();
  for(int ipt=npt-2;ipt>=0;ipt--){ 
   hmpt_vs_alpha_ref[ipt][icent]->Draw("same");
   hmpt_vs_alpha_ref[ipt][icent]->Draw("same hist");
  }
  hmpt_vs_alpha_ref[npt-1][icent]->Draw("same");

  
  if(icent==0)leg2->AddEntry(hmpt_vs_alpha_ref[npt-1][icent],">0.5","p");
  hmpt_vs_alpha_ref[npt-1][icent]->Draw("same");
  drawText(Form("%d-%d %%",(int)(0.5*cent_min[icent]),(int)(0.5*cent_max[icent])),0.22,0.85);
  drawText(Form("%d-%d %%",(int)(0.5*cent_min[icent]),(int)(0.5*cent_max[icent])),0.22,0.85);
   zeroLine_p->Draw("same");
   c1->cd(ncent-icent)->RedrawAxis();
   if(!doMC){
     if(icent==1) drawText("PbPb",0.3,0.9);
     else drawText("PbPb",0.22,0.9);
   }else drawText("PYTHIA+HYDJET,gen",0.22,0.9);

 } 
  
 c1->SaveAs(Form("systematics/systematics_eta%d_doMC%d_%s_doIntegrate%d_doGenJet%d_kind%d_Aj%d_%d_ncent%d.png",(int)(etadijet*10),doMC,svariable[index_var].Data(),doIntegrate,doGenJet,gensigbkgreco,(int)(Ajmin[iAj]*10),(int)(Ajmax[iAj]*10),ncent));
 c1->SaveAs(Form("systematics/systematics_eta%d_doMC%d_%s_doIntegrate%d_doGenJet%d_kind%d_Aj%d_%d_ncent%d.pdf",(int)(etadijet*10),doMC,svariable[index_var].Data(),doIntegrate,doGenJet,gensigbkgreco,(int)(Ajmin[iAj]*10),(int)(Ajmax[iAj]*10),ncent));

 TFile *outf = new TFile(Form("systematics/outf_Aj%d.root",iAj),"recreate");
 for(int icent=0;icent<ncent;icent++){ 
  hmpt_vs_alpha_diff[npt-1][icent]->Write();
 }
 outf->Close();
 
}