#include <iostream>
#include <vector>
#include <cmath>
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
		 
		 if(i==0 && j==1){
          pad[i][j] = new TPad(padName.Data(),padName.Data(),Xlow[i],Ylow[j],Xup[i]*0.837,Yup[j]*0.935);
		 }
		 else if(i==0 && j==0){
		  pad[i][j] = new TPad(padName.Data(),padName.Data(),Xlow[i],Ylow[j]*0.935,Xup[i],Yup[j]);
		 }
		 else if(i==1 && j==1){
		  pad[i][j] = new TPad(padName.Data(),padName.Data(),0.837*Xlow[i],Ylow[j],Xup[i],Yup[j]);
		 }
         else pad[i][j] = new TPad(padName.Data(),padName.Data(),
            Xlow[i],Ylow[j],Xup[i],Yup[j]);
         if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
		 else if(i==1 && j==1) pad[i][j]->SetLeftMargin(0.57*PadWidth);
         else pad[i][j]->SetLeftMargin(0);

         if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
		   else pad[i][j]->SetRightMargin(0);

         if(j==0){
          if(i==0)pad[i][j]->SetTopMargin(edge);
          else pad[i][j]->SetTopMargin(edge*0.9);
         }else pad[i][j]->SetTopMargin(0);

         if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
		    else if(i==0 && j==0) pad[i][j]->SetBottomMargin(0.17*PadHeight);
         else pad[i][j]->SetBottomMargin(0);

         pad[i][j]->Draw();
         pad[i][j]->cd();
         pad[i][j]->SetNumber(columns*j+i+1);
      }
   }
   pad[0][0]->cd();

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
 bool doGenJet=true;
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
 double systematics[] = { 0.5,1};
 double systematics_firstbin[] = { 3, 3};
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
 double mpt_ref[npt][nalpha+1];
 double mpterr_ref[npt][nalpha+1];
 double mptsum_ref[nalpha+1];
 double mptsumpos_ref[nalpha+1];
 double mptsumneg_ref[nalpha+1];
 TH1D* h_ref[npt][nalpha+1];
 TH1D* h_ref_stat[npt][nalpha+1];
 TH1D *hsum_ref[nalpha+1];
 TH1D * hmpt_vs_alpha_ref[npt];
 
 double mpt_cum;
 double mpterr_cum;  
 TH1D* h_cum[nalpha+1][ncent];
 TH1D* h_stat_cum[nalpha+1][ncent];
 TH1D * hmpt_vs_alpha_cum[ncent];

 double mpt_ref_cum[nalpha+1];
 double mpterr_ref_cum[nalpha+1];
 TH1D* h_ref_cum[nalpha+1];
 TH1D* h_ref_stat_cum[nalpha+1];
 TH1D *hsum_ref_cum[nalpha+1];
 TH1D * hmpt_vs_alpha_ref_cum;
 
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
   
 TLegend *leg,*leg2, *leg3,*leg4;
 leg= new TLegend(0.45,0.1,0.99,0.4);
 leg->SetBorderSize(0); 
 leg->SetFillStyle(0);
 leg->SetTextFont(43);
 leg->SetTextSize(18);
 leg3 = new TLegend(0.5,0.1,0.95,0.4);
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
   if(ipt==0 || ipt==npt-1) f[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_PbPb_20140428_lastptbinoldcorr//full_ntuple_HIRun2011-14Mar2014-v2-6lumi-jet80-forest-v4-merged_pt%d_%d_akVs3Calo_v2.root",(int)ptmin[ipt],(int)ptmax[ipt])); 
   else f[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_PbPb_20140428//full_ntuple_HIRun2011-14Mar2014-v2-6lumi-jet80-forest-v4-merged_pt%d_%d_akVs3Calo_v2.root",(int)ptmin[ipt],(int)ptmax[ipt])); 
   f_ref[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_pp_20140428//full_ntuple_HiForest_pp_Jet80_v8_PP2013_HiForest_PromptReco_JsonPP_Jet80_PPReco_merged_forest_0_pt%d_%d_ak3Calo.root",(int)ptmin[ipt],(int)ptmax[ipt])); 
  }
  
  if(doMC){
   f[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_MC_20140419/full_ntuple_HydjetDrum_Pyquen_Dijet_FOREST_Track8_Jet24_FixedPtHatJES_v0_pt%d_%d_akVs3Calo_doGenJet%d.root",(int)ptmin[ipt],(int)ptmax[ipt],doGenJet)); 
   f_ref[ipt] = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_PYTHIA_20140423/full_ntuple_HiForest_pt80_PYTHIA_ppReco_JECv85_merged_forest_0_pt%d_%d_ak3Calo_doGenJet%d.root",(int)ptmin[ipt],(int)ptmax[ipt],doGenJet)); 
  }
 
  // cout<<Form("nt_%s%s_%s%s",smpt[domp].Data(),skind[gensigbkgreco].Data(),svariable[index_var].Data(),strave.Data())<<endl; 

  f[ipt]->cd();
  t_jet[ipt]=(TTree*)f[ipt]->Get("nt_jet");
  f[ipt]->cd();
  if(index_var==0) t[ipt]=(TTree*)f[ipt]->Get(Form("nt_%s%s%s",smpt[domp].Data(),skind[gensigbkgreco].Data(),strave.Data()));
  else t[ipt]=(TTree*)f[ipt]->Get(Form("nt_%s%s_%s%s",smpt[domp].Data(),skind[gensigbkgreco].Data(),svariable[index_var].Data(),strave.Data()));
  t[ipt]->AddFriend(t_jet[ipt]);

 // cout<<Form("nt_%s%s_%s%s",smpt[domp].Data(),skind[gensigbkgreco].Data(),svariable[index_var].Data(),strave.Data())<<endl; 

  // cout<<Form("nt_%s%s%s",smpt[domp].Data(),skind[gensigbkgreco].Data(),strave.Data())<<endl;
  if(index_var==0) t_ref[ipt]=(TTree*)f_ref[ipt]->Get(Form("nt_%s%s%s",smpt[domp].Data(),skind[gensigbkgreco].Data(),strave.Data()));
  else t_ref[ipt]=(TTree*)f_ref[ipt]->Get(Form("nt_%s%s_%s%s",smpt[domp].Data(),skind[gensigbkgreco].Data(),svariable[index_var].Data(),strave.Data()));
  t_jet_ref[ipt]=(TTree*)f_ref[ipt]->Get("nt_jet");
  t_ref[ipt]->AddFriend(t_jet_ref[ipt]);
  
  hmpt_vs_alpha_ref[ipt]= new TH1D(Form("hmpt_ref_%d",ipt),"",nalpha/2+1,frac);
  if(ipt==npt-1)hmpt_vs_alpha_ref_cum= new TH1D(Form("hmpt_ref_cum_%d",ipt),"",nalpha/2+1,frac);
   
  for(int ialpha=0;ialpha<nalpha/2+1;ialpha++){
   if(ipt==npt-1){
    if(ialpha<nalpha/2)s[ipt]=Form("-%s_%d",svariable[index_var].Data(),2*ialpha+1);
    if(ialpha<nalpha/2)sintegrate_stat[ialpha][ipt]=s[ipt];
    else  sintegrate_stat[ialpha][ipt]=Form("-%s%s%s",smpt[domp].Data(),skind[gensigbkgreco].Data(),strave.Data());
	
    if(ialpha<nalpha/2 && ialpha>0) sintegrate[ialpha][ipt]=Form("-%s_%d+%s_%d",svariable[index_var].Data(),2*ialpha+1,svariable[index_var].Data(),2*(ialpha-1)+1);
    else if(ialpha==0)sintegrate[ialpha][ipt]=Form("-%s_%d",svariable[index_var].Data(),ialpha+1);
    else  sintegrate[ialpha][ipt]=Form("-%s%s%s+%s_%d",smpt[domp].Data(),skind[gensigbkgreco].Data(),strave.Data(),svariable[index_var].Data(),nalpha-1);
   }else{
    if(ialpha<nalpha/2 && ialpha>0) sintegrate[ialpha][ipt]=Form("-%s_%d+%s_%d",svariable[index_var].Data(),2*ialpha+1,svariable[index_var].Data(),2*(ialpha-1)+1);
    else if(ialpha==0)sintegrate[ialpha][ipt]=Form("-%s_%d",svariable[index_var].Data(),ialpha+1);
    else  sintegrate[ialpha][ipt]=Form("-%s%s%s+%s_%d",smpt[domp].Data(),skind[gensigbkgreco].Data(),strave.Data(),svariable[index_var].Data(),nalpha-1);
   }
   h_ref[ipt][ialpha]=new TH1D(Form("h_ref_%d_%d",ialpha,ipt),"",4000,-2000,2000);
   t_ref[ipt]->Draw(Form("%s>>h_ref_%d_%d",sintegrate[ialpha][ipt].Data(),ialpha,ipt),Form("pt1>%d && pt2>%d && abs(eta1)<%.1f && abs(eta2)<%.1f && dphi>(5*TMath::Pi()/6)&& ((pt1-pt2)/(pt2+pt1))>=%.2f && ((pt1-pt2)/(pt2+pt1))<%.2f && abs(vz)<15",jetpt1,jetpt2,etadijet,etadijet,Ajmin[iAj],Ajmax[iAj]),"");
   
   mpt_ref[ipt][ialpha]=h_ref[ipt][ialpha]->GetMean();
   mpterr_ref[ipt][ialpha]=h_ref[ipt][ialpha]->GetMeanError();
   if(ipt==npt-1){
    h_ref_stat[ipt][ialpha]=new TH1D(Form("h_ref_stat_%d_%d",ialpha,ipt),"",4000,-2000,2000);
    t_ref[ipt]->Draw(Form("%s>>h_ref_stat_%d_%d",sintegrate_stat[ialpha][ipt].Data(),ialpha,ipt),Form("pt1>%d && pt2>%d && abs(eta1)<%.1f && abs(eta2)<%.1f && dphi>(5*TMath::Pi()/6)&& ((pt1-pt2)/(pt2+pt1))>=%.2f && ((pt1-pt2)/(pt2+pt1))<%.2f && abs(vz)<15",jetpt1,jetpt2,etadijet,etadijet,Ajmin[iAj],Ajmax[iAj]),"");
    mpt_ref_cum[ialpha]=h_ref_stat[ipt][ialpha]->GetMean();
    mpterr_ref_cum[ialpha]=h_ref[ipt][ialpha]->GetMeanError();
   }
   
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
     hmpt_vs_alpha_ref[ipt]->SetMarkerStyle(25);
     hmpt_vs_alpha_ref_cum->SetBinContent(ialpha+1, mpt_ref_cum[ialpha]);
     hmpt_vs_alpha_ref_cum->SetBinError(ialpha+1, mpterr_ref_cum[ialpha]);
   }
  }
  
  for(int icent=0;icent<ncent;icent++){
  hmpt_vs_alpha[ipt][icent]= new TH1D(Form("hmpt_%d_%d",ipt,icent),"",nalpha/2+1,frac);
  hmpt_vs_alpha_diff[ipt][icent]= new TH1D(Form("hmpt_diff_%d_%d",ipt,icent),"",nalpha/2+1,frac);
  if(ipt==npt-1){
   hmpt_vs_alpha_cum[icent]= new TH1D(Form("hmpt_ref_%d_%d",ipt,icent),"",nalpha/2+1,frac);
  } 
  for(int ialpha=0;ialpha<nalpha/2+1;ialpha++ ){
   h[ipt][ialpha][icent]=new TH1D(Form("h_%d_%d_%d",ialpha,ipt,icent),"",4000,-2000,2000);
   h_stat[ipt][ialpha][icent]=new TH1D(Form("h_stat_%d_%d_%d",ialpha,ipt,icent),"",4000,-2000,2000);
   if(!doMC)t[ipt]->Draw(Form("%s>>h_%d_%d_%d",sintegrate[ialpha][ipt].Data(),ialpha,ipt,icent),Form("pt1>%d && pt2>%d && abs(eta1)<%.1f && abs(eta2)<%.1f && dphi>(5*TMath::Pi()/6) && cent>=%d && cent<%d && ((pt1-pt2)/(pt2+pt1))>=%.2f && ((pt1-pt2)/(pt2+pt1))<%.2f && abs(vz)<15",jetpt1,jetpt2,etadijet,etadijet, cent_min[icent],cent_max[icent],Ajmin[iAj],Ajmax[iAj]),"");
   else t[ipt]->Draw(Form("%s>>h_%d_%d_%d",sintegrate[ialpha][ipt].Data(),ialpha,ipt,icent),Form("weight*cent_weight*(pt1>%d && pt2>%d && abs(eta1)<%.1f && abs(eta2)<%.1f && dphi>(5*TMath::Pi()/6) && cent>=%d && cent<%d && ((pt1-pt2)/(pt2+pt1))>=%.2f && ((pt1-pt2)/(pt2+pt1))<%.2f && abs(vz)<15)",jetpt1,jetpt2,etadijet,etadijet, cent_min[icent],cent_max[icent],Ajmin[iAj],Ajmax[iAj]),"");
    mpt=h[ipt][ialpha][icent]->GetMean(); 
    mpterr=h[ipt][ialpha][icent]->GetMeanError();
   
   mpt=h[ipt][ialpha][icent]->GetMean(); 
   mpterr=h[ipt][ialpha][icent]->GetMeanError();
   
   if(ipt==npt-1){
    if(!doMC)t[ipt]->Draw(Form("%s>>h_stat_%d_%d_%d",sintegrate_stat[ialpha][ipt].Data(),ialpha,ipt,icent),Form("pt1>%d && pt2>%d && abs(eta1)<%.1f && abs(eta2)<%.1f && dphi>(5*TMath::Pi()/6) && cent>=%d && cent<%d && ((pt1-pt2)/(pt2+pt1))>=%.2f && ((pt1-pt2)/(pt2+pt1))<%.2f && abs(vz)<15",jetpt1,jetpt2,etadijet,etadijet, cent_min[icent],cent_max[icent],Ajmin[iAj],Ajmax[iAj]),"");
    else t[ipt]->Draw(Form("%s>>h_stat_%d_%d_%d",sintegrate_stat[ialpha][ipt].Data(),ialpha,ipt,icent),Form("weight*cent_weight*(pt1>%d && pt2>%d && abs(eta1)<%.1f && abs(eta2)<%.1f && dphi>(5*TMath::Pi()/6) && cent>=%d && cent<%d && ((pt1-pt2)/(pt2+pt1))>=%.2f && ((pt1-pt2)/(pt2+pt1))<%.2f && abs(vz)<15)",jetpt1,jetpt2,etadijet,etadijet, cent_min[icent],cent_max[icent],Ajmin[iAj],Ajmax[iAj]),""); 
    mpterr=h[ipt][ialpha][icent]->GetMeanError();
    mpterr_cum=h[ipt][ialpha][icent]->GetMeanError();
    mpt_cum=h_stat[ipt][ialpha][icent]->GetMean();
   }
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
     hmpt_vs_alpha_cum[icent]->SetBinContent(ialpha+1,mpt_cum);
     hmpt_vs_alpha_cum[icent]->SetBinError(ialpha+1,mpterr_cum);
   }   
   hmpt_vs_alpha[ipt][icent]->SetMarkerStyle(5);

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
 
 TFile *syst_PbPb = new TFile(Form("systematics/outf_Aj%d.root",iAj));
 TFile *syst_pp = new TFile(Form("systematics_pp/outf_Aj%d.root",iAj));
 // TFile *syst_PbPb_pp = new TFile(Form("systematics/outf_Aj%d.root",iAj),"recreate");

 TH1D * hsyst[ncent];
 TH1D * hsyst_PbPb_pp[ncent];
 TH1D * hsyst_pp=(TH1D *)syst_pp->Get(Form("hmpt_diff_%d_%d",npt-1,0));
 
 for(int icent=0;icent<ncent;icent++){
  hsyst[icent]=(TH1D *)syst_PbPb->Get(Form("hmpt_diff_%d_%d",npt-1,icent));
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
  
  TCanvas *c1 = new TCanvas("c1","",(ncent+1)*300,700);
  makeMultiPanelCanvas(c1,ncent+1,2,0.0,0.0,0.2,0.2,0.02);
 TH1D * empty=new TH1D("empty",Form(";%s;<#slash{p}_{T}^{#parallel}> (GeV/c)",axistitle[index_var].Data()),nalpha/2+1,frac);
 TH1D * empty2=new TH1D("empty2",Form(";%s;<#slash{p}_{T}^{#parallel}> (GeV/c)",axistitle[index_var].Data()),nalpha/2+1,frac);
 empty->Fill(0.5,1000); 
 empty2->Fill(0.5,1000); 
 if(doIntegrate){
  if(index_var==0){
   empty->SetMaximum(20); 
   empty->SetMinimum(-40); 
  
  }else{
   empty->SetMaximum(20); 
   empty->SetMinimum(-40); 
  }
 }else{
  empty->SetMaximum(20); 
  empty->SetMinimum(-40); 
  // empty->SetMaximum(10); 
  // empty->SetMinimum(-10); 
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
   
   
 c1->cd(ncent+2); 
 leg2->Draw("same");
 leg->Draw("same");

 if(Ajmin[iAj]==0 && Ajmax[iAj]==1) drawText(Form("|#eta_{1}|,|#eta_{2}|<%.2f",etadijet),0.05,0.54);
 else if(Ajmin[iAj]==0) drawText(Form("|#eta_{1}|,|#eta_{2}|<%.2f, A_{J} < %.2f",etadijet,Ajmax[iAj]),0.05,0.54);
 else if(Ajmax[iAj]==1)drawText(Form("|#eta_{1}|,|#eta_{2}|<%.2f,  A_{J} > %.2f",etadijet,Ajmin[iAj]),0.05,0.54);
 else drawText(Form("|#eta_{1}|,|#eta_{2}|<%.2f, %.2f < A_{J} < %.0f",etadijet,Ajmin[iAj],Ajmax[iAj]),0.05,0.54);
 drawText(Form("p_{T,1}>%d, p_{T,2}>%d GeV/c",jetpt1,jetpt2),0.05,0.9);
 drawText(Form("|#eta_{jet}|<2, #Delta#phi_{1,2}>5#pi/6"),0.05,0.78);
 drawText(Form("anti-kt Vs Calo R=0.3"),0.05,0.66);
 drawText(Form("|#eta_{trk}|<2.4"),0.05,0.42);
 if(!doMC){
   // drawText(Form("PbPb #sqrt{s_{NN}}=2.76 TeV 150/#mub"),0.05,0.66);
   // drawText(Form("pp #sqrt{s_{NN}}=2.76 TeV 5.3/pb"),0.05,0.54);
 }
 if(doMC){
  if(!doGenJet)drawText(Form("akVs3Calo, ak3Calo"),0.05,0.66);
  else drawText(Form("gen jet"),0.05,0.66);
  drawText(Form("%s",spart[gensigbkgreco].Data()),0.05,0.54);
 } 
 for(int icent=0;icent<ncent;icent++){ 
 c1->cd((2*(ncent+1))-icent);
 
 
 if(doIntegrate) empty->Draw();
   else empty2->Draw();
  for(int ipt=npt-2;ipt>=0;ipt--){ 
   hmpt_vs_alpha_diff[ipt][icent]->Draw("same");
   hmpt_vs_alpha_diff[ipt][icent]->Draw("same hist");
  }
  if(index_var==3){
  // if(0){
  for(int ialpha = 0 ;ialpha<hmpt_vs_alpha_diff[npt-1][icent]->GetNbinsX();ialpha++){
   double  y = hmpt_vs_alpha_diff[npt-1][icent]->GetBinContent(ialpha+1);
   double yerr;
   if(ialpha==0) yerr = 3;
   else  yerr = systematics[icent];
   double x=hmpt_vs_alpha_diff[npt-1][icent]->GetBinCenter(ialpha+1);
   
    TLine * l1= new TLine(x-0.035,y-yerr,x+0.035,y-yerr);
    TLine * l2= new TLine(x-0.035,y-yerr,x+0.035,y-yerr);
    TLine * l3= new TLine(x-0.035,y+yerr,x+0.035,y+yerr);
    TLine * l4= new TLine(x-0.035,y+yerr,x+0.035,y+yerr);
    TLine * l5= new TLine(x-0.035,y-yerr,x-0.035,y-yerr+0.25);
    TLine * l6= new TLine(x-0.035,y+yerr,x-0.035,y+yerr-0.25);
    TLine * l7= new TLine(x+0.035,y-yerr,x+0.035,y-yerr+0.25);
    TLine * l8= new TLine(x+0.035,y+yerr,x+0.035,y+yerr-0.25);
  
    l1->Draw("same");
    l2->Draw("same");
    l3->Draw("same");
    l4->Draw("same");
    l5->Draw("same");
    l6->Draw("same");
    l7->Draw("same");
    l8->Draw("same");
   }
  }
  hmpt_vs_alpha_diff[npt-1][icent]->Draw("same");
  zeroLine_p->Draw("same");
  if(icent==1)drawText(Form("%d-%d %%",(int)(0.5*cent_min[icent]),(int)(0.5*cent_max[icent])),0.2,0.85);
  else drawText(Form("%d-%d %%",(int)(0.5*cent_min[icent]),(int)(0.5*cent_max[icent])),0.05,0.85);
   c1->cd((2*(ncent+1))-icent)->RedrawAxis();
   if(!doMC){
    if(icent==1)  drawText("PbPb-pp",0.2,0.9);
    else  drawText("PbPb-pp",0.05,0.9);
   }else drawText("(PYT.+HYD.)-PYT.",0.22,0.9);
   if(icent!=ncent-1)drawPatch(0.0,0.0,0.05,0.19);
   if(icent!=0)drawPatch(0.95,0.0,1,0.19);
   if(icent!=0)drawPatch(0.0,0.95,0.16,1);
 }
 // leg4->Draw("same");
 c1->cd(1);
 empty->Draw();
 leg3->AddEntry(hmpt_vs_alpha_ref[npt-1],"pp","p");
 leg3->AddEntry(hmpt_vs_alpha[npt-1][0],"PbPb","p");
 leg3->AddEntry(hmpt_vs_alpha_diff[npt-1][0],"PbPb-pp","p");

 leg3->Draw("same");
 for(int ipt=npt-2;ipt>=0;ipt--){
  hmpt_vs_alpha_ref[ipt]->Draw("same");
  hmpt_vs_alpha_ref[ipt]->Draw("same hist");
 }
 hmpt_vs_alpha_ref[npt-1]->Draw("same");
 hmpt_vs_alpha_ref_cum->SetMarkerSize(0);
 hmpt_vs_alpha_ref_cum->SetLineStyle(2);
 hmpt_vs_alpha_ref_cum->Draw("same C HIST");
 if(index_var==3){
  for(int ialpha = 0 ;ialpha<hmpt_vs_alpha_ref[npt-1]->GetNbinsX();ialpha++){
   double  y = hmpt_vs_alpha_ref[npt-1]->GetBinContent(ialpha+1);
   double yerr = fabs(hsyst_pp->GetBinContent(ialpha+1));
   double x=hmpt_vs_alpha_ref[npt-1]->GetBinCenter(ialpha+1);
   
    TLine * l1= new TLine(x-0.035,y-yerr,x+0.035,y-yerr);
    TLine * l2= new TLine(x-0.035,y-yerr,x+0.035,y-yerr);
    TLine * l3= new TLine(x-0.035,y+yerr,x+0.035,y+yerr);
    TLine * l4= new TLine(x-0.035,y+yerr,x+0.035,y+yerr);
    TLine * l5= new TLine(x-0.035,y-yerr,x-0.035,y-yerr+0.25);
    TLine * l6= new TLine(x-0.035,y+yerr,x-0.035,y+yerr-0.25);
    TLine * l7= new TLine(x+0.035,y-yerr,x+0.035,y-yerr+0.25);
    TLine * l8= new TLine(x+0.035,y+yerr,x+0.035,y+yerr-0.25);
  
    l1->Draw("same");
    l2->Draw("same");
    l3->Draw("same");
    l4->Draw("same");
    l5->Draw("same");
    l6->Draw("same");
    l7->Draw("same");
    l8->Draw("same");
   }
  }
 
 zeroLine_p->Draw("same");
 c1->cd(1)->RedrawAxis();
 if(!doMC){
  drawText("pp #sqrt{s_{NN}}=2.76 TeV 5.3 pb^{-1}",0.3,0.9);
  drawText("CMS Preliminary",0.5,0.84);
 }else drawText("PYTHIA",0.3,0.9);

 for(int icent=0;icent<ncent;icent++){
 c1->cd(ncent+1-icent);
 empty->Draw();
  for(int ipt=npt-2;ipt>=0;ipt--){ 
   if(icent==0){
    if(ipt>2) leg2->AddEntry(hmpt_vs_alpha[ipt][icent],Form("%.1f - %.1f",ptmin[ipt],ptmax[ipt]),"f");
    else leg->AddEntry(hmpt_vs_alpha[ipt][icent],Form("%.1f - %.1f",ptmin[ipt],ptmax[ipt]),"f");
   }
   hmpt_vs_alpha[ipt][icent]->Draw("same");
   hmpt_vs_alpha[ipt][icent]->Draw("same hist");
  }
  hmpt_vs_alpha[npt-1][icent]->Draw("same");
  hmpt_vs_alpha_cum[icent]->SetMarkerSize(0);
  hmpt_vs_alpha_cum[icent]->Draw("same C HIST");
  hmpt_vs_alpha_ref_cum->SetLineStyle(2);
  hmpt_vs_alpha_ref_cum->SetMarkerSize(0);
  hmpt_vs_alpha_ref_cum->Draw("same C HIST");
  if(index_var==3){
   for(int ialpha = 0 ;ialpha<hmpt_vs_alpha[npt-1][icent]->GetNbinsX();ialpha++){
    double  y = hmpt_vs_alpha[npt-1][icent]->GetBinContent(ialpha+1);
    double yerr = fabs(hsyst[icent]->GetBinContent(ialpha+1));
    double x=hmpt_vs_alpha[npt-1][icent]->GetBinCenter(ialpha+1);
   
     TLine * l1= new TLine(x-0.035,y-yerr,x+0.035,y-yerr);
     TLine * l2= new TLine(x-0.035,y-yerr,x+0.035,y-yerr);
     TLine * l3= new TLine(x-0.035,y+yerr,x+0.035,y+yerr);
     TLine * l4= new TLine(x-0.035,y+yerr,x+0.035,y+yerr);
     TLine * l5= new TLine(x-0.035,y-yerr,x-0.035,y-yerr+TMath::Min(0.25,yerr));
     TLine * l6= new TLine(x-0.035,y+yerr,x-0.035,y+yerr-TMath::Min(0.25,yerr));
     TLine * l7= new TLine(x+0.035,y-yerr,x+0.035,y-yerr+TMath::Min(0.25,yerr));
     TLine * l8= new TLine(x+0.035,y+yerr,x+0.035,y+yerr-TMath::Min(0.25,yerr));
  
     l1->Draw("same");
     l2->Draw("same");
     l3->Draw("same");
     l4->Draw("same");
     l5->Draw("same");
     l6->Draw("same");
     l7->Draw("same");
     l8->Draw("same");
    }
   }
  if(icent==0)leg2->AddEntry(hmpt_vs_alpha[npt-1][icent],">0.5","p");
   drawText(Form("%d-%d %%",(int)(0.5*cent_min[icent]),(int)(0.5*cent_max[icent])),0.05,0.84);
   zeroLine_p->Draw("same");
   c1->cd(ncent+1-icent)->RedrawAxis();
   if(!doMC){
     drawText("PbPb #sqrt{s_{NN}}=2.76 TeV 150 #mub^{-1}",0.05,0.9);
     if(icent==1)drawText(" ",0.4,0.82);
   }else drawText("PYTHIA+HYDJET",0.22,0.9);

 } 
  
 c1->SaveAs(Form("cumulative_eta%d_doMC%d_%s_doIntegrate%d_doGenJet%d_kind%d_Aj%d_%d_ncent%d.png",(int)(etadijet*10),doMC,svariable[index_var].Data(),doIntegrate,doGenJet,gensigbkgreco,(int)(Ajmin[iAj]*10),(int)(Ajmax[iAj]*10),ncent));
 c1->SaveAs(Form("cumulative_eta%d_doMC%d_%s_doIntegrate%d_doGenJet%d_kind%d_Aj%d_%d_ncent%d.pdf",(int)(etadijet*10),doMC,svariable[index_var].Data(),doIntegrate,doGenJet,gensigbkgreco,(int)(Ajmin[iAj]*10),(int)(Ajmax[iAj]*10),ncent));
 

 
 
 
 // [7:09:13 AM] Yen-Jie Lee: 
}