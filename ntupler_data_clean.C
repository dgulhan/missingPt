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
#include "TF1.h"

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
#include "ntupler/trackTree.C"

Float_t getFlippedPhi(Float_t inPhi)
{
  Float_t outPhi;

  if(TMath::Abs(inPhi) > TMath::Pi()){
    // std::cout << "getFlippedPhi: inPhi is outside accepted range, return -10" << std::endl;
    return -10;
  }
  else if(inPhi > 0)
    outPhi = inPhi - TMath::Pi();
  else
    outPhi = inPhi + TMath::Pi();

  return outPhi;
}

Bool_t sameSign(Float_t num1, Float_t num2){
  if((num1 > 0 && num2 > 0) || (num1 < 0 && num2 < 0)) return true;

  return false;
}

Float_t getAvePhi(Float_t inLeadPhi, Float_t inSubLeadPhi)
{
  Float_t flipPhi = getFlippedPhi(inSubLeadPhi);
  Float_t avePhi;

  if(sameSign(inLeadPhi, flipPhi) || (TMath::Abs(inLeadPhi) < TMath::Pi()/2 && TMath::Abs(flipPhi) < TMath::Pi()/2))
    avePhi = (flipPhi + inLeadPhi)/2;
  else if(TMath::Abs(inLeadPhi) > TMath::Pi()/2 && TMath::Abs(flipPhi) > TMath::Pi()/2){
    avePhi = (flipPhi + inLeadPhi)/2;
    if(avePhi > 0)
      avePhi = TMath::Pi() - avePhi;
    else
      avePhi = -TMath::Pi() - avePhi;
  }
  else{
    avePhi = 0.;
  }

  return avePhi;
}


void ntupler_data_clean(double ptmin_trk=8,double ptmax_trk=300){ 
 TH1D::SetDefaultSumw2();
 TString algo="akVs3Calo"; 
 bool doSaveTrackInfo=false;
 //input file  
 // cout<<"ptmin= "<<ptmin_trk<<" ptmax= "<<ptmax_trk<<endl;
 // TString directory="root://eoscms//eos/cms/store/group/phys_heavyions/dgulhan";
 // TString infname="HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_0"; 
  
 TString directory="/d100/JetSample/hiForest";
 TString infname="hiForest_Jet80or95_GR_R_53_LV6_12Mar2014_0000CET_Track8_Jet21_0"; 
    
 trackTree * ftrk = new trackTree(Form("%s/%s.root",directory.Data(),infname.Data()));
 HiTree * fhi = new HiTree(Form("%s/%s.root",directory.Data(),infname.Data()));
 skimTree * fskim = new skimTree(Form("%s/%s.root",directory.Data(),infname.Data()));
 t * fjet = new t(Form("%s/%s.root",directory.Data(),infname.Data()),algo.Data());
 hi * fgen = new hi(Form("%s/%s.root",directory.Data(),infname.Data()));
 // hi * fgen = new hi(Form("%s/%s.root",directory.Data(),infname.Data()));

 
 //getting unfolding profile 
 TFile * f_unfold=new TFile("binbybin_pt.root");
 TH1D *hweight = (TH1D*)f_unfold->Get("hweight");
 TH1D *h_res = (TH1D*)f_unfold->Get("h_res");
 TF1 * fgaus=new TF1("fgaus","gaus(0)",-20,20);
 fgaus->SetParameters(1,0,1);
 
 TFile *f_cent = new TFile("/afs/cern.ch/user/d/dgulhan/workDir/JetTrack/centralityWeight/f_cent_weight.root");
 
 //pt bins for track efficiency correction
  int npt=14; 
 double ptmin[]={0.4,0.4,0.4,0.4,0.4, 1, 1, 1,  1,  1, 3, 3,  3,  8};
 double ptmax[]={  1,  1,  1,  1,  1, 3, 3, 3,  3,  3, 8, 8,  8,300};
 
 int cent_min[]={  0, 20, 40, 60,100, 0,20,40, 60,100, 0,20, 40,  0};
 int cent_max[]={ 20, 40, 60,100,200,20,40,60,100,200,20,40,200,200};
 

 double frac[]={0,0.01,0.02,0.03,0.04,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.925,0.95,0.96,0.97,0.98,0.99,1.01};
  double dR_upperbound[]={0,0.01,0.02,0.03,0.04,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,1,1.1,1.25,1.5,2,2.5,3,5};
  
 TFile *f_eff[npt];
 TProfile *p_eff_cent[npt]; 
 TProfile2D *p_eff_accept[npt]; 
 TProfile *p_eff_pt[npt]; 
 TProfile *p_eff_rmin[npt]; 
 for(int ipt=0; ipt<npt;ipt++){
   if(ipt<13)f_eff[ipt]= new TFile(Form("/afs/cern.ch/user/d/dgulhan/workDir/trackCorrection/akVs3Calo/eff/eff_pt%d_%d_cent%d_%d_step_cent4accept4pt4rmin3.root",(int)ptmin[ipt],(int)ptmax[ipt],(int)(0.5*cent_min[ipt]),(int)(0.5*cent_max[ipt])));
   else f_eff[ipt]= new TFile(Form("/afs/cern.ch/user/d/dgulhan/workDir/trackCorrection/akVs3Calo/eff/eff_pt%d_%d_cent%d_%d_step_cent3accept3pt3rmin3.root",(int)ptmin[ipt],(int)ptmax[ipt],(int)(0.5*cent_min[ipt]),(int)(0.5*cent_max[ipt])));
   p_eff_cent[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_cent");
   p_eff_pt[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_pt");
   p_eff_accept[ipt]=(TProfile2D*)f_eff[ipt]->Get("p_eff_acceptance");
   p_eff_rmin[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_rmin");
 }

 TFile *f_fake[npt];
 TProfile *p_fake_cent[npt]; 
 TProfile2D *p_fake_accept[npt]; 
 TProfile *p_fake_pt[npt]; 
 TProfile *p_fake_rmin[npt]; 
 for(int ipt=0; ipt<npt;ipt++){
   // if(ipt==0)f_fake[ipt]= new TFile(Form("/afs/cern.ch/work/d/dgulhan/trackFake/final_hists_akVs3Calo/fake_pt%d_%d_cent%d_%d_step_cent4accept4pt4rmin3_%s_dogenjet0.root",(int)ptmin[ipt],(int)ptmax[ipt],cent_min[ipt],cent_max[ipt],algo.Data()));
   if(ipt<4)f_fake[ipt]= new TFile(Form("/afs/cern.ch/user/d/dgulhan/workDir/trackFake/final_hists_akVs3Calo_20140312/fake_pt%d_%d_cent%d_%d_step_cent5accept5pt5rmin4_%s_dogenjet0.root",(int)ptmin[ipt],(int)ptmax[ipt],cent_min[ipt],cent_max[ipt],algo.Data()));
   else if(ipt<13)f_fake[ipt]= new TFile(Form("/afs/cern.ch/user/d/dgulhan/workDir/trackFake/final_hists_akVs3Calo_20140312/fake_pt%d_%d_cent%d_%d_step_cent4accept4pt4rmin3_%s_dogenjet0.root",(int)ptmin[ipt],(int)ptmax[ipt],cent_min[ipt],cent_max[ipt],algo.Data()));
   else f_fake[ipt]= new TFile(Form("/afs/cern.ch/user/d/dgulhan/workDir/trackFake/final_hists_akVs3Calo_20140312/fake_pt%d_%d_cent%d_%d_step_cent3accept3pt3rmin3_%s_dogenjet0.root",(int)ptmin[ipt],(int)ptmax[ipt],cent_min[ipt],cent_max[ipt],algo.Data()));
   p_fake_cent[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_cent");
   p_fake_pt[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_pt");
   p_fake_accept[ipt]=(TProfile2D*)f_fake[ipt]->Get("p_fake_acceptance");
   p_fake_rmin[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_rmin");
 }
 
 //output file and tree
 TFile *outf= new TFile(Form("ntuples_PbPb_mptonly_20140322/full_ntuple_%s_pt%d_%d_%s.root",infname.Data(),(int)ptmin_trk,(int)ptmax_trk,algo.Data()),"recreate");
 std::string partVars="pt:eta:phi:rmin:pNRec:smeared_pt:cent:matchedpt:eff:trackselect:pt1:eta1:phi1:pt2:eta2:phi2:pt3:eta3:phi3:dphi:ptratio:jetspz:jetspt:jetseta:jetsrap:jetspz_12:jetspt_12:jetsrap_12:alpha:mpt_track:vz";
 std::string trackVars="trackselect:eff:fake:trkfake:trkstatus:weight_unfold:pt:eta:phi:rmin:cent:cent_weight:pt1:eta1:phi1:pt2:eta2:phi2:pt3:eta3:phi3:dphi:ptratio:jetspz:jetspt:jetseta:jetsrap:jetspz_12:jetspt_12:jetsrap_12:alpha:mpt_track:vz";
 
 
 TNtuple *nt_track = new TNtuple("nt_track","",trackVars.data()); 
 
 std:string jetVars="cent:pt1:eta1:phi1:pt2:eta2:phi2:pt3:eta3:phi3:dphi:ptratio:jetspz:jetspt:jetseta:jetsrap:jetsphi:jetspz_12:jetspt_12:jetseta_12:jetsrap_12:jetsphi_12:trkrap:mpt_tracks:mpt_tracks_ls:mpt_tracks_uncorr:mpt_tracks_ave:mpt_tracks_uncorr_ave:hfp:hfm:Npassingtrk:vz";
 
 
 
 TNtuple *nt_jet = new TNtuple("nt_jet","",jetVars.data());

 std::string mptVars="alpha_0:alpha_1:alpha_2:alpha_3:alpha_4:alpha_5:alpha_6:alpha_7:alpha_8:alpha_9:alpha_10:alpha_11:alpha_12:alpha_13:alpha_14:alpha_15:alpha_16:alpha_17:alpha_18:alpha_19:alpha_20:alpha_21:alpha_22:alpha_23:alpha_24:alpha_25:alpha_26:alpha_27:alpha_28:alpha_29";
 TNtuple *nt_mpt_tracks=new TNtuple("nt_mpt_tracks","",mptVars.data());
 TNtuple *nt_mpt_tracks_uncorr=new TNtuple("nt_mpt_tracks_uncorr","",mptVars.data());
 
 std::string mptVars_dR="dR_0:dR_1:dR_2:dR_3:dR_4:dR_5:dR_6:dR_7:dR_8:dR_9:dR_10:dR_11:dR_12:dR_13:dR_14:dR_15:dR_16:dR_17:dR_18:dR_19:dR_20:dR_21:dR_22:dR_23:dR_24:dR_25:dR_26:dR_27:dR_28:dR_29:dR_30:dR_31:dR_32:dR_33:dR_34:dR_35";
 TNtuple *nt_mpt_tracks_dR=new TNtuple("nt_mpt_tracks_dR","",mptVars_dR.data());
 TNtuple *nt_mpt_tracks_uncorr_dR=new TNtuple("nt_mpt_tracks_uncorr_dR","",mptVars_dR.data());
 TNtuple *nt_mpt_tracks_dR_ave=new TNtuple("nt_mpt_tracks_dR_ave","",mptVars_dR.data());
 TNtuple *nt_mpt_tracks_uncorr_dR_ave=new TNtuple("nt_mpt_tracks_uncorr_dR_ave","",mptVars_dR.data());

 
 cout<<"2"<<endl;
 //loop over events
 int nentries = ftrk->GetEntriesFast();
 for(int jentry=0;jentry<nentries;jentry++){
 // for(int jentry=0;jentry<1;jentry++){
  if((jentry%1000)==0) std::cout<<jentry<<"/"<<nentries<<std::endl;
  fskim->GetEntry(jentry);
  ftrk->GetEntry(jentry);
  fhi->GetEntry(jentry);
  fjet->GetEntry(jentry);
  
  // cout<<"pcoll "<<fskim->pcollisionEventSelection<<endl;
  // cout<<"noisefilter "<<fskim->pHBHENoiseFilter<<endl;
   if(!(fskim->pcollisionEventSelection && fskim->pHBHENoiseFilter)) continue;
  
  // cout<<"passed selection"<<endl;
  float cent=fhi->hiBin;
  float vz = fhi->vz;
  float hfp = fhi->hiHFplusEta4;
  float hfm = fhi->hiHFminusEta4;
  float pt1=-99;
  float phi1=-99;
  float eta1=-99;
  float neutralMax1=-99;
  float trackMax1=-99;
  float trackN1=-99;
  float pt2=-99;
  float phi2=-99;
  float eta2=-99;
  float neutralMax2=-99;
  float trackMax2=-99;
  float trackN2=-99;
  float pt3=-99;
  float phi3=-99;
  float eta3=-99;
  float neutralMax3=-99;
  float trackMax3=-99;
  float trackN3=-99; 
  float dphi=-99;
  float ptratio=-99;
  std::vector<std::pair<double, std::pair<double,double> >  >jets;
  
  int njet=0;
  TLorentzVector vjets;
  for(int ijet=0;ijet<fjet->nref;ijet++){
   if(fabs(fjet->jteta[ijet])>2 ||fjet->jtpt[ijet]<30) continue;
   jets.push_back(std::make_pair(fjet->jtpt[ijet],std::make_pair(fjet->jteta[ijet],fjet->jtphi[ijet])));
   TLorentzVector vjet;
   vjet.SetPtEtaPhiM(fjet->jtpt[ijet],fjet->jteta[ijet],fjet->jtphi[ijet],0);
   vjets+=vjet;
   njet++;
  } 
 

  double jetspx=-99;
  double jetspy=-99; 
  double jetspz=-99;
  double jetspt=-99;
  double jetseta=-99; 
  double jetsrap=-99;
  double jetsphi=-99;
  TLorentzVector v1,v2,v3;

  TLorentzVector vls;
  std::sort(jets.begin(),jets.end());
  if(njet>0){
   jetspx =vjets.Px();
   jetspy=vjets.Py();
   jetspz=vjets.Pz();
   jetspt =vjets.Pt();
   jetseta= vjets.Eta();
   jetsrap =vjets.Rapidity();
   jetsphi =vjets.Phi();
   pt1= jets[njet-1].first;
   eta1= jets[njet-1].second.first;
   phi1= jets[njet-1].second.second;
   v1.SetPtEtaPhiM(pt1,eta1,phi1,0);
   vls+=v1;
   if(njet>1){
    pt2=jets[njet-2].first;
    eta2=jets[njet-2].second.first;
    phi2=jets[njet-2].second.second;
    v2.SetPtEtaPhiM(pt2,eta2,phi2,0);
    dphi=acos(cos(phi1-phi2));
    ptratio=pt2/pt1;
	  vls+=v2;

    if(njet>2){
     pt3=jets[njet-3].first;
     eta3=jets[njet-3].second.first;
     phi3=jets[njet-3].second.second;
     v3.SetPtEtaPhiM(pt2,eta2,phi2,0);
    }
   }
  }
  
  double phiave =getAvePhi(phi1,phi2);
  double jetspt_12=-99;
  double jetspx_12=-99;
  double jetspy_12=-99;
  double jetspz_12=-99;
  double jetseta_12=-99;
  double jetsphi_12=-99;
  double jetsrap_12=-99;
  
  if(njet>1){ 
    jetspt_12=vls.Pt();
    jetspx_12=vls.Px();
    jetspy_12=vls.Py();
    jetspz_12=vls.Pz();
    jetseta_12=vls.Eta();
    jetsphi_12=vls.Phi();
    jetsrap_12=vls.Rapidity();
  }
  
  
  float mpt_dR_tracks_uncorr_ave[35];
  float mpt_dR_tracks_ave[35];

  float mpt_tracks_ave=0;
  float mpt_tracks_uncorr_ave=0;
  float mpt_tracks=0;
  float mpt_tracks_ls=0;
  float mpt_tracks_uncorr=0;
  float mpt_alpha_tracks[30];  
  float mpt_alpha_tracks_uncorr[30];
  float mpt_dR_tracks[35];
  float mpt_dR_tracks_uncorr[35];
  
  
  for(int i=0;i<30;i++){
   mpt_alpha_tracks[i]=0;  
   mpt_alpha_tracks_uncorr[i]=0; 
  }
  for(int i=0;i<35;i++){
   mpt_dR_tracks[i]=0;  
   mpt_dR_tracks_uncorr[i]=0;  
      
   mpt_dR_tracks_ave[i]=0;  
   mpt_dR_tracks_uncorr_ave[i]=0;  
  }
  
  TLorentzVector vtracks;
  
  int Npassingtrk=0;
  for(int itrk=0;itrk<ftrk->nTrk;itrk++){
   float eta=ftrk->trkEta[itrk];
   float pt=ftrk->trkPt[itrk];
   float phi=ftrk->trkPhi[itrk];
   
   float trkfake=ftrk->trkFake[itrk];
   float trkstatus=ftrk->trkStatus[itrk];
   float eff_pt,eff_cent,eff_accept,eff_rmin;
   eff_pt=eff_cent=eff_accept=eff_rmin=1;
   float fake_pt,fake_cent,fake_accept,fake_rmin;
   fake_pt=fake_cent=fake_accept=fake_rmin=0;
   
   if(fabs(eta)>2.4) continue;
   if(pt<ptmin_trk || pt>ptmax_trk) continue; 
   float trackselect=(ftrk->highPurity[itrk] && fabs(ftrk->trkDxy1[itrk]/ftrk->trkDxyError1[itrk])<3.0 && fabs(ftrk->trkDz1[itrk]/ftrk->trkDzError1[itrk])<3 && (ftrk->trkPtError[itrk]/ftrk->trkPt[itrk])<0.1);
   
   float weight_unfold=hweight->GetBinContent(hweight->FindBin(pt));
   
   float rmin=100;
 
   //find rmin; 
   for(int ijet=0;ijet<fjet->nref;ijet++){
    if(fabs(fjet->jteta[ijet])>2 || fjet->jtpt[ijet]<30) continue;
    float r_reco=sqrt(pow(eta-fjet->jteta[ijet],2)+pow(acos(cos(phi-fjet->jtphi[ijet])),2));
    if(r_reco<rmin)rmin=r_reco;
   }
   
   for(int ipt=0;ipt<npt;ipt++){
    if(pt>=ptmin[ipt] && pt<ptmax[ipt] && cent>=cent_min[ipt] && cent<cent_max[ipt]){
      eff_pt=p_eff_pt[ipt]->GetBinContent(p_eff_pt[ipt]->FindBin(pt));
      eff_cent=p_eff_cent[ipt]->GetBinContent(p_eff_cent[ipt]->FindBin(cent));
      eff_accept=p_eff_accept[ipt]->GetBinContent(p_eff_accept[ipt]->GetXaxis()->FindBin(phi),p_eff_accept[ipt]->GetYaxis()->FindBin(eta));
      if(rmin<5)eff_rmin=p_eff_rmin[ipt]->GetBinContent(p_eff_rmin[ipt]->FindBin(rmin));
     }     
   } 
   
   for(int ipt=0;ipt<npt;ipt++){
    if(pt>=ptmin[ipt] && pt<ptmax[ipt] && cent>=cent_min[ipt] && cent<cent_max[ipt]){
      fake_pt=p_fake_pt[ipt]->GetBinContent(p_fake_pt[ipt]->FindBin(pt));
      fake_cent=p_fake_cent[ipt]->GetBinContent(p_fake_cent[ipt]->FindBin(cent));
      fake_accept=p_fake_accept[ipt]->GetBinContent(p_fake_accept[ipt]->GetXaxis()->FindBin(phi),p_fake_accept[ipt]->GetYaxis()->FindBin(eta));
      if(rmin<5) fake_rmin=p_fake_rmin[ipt]->GetBinContent(p_fake_rmin[ipt]->FindBin(rmin));
     }     
   }
  
   float eff=1;
   eff=eff_accept*eff_cent*eff_pt*eff_rmin;
   if(eff==0){
    if(pt>100)eff=0.8;
	else eff=1;
   }
   float fake=0;
   if(pt<100)fake=fake_accept+fake_cent+fake_pt+fake_rmin;
   
   TLorentzVector vp;
   vp.SetPtEtaPhiM(pt,eta,phi,0.13957018);
   vtracks+=vp;
   
   float alpha=-99;
   float alpha_sub=-99;
    
   alpha=vp.Angle(v1.Vect());
   alpha_sub=vp.Angle(v2.Vect());
   float mpt_track=pt*cos(phi-phi1);
   if(trackselect){
    mpt_tracks+=((1-fake)/eff)*pt*cos(phi-phi1);
    mpt_tracks_ls+=((1-fake)/eff)*pt*(cos(phi-phi1)-cos(phi-phi2));
    mpt_tracks_uncorr+=pt*cos(phi-phi1);
    mpt_tracks_ave+=((1-fake)/eff)*pt*cos(phi-phiave);
    mpt_tracks_uncorr_ave+=pt*cos(phi-phiave);
   }
   if(njet>1){
   
    float entry[]={trackselect,eff,fake,trkfake,trkstatus,weight_unfold,pt,eta,phi,rmin,cent,pt1,eta1,phi1,pt2,eta2,phi2,pt3,eta3,phi3,dphi,ptratio,jetspz,jetspt,jetseta,jetsrap,jetspz_12,jetspt_12,jetsrap_12,alpha,mpt_track,vz};
    if(doSaveTrackInfo) nt_track->Fill(entry);
    
	if(!trackselect) continue;
    Npassingtrk++;
    double mpt_alpha_t[30];    
    double mpt_dR_t[35];
    double mpt_dR_t_ave[35];
    
   for(int i=0;i<30;i++){
      mpt_alpha_t[i]=0;
	  if(alpha <(frac[i+1]*TMath::Pi()) && alpha >=(frac[i]*TMath::Pi())) mpt_alpha_t[i] =pt*cos(phi -phi1);
      mpt_alpha_tracks[i]+=((1-fake)/eff)*mpt_alpha_t[i];
      mpt_alpha_tracks_uncorr[i]+=mpt_alpha_t[i];
    }
    for(int i=0;i<35;i++){
      mpt_dR_t[i]=0;
      mpt_dR_t_ave[i]=0;
	  
	  double dR1=sqrt(pow(acos(cos(phi1-phi)),2)+pow(eta1-eta,2));
	  double dR2=sqrt(pow(acos(cos(phi2-phi)),2)+pow(eta2-eta,2));
	  
	  if(((dR1 < dR_upperbound[i+1]) && dR1 >=dR_upperbound[i]) || ((dR2 < dR_upperbound[i+1]) && dR2 >=dR_upperbound[i])){
 	   mpt_dR_t[i] =pt*cos(phi -phi1);
	   mpt_dR_t_ave[i]=pt*cos(phi-phiave);
	  }
    
      mpt_dR_tracks[i]+=((1-fake)/eff)*mpt_dR_t[i];
      mpt_dR_tracks_ave[i]+=((1-fake)/eff)*mpt_dR_t_ave[i];
      mpt_dR_tracks_uncorr[i]+=mpt_dR_t[i];
      mpt_dR_tracks_uncorr_ave[i]+=mpt_dR_t_ave[i];
	 }
   }
   //fill in the output tree
  } 
  
  double trkrap=vtracks.Rapidity();
  nt_mpt_tracks->Fill(mpt_alpha_tracks);
  nt_mpt_tracks_uncorr->Fill(mpt_alpha_tracks_uncorr);
  nt_mpt_tracks_dR->Fill(mpt_dR_tracks);
  nt_mpt_tracks_uncorr_dR->Fill(mpt_dR_tracks_uncorr);
  nt_mpt_tracks_dR_ave->Fill(mpt_dR_tracks_ave);
  nt_mpt_tracks_uncorr_dR_ave->Fill(mpt_dR_tracks_uncorr_ave);
  
  float jtentry[]={cent,pt1,eta1,phi1,pt2,eta2,phi2,pt3,eta3,phi3,dphi,ptratio,jetspz,jetspt,jetseta,jetsrap,jetsphi,jetspz_12,jetspt_12,jetseta_12,jetsrap_12,jetsphi_12,trkrap,mpt_tracks,mpt_tracks_ls,mpt_tracks_uncorr,mpt_tracks_ave,mpt_tracks_uncorr_ave,hfp,hfm,Npassingtrk,vz};
  
 
  nt_jet->Fill(jtentry);
 }
 

 
 outf->cd();
 if(doSaveTrackInfo){ 
  nt_track->Write();
 }
 nt_jet->Write();
 
  nt_mpt_tracks->Write();
  nt_mpt_tracks_dR->Write();
  nt_mpt_tracks_uncorr->Write();
  nt_mpt_tracks_uncorr_dR->Write();
  nt_mpt_tracks_dR_ave->Write();
  nt_mpt_tracks_uncorr_dR_ave->Write();
  // nt_track->Write();
 outf->Close();
 }
