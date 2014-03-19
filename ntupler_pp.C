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
#include "ntupler/ppTrack.C"

void ntupler_pp(double ptmin_trk=0.5,double ptmax_trk=300){ 
 TH1D::SetDefaultSumw2();

 TString algo="ak3Calo"; 
 //input file  
 cout<<"ptmin= "<<ptmin_trk<<" ptmax= "<<ptmax_trk<<endl;
 TString directory="root://eoscms//eos/cms/store/caf/user/yjlee//pp2013/promptReco/";
 TString infname="PP2013_HiForest_PromptReco_JsonPP_Jet80_PPReco_forestv82"; 
    
 ppTrack * ftrk = new ppTrack(Form("%s/%s.root",directory.Data(),infname.Data()));
 HiTree * fhi = new HiTree(Form("%s/%s.root",directory.Data(),infname.Data()));
 HltTree_pp * fhlt = new HltTree_pp(Form("%s/%s.root",directory.Data(),infname.Data()));
 skimTree_pp * fskim = new skimTree_pp(Form("%s/%s.root",directory.Data(),infname.Data()));
 t * fjet = new t(Form("%s/%s.root",directory.Data(),infname.Data()),algo);
 //getting unfolding profile 
 TFile * f_unfold=new TFile("binbybin_pt.root");
 TH1D *hweight = (TH1D*)f_unfold->Get("hweight"); 
 TH1D *h_res = (TH1D*)f_unfold->Get("h_res");
 TF1 * fgaus=new TF1("fgaus","gaus(0)",-20,20);
 fgaus->SetParameters(1,0,1);
  
 //pt bins for track efficiency correction
 int npt=4; 
 double ptmin[]={0.4, 1, 3,   8};
 double ptmax[]={  1, 3,  8,300};
 

 double frac[]={0,0.01,0.02,0.03,0.04,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.925,0.95,0.96,0.97,0.98,0.99,1.01};
 double dR_upperbound[]={0,0.01,0.02,0.03,0.04,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,1,1.1,1.34};
 // double dR_upperbound[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.01};

 //getting histograms for track efficiency correction 
 TFile *f_eff[npt];
 TProfile2D *p_eff_accept[npt]; 
 TProfile *p_eff_pt[npt]; 
 TProfile *p_eff_rmin[npt]; 
 for(int ipt=0; ipt<npt;ipt++){
   if(ipt<npt-1)f_eff[ipt]= new TFile(Form("/afs/cern.ch/user/d/dgulhan/workDir/trackEfficiency/final_hists_%s/eff_pt%d_%d_accept4pt4rmin3_%s_dogenjet0.root",algo.Data(),(int)ptmin[ipt],(int)ptmax[ipt],algo.Data()));
   else f_eff[ipt]= new TFile(Form("/afs/cern.ch/user/d/dgulhan/workDir/trackEfficiency/final_hists_%s/eff_pt%d_%d_accept3pt3rmin3_%s_dogenjet0.root",algo.Data(),(int)ptmin[ipt],(int)ptmax[ipt],algo.Data()));
   p_eff_pt[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_pt");
   p_eff_accept[ipt]=(TProfile2D*)f_eff[ipt]->Get("p_eff_acceptance");
   p_eff_rmin[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_rmin");
 }
 TFile *f_fake[npt];
 TProfile2D *p_fake_accept[npt]; 
 TProfile *p_fake_pt[npt];  
 TProfile *p_fake_rmin[npt]; 
 for(int ipt=0; ipt<npt;ipt++){
   f_fake[ipt]= new TFile(Form("/afs/cern.ch/user/d/dgulhan/workDir/trackFake_pp/final_hists_%s/fake_pt%d_%d_step_accept5pt5rmin4_%s_dogenjet0.root",algo.Data(),(int)ptmin[ipt],(int)ptmax[ipt],algo.Data()));
   p_fake_pt[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_pt");
   p_fake_accept[ipt]=(TProfile2D*)f_fake[ipt]->Get("p_fake_acceptance");
   p_fake_rmin[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_rmin");
 }
 //output file and tree
 // TFile *outf= new TFile(Form("track_ntuple_%s_jetscut_trackercut_jetframe_v2_pt%d_%d_nogenptcut.root",infname.Data(),(int)ptmin_trk,(int)ptmax_trk),"recreate");
 TFile *outf= new TFile(Form("ntuples_data_pp/hastracks_full_ntuple_HiForest_PP2013_HiForest_PromptReco_JsonPP_Jet80_PPReco_forestv82_pt%d_%d_%s.root",(int)ptmin_trk,(int)ptmax_trk,algo.Data()),"recreate");
 
 std::string trackVars="trackselect:eff:fake:trkfake:weight_unfold:pt:eta:phi:rmin:pt_boosted:p_boosted:eta_boosted:phi_boosted:pt1:eta1:phi1:pt1_boosted:eta1_boosted:phi1_boosted:pt2:eta2:phi2:pt2_boosted:eta2_boosted:phi2_boosted:pt3:eta3:phi3:pt3_boosted:eta3_boosted:phi3_boosted:dphi:ptratio:jetspz:jetspt:jetseta:jetsrap:jetspz_12:jetspt_12:jetsrap_12:alpha:alpha_boosted:mpt_track:mpt_boosted_track:hfp:hfm:vz";
 

 TNtuple *nt_track = new TNtuple("nt_track","",trackVars.data()); 
 
 std:string jetVars="pt1:eta1:phi1:pt1_boosted:eta1_boosted:phi1_boosted:pt2:eta2:phi2:pt2_boosted:eta2_boosted:phi2_boosted:pt3:eta3:phi3:pt3_boosted:eta3_boosted:phi3_boosted:dphi:ptratio:jetspz:jetspt:jetseta:jetsrap:jetsphi:jetspz_12:jetspt_12:jetseta_12:jetsrap_12:jetsphi_12:trkrap:mpt_tracks:mpt_boosted_tracks:mp_tracks:mp_boosted_tracks:mpt_tracks_uncorr:mpt_boosted_tracks_uncorr:mp_tracks_uncorr:mp_boosted_tracks_uncorr:hfp:hfm:vz";
 
 TNtuple *nt_jet = new TNtuple("nt_jet","",jetVars.data());

 std::string mptVars="alpha_0:alpha_1:alpha_2:alpha_3:alpha_4:alpha_5:alpha_6:alpha_7:alpha_8:alpha_9:alpha_10:alpha_11:alpha_12:alpha_13:alpha_14:alpha_15:alpha_16:alpha_17:alpha_18:alpha_19:alpha_20:alpha_21:alpha_22:alpha_23:alpha_24:alpha_25:alpha_26:alpha_27:alpha_28:alpha_29";
 

 TNtuple *nt_mpt_tracks=new TNtuple("nt_mpt_tracks","",mptVars.data());
 TNtuple *nt_mpt_boosted_tracks=new TNtuple("nt_mpt_boosted_tracks","",mptVars.data());
 TNtuple *nt_mp_tracks=new TNtuple("nt_mp_tracks","",mptVars.data());
 TNtuple *nt_mp_boosted_tracks=new TNtuple("nt_mp_boosted_tracks","",mptVars.data());

 TNtuple *nt_mpt_tracks_uncorr=new TNtuple("nt_mpt_tracks_uncorr","",mptVars.data());
 TNtuple *nt_mpt_boosted_tracks_uncorr=new TNtuple("nt_mpt_boosted_tracks_uncorr","",mptVars.data());
 TNtuple *nt_mp_tracks_uncorr=new TNtuple("nt_mp_tracks_uncorr","",mptVars.data());
 TNtuple *nt_mp_boosted_tracks_uncorr=new TNtuple("nt_mp_boosted_tracks_uncorr","",mptVars.data());
 
 std::string mptVars_dR="dR_0:dR_1:dR_2:dR_3:dR_4:dR_5:dR_6:dR_7:dR_8:dR_9:dR_10:dR_11:dR_12:dR_13:dR_14:dR_15:dR_16:dR_17:dR_18:dR_19:dR_20:dR_21:dR_22:dR_23:dR_24:dR_25:dR_26:dR_27:dR_28:dR_29";
 
 TNtuple *nt_mpt_tracks_dR=new TNtuple("nt_mpt_tracks_dR","",mptVars.data());
 TNtuple *nt_mpt_boosted_tracks_dR=new TNtuple("nt_mpt_boosted_tracks_dR","",mptVars.data());
 TNtuple *nt_mp_tracks_dR=new TNtuple("nt_mp_tracks_dR","",mptVars.data());
 TNtuple *nt_mp_boosted_tracks_dR=new TNtuple("nt_mp_boosted_tracks_dR","",mptVars_dR.data());
 
 TNtuple *nt_mpt_tracks_uncorr_dR=new TNtuple("nt_mpt_tracks_uncorr_dR","",mptVars.data());
 TNtuple *nt_mpt_boosted_tracks_uncorr_dR=new TNtuple("nt_mpt_boosted_tracks_uncorr_dR","",mptVars.data());
 TNtuple *nt_mp_tracks_uncorr_dR=new TNtuple("nt_mp_tracks_uncorr_dR","",mptVars.data());
 TNtuple *nt_mp_boosted_tracks_uncorr_dR=new TNtuple("nt_mp_boosted_tracks_uncorr_dR","",mptVars_dR.data());
 
 cout<<"2"<<endl;
 //loop over events
 int nentries = ftrk->GetEntriesFast(); 
 for(int jentry=0;jentry<nentries;jentry++){
 // for(int jentry=0;jentry<100;jentry++){
  if((jentry%1000)==0) std::cout<<jentry<<"/"<<nentries<<std::endl;

  ftrk->GetEntry(jentry);
  fhi->GetEntry(jentry);
  fskim->GetEntry(jentry);
  fjet->GetEntry(jentry);
  //fgen->GetEntry(jentry);
  if(!(fskim->pPAcollisionEventSelectionPA && fskim->pHBHENoiseFilter))continue;
  
  float hfp = fhi->hiHFplusEta4;
  float hfm = fhi->hiHFminusEta4;
  float vz = fhi->vz;
  float pt1=-99;
  float phi1=-99;
  float eta1=-99;
  float pt1_boosted=-99;
  float phi1_boosted=-99;
  float eta1_boosted=-99;
  float refpt1=-99;
  float refeta1=-99;
  float refphi1=-99;
  float neutralMax1=-99;
  float trackMax1=-99;
  float trackN1=-99;
  float pt2=-99;
  float phi2=-99;
  float eta2=-99;
  float pt2_boosted=-99;
  float phi2_boosted=-99;
  float eta2_boosted=-99;
  float refpt2=-99;
  float refphi2=-99;
  float refeta2=-99;
  float neutralMax2=-99;
  float trackMax2=-99;
  float trackN2=-99;
  float pt3=-99;
  float phi3=-99;
  float eta3=-99;
  float pt3_boosted=-99;
  float phi3_boosted=-99;
  float eta3_boosted=-99;
  float refpt3=-99;
  float refeta3=-99;
  float refphi3=-99;
  float neutralMax3=-99;
  float trackMax3=-99;
  float trackN3=-99; 
  float dphi=-99;
  float ptratio=-99;
  std::vector<std::pair<double, std::pair<double,double> >  >jets;
  
  int njet=0;
  TLorentzVector vjets;
  for(int ijet=0;ijet<fjet->nref;ijet++){
   if(fabs(fjet->jteta[ijet])>1.6 ||fjet->jtpt[ijet]<30) continue;
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
  TLorentzVector v1_boosted,v2_boosted,v3_boosted;

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
   v1_boosted.SetPtEtaPhiM(pt1,eta1,phi1,0);
   vls+=v1;
   if(njet>1){
    pt2=jets[njet-2].first;
    eta2=jets[njet-2].second.first;
    phi2=jets[njet-2].second.second;
    v2.SetPtEtaPhiM(pt2,eta2,phi2,0);
    v2_boosted.SetPtEtaPhiM(pt2,eta2,phi2,0);
    dphi=acos(cos(phi1-phi2));
    ptratio=pt2/pt1;
	  vls+=v2;

    if(njet>2){
     pt3=jets[njet-3].first;
     eta3=jets[njet-3].second.first;
     phi3=jets[njet-3].second.second;
     v3.SetPtEtaPhiM(pt2,eta2,phi2,0);
     v3_boosted.SetPtEtaPhiM(pt2,eta2,phi2,0);
    }
   }
  }
  
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
    v1_boosted.Boost(0,0,-tanh(jetsrap_12));
    pt1_boosted=v1_boosted.Pt();
    eta1_boosted=v1_boosted.Eta();
    phi1_boosted=v1_boosted.Phi();
  	v2_boosted.Boost(0,0,-tanh(jetsrap_12));
    pt2_boosted=v2_boosted.Pt();
    eta2_boosted=v2_boosted.Eta();
    phi2_boosted=v2_boosted.Phi();
	v3_boosted.Boost(0,0,-tanh(jetsrap_12));
    pt3_boosted=v3_boosted.Pt();
    eta3_boosted=v3_boosted.Eta();
    phi3_boosted=v3_boosted.Phi();
  }
  
  
  
  float mpt=0;
  float mp=0;
  float mpt_boosted=0;
  float mp_boosted=0;
  float mpt_tracks=0;
  float mp_tracks=0;
  float mpt_boosted_tracks=0;
  float mp_boosted_tracks=0;
  
  float mp_alpha_boosted_tracks[30];
  float mp_alpha_tracks[30];
  float mpt_alpha_boosted_tracks[30];
  float mpt_alpha_tracks[30];
  
  float mp_dR_boosted_tracks[30];
  float mp_dR_tracks[30];
  float mpt_dR_boosted_tracks[30];
  float mpt_dR_tracks[30];
  
  float mpt_tracks_uncorr=0;
  float mp_tracks_uncorr=0;
  float mpt_boosted_tracks_uncorr=0;
  float mp_boosted_tracks_uncorr=0;
  
  float mp_alpha_boosted_tracks_uncorr[30];
  float mp_alpha_tracks_uncorr[30];
  float mpt_alpha_boosted_tracks_uncorr[30];
  float mpt_alpha_tracks_uncorr[30];
  
  float mp_dR_boosted_tracks_uncorr[30];
  float mp_dR_tracks_uncorr[30];
  float mpt_dR_boosted_tracks_uncorr[30];
  float mpt_dR_tracks_uncorr[30];
  
  for(int i=0;i<30;i++){
   mp_alpha_tracks[i]=0;  
   mp_alpha_boosted_tracks[i]=0;  
   mpt_alpha_tracks[i]=0;  
   mpt_alpha_boosted_tracks[i]=0;  
   mp_dR_tracks[i]=0;  
   mp_dR_boosted_tracks[i]=0;  
   mpt_dR_tracks[i]=0;  
   mpt_dR_boosted_tracks[i]=0;  
   
   mp_alpha_tracks_uncorr[i]=0;  
   mp_alpha_boosted_tracks_uncorr[i]=0;  
   mpt_alpha_tracks_uncorr[i]=0;  
   mpt_alpha_boosted_tracks_uncorr[i]=0;  
   mp_dR_tracks_uncorr[i]=0;  
   mp_dR_boosted_tracks_uncorr[i]=0;  
   mpt_dR_tracks_uncorr[i]=0;  
   mpt_dR_boosted_tracks_uncorr[i]=0;  
  }
  
  TLorentzVector vtracks;
  
  for(int itrk=0;itrk<ftrk->nTrk;itrk++){
   float eta=ftrk->trkEta[itrk];
   float pt=ftrk->trkPt[itrk];
   float phi=ftrk->trkPhi[itrk];
   
   float trkfake=ftrk->trkFake[itrk];
   float eff_pt,eff_accept,eff_rmin;
   eff_pt=eff_accept=eff_rmin=1;
   float fake_pt,fake_accept,fake_rmin;
   fake_pt=fake_accept=fake_rmin=0;
   if(pt<ptmin_trk || pt>ptmax_trk) continue; 
   if(abs(eta)>2.4) continue;
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
    if(pt>=ptmin[ipt] && pt<ptmax[ipt]){
      eff_pt=p_eff_pt[ipt]->GetBinContent(p_eff_pt[ipt]->FindBin(pt));
      eff_accept=p_eff_accept[ipt]->GetBinContent(p_eff_accept[ipt]->GetXaxis()->FindBin(phi),p_eff_accept[ipt]->GetYaxis()->FindBin(eta));
      if(rmin<5)eff_rmin=p_eff_rmin[ipt]->GetBinContent(p_eff_rmin[ipt]->FindBin(rmin));
     }     
   } 
   
   for(int ipt=0;ipt<npt;ipt++){
    if(pt>=ptmin[ipt] && pt<ptmax[ipt]){
      fake_pt=p_fake_pt[ipt]->GetBinContent(p_fake_pt[ipt]->FindBin(pt));
      fake_accept=p_fake_accept[ipt]->GetBinContent(p_fake_accept[ipt]->GetXaxis()->FindBin(phi),p_fake_accept[ipt]->GetYaxis()->FindBin(eta));
      if(rmin<5) fake_rmin=p_fake_rmin[ipt]->GetBinContent(p_fake_rmin[ipt]->FindBin(rmin));
     }     
   }
  
   float eff=1;
   eff=eff_accept*eff_pt*eff_rmin;
   if(eff==0){
    if(pt>100)eff=0.8;
	else eff=1;
   }
   float fake=0;
   if(pt<100)fake=fake_accept+fake_pt+fake_rmin;
   
   TLorentzVector vp;
   vp.SetPtEtaPhiM(pt,eta,phi,0.13957018);
   vtracks+=vp;
   
   float p_boosted=-99;
   float pt_boosted=-99;
   float eta_boosted=-99;
   float phi_boosted=-99;
   float alpha=-99;
   float alpha_boosted=-99;
   float mpt_boosted_track=-99;
    
   float mpt_track=pt*cos(phi-phi1);
   alpha=vp.Angle(v1.Vect());
   if(trackselect){
    // cout<<"eff="<<eff<<" fake="<<fake<<" (1-fake)/eff="<<(1-fake)/eff<<" pt="<<pt<<" eta="<<eta<<" phi="<<phi<<endl;
    mpt_tracks+=((1-fake)/eff)*pt*cos(phi-phi1);
    mpt_tracks_uncorr+=pt*cos(phi-phi1);
    mp_tracks+=((1-fake)/eff)*vp.P()*cos(alpha);
    mp_tracks_uncorr+=vp.P()*cos(alpha);
   }
   if(njet>1){
    TLorentzVector vp_boosted;
    vp_boosted.SetPtEtaPhiM(pt,eta,phi,0.13957018);
    vp_boosted.Boost(0,0,-tanh(jetsrap_12)); 
	
   	pt_boosted=vp_boosted.Pt();
  	p_boosted=vp_boosted.P();
  	eta_boosted=vp_boosted.Eta();
  	phi_boosted=vp_boosted.Phi();
	
	  mpt_boosted_track=pt_boosted*cos(phi_boosted-phi1_boosted);
    mpt_boosted_tracks+=pt_boosted*cos(phi_boosted-phi1_boosted);
    alpha_boosted=vp_boosted.Angle(v1_boosted.Vect());
    mp_boosted_tracks+=vp_boosted.P()*cos(alpha_boosted);
	  
	  
    float entry[]={trackselect,eff,fake,trkfake,weight_unfold,pt,eta,phi,rmin,pt_boosted,p_boosted,eta_boosted,phi_boosted,pt1,eta1,phi1,pt1_boosted,eta1_boosted,phi1_boosted,pt2,eta2,phi2,pt2_boosted,eta2_boosted,phi2_boosted,pt3,eta3,phi3,pt3_boosted,eta3_boosted,phi3_boosted,dphi,ptratio,jetspz,jetspt,jetseta,jetsrap,jetspz_12,jetspt_12,jetsrap_12,alpha,alpha_boosted,mpt_track,mpt_boosted_track,hfp,hfm,vz};
    nt_track->Fill(entry);

    if(!(trackselect)) continue;

    double mpt_alpha_t[30];
    double mpt_alpha_boosted_t[30];
    double mp_alpha_t[30];
    double mp_alpha_boosted_t[30];
    
    double mpt_dR_t[30];
    double mpt_dR_boosted_t[30];
    double mp_dR_t[30];
    double mp_dR_boosted_t[30];
    
   for(int i=0;i<30;i++){
      mpt_alpha_t[i]=0;
      mpt_alpha_boosted_t[i]=0;
	  if(alpha <(frac[i+1]*TMath::Pi()) && alpha >=(frac[i]*TMath::Pi())) mpt_alpha_t[i] =pt*cos(phi -phi1);
      mpt_alpha_tracks[i]+=((1-fake)/eff)*mpt_alpha_t[i];
      mpt_alpha_tracks_uncorr[i]+=mpt_alpha_t[i];
	  if(alpha_boosted <(frac[i+1]*TMath::Pi()) && alpha_boosted >=(frac[i-1]*TMath::Pi())) mpt_alpha_boosted_t[i] =pt *cos(phi -phi1);
      mpt_alpha_boosted_tracks[i]+=((1-fake)/eff)*mpt_alpha_boosted_t[i];
      mpt_alpha_boosted_tracks_uncorr[i]+=mpt_alpha_boosted_t[i];
    
      mp_alpha_t[i]=0;
      mp_alpha_boosted_t[i]=0;
	  if(alpha <(frac[i+1]*TMath::Pi()) && alpha >=(frac[i]*TMath::Pi())) mp_alpha_t[i] =vp.P()*cos(alpha);
      mp_alpha_tracks[i]+=((1-fake)/eff)*mp_alpha_t[i];
      mp_alpha_tracks_uncorr[i]+=mp_alpha_t[i];
	  if(alpha_boosted <(frac[i+1]*TMath::Pi()) && alpha_boosted >=(frac[i]*TMath::Pi())) mp_alpha_boosted_t[i] =vp_boosted.P()*cos(alpha_boosted);
      mp_alpha_boosted_tracks[i]+=((1-fake)/eff)*mp_alpha_boosted_t[i];
      mp_alpha_boosted_tracks_uncorr[i]+=mp_alpha_boosted_t[i];
	}
	
      for(int i=0;i<30;i++){
        mpt_dR_t[i]=0;
        mpt_dR_boosted_t[i]=0;
	    double dR1=sqrt(pow(acos(cos(phi1-phi)),2)+pow(eta1-eta,2));
	    double dR2=sqrt(pow(acos(cos(phi2-phi)),2)+pow(eta2-eta,2));
	    double dR1_boosted=sqrt(pow(acos(cos(phi1_boosted-phi_boosted)),2)+pow(eta1_boosted-eta_boosted,2));
	    double dR2_boosted=sqrt(pow(acos(cos(phi2_boosted-phi_boosted)),2)+pow(eta2_boosted-eta_boosted,2));
	    if(((dR1 < dR_upperbound[i+1]) && dR1 >=dR_upperbound[i]) || ((dR2 < dR_upperbound[i+1])&& dR2 >=dR_upperbound[i])) mpt_dR_t[i] =pt*cos(phi -phi1);
        mpt_dR_tracks[i]+=((1-fake)/eff)*mpt_dR_t[i];
        mpt_dR_tracks_uncorr[i]+=mpt_dR_t[i];
	    if(((dR1 < dR_upperbound[i+1]) && dR1 >=dR_upperbound[i]) || ((dR2 < dR_upperbound[i+1]) && dR2 >=dR_upperbound[i])) mpt_dR_boosted_t[i] =pt *cos(phi -phi1);
        mpt_dR_boosted_tracks[i]+=((1-fake)/eff)*mpt_dR_boosted_t[i];
        mpt_dR_boosted_tracks_uncorr[i]+=mpt_dR_boosted_t[i];
        mp_dR_t[i]=0;
        mp_dR_boosted_t[i]=0;
	    if(((dR1_boosted < dR_upperbound[i+1]) && dR1_boosted >=dR_upperbound[i]) || ((dR2_boosted < dR_upperbound[i+1])&& dR2_boosted >=dR_upperbound[i])) mp_dR_t[i] =vp.P()*cos(alpha);
        mp_dR_tracks[i]+=((1-fake)/eff)*mp_dR_t[i];
        mp_dR_tracks_uncorr[i]+=mp_dR_t[i];
	    if(((dR1_boosted < dR_upperbound[i+1]) && dR1_boosted >=dR_upperbound[i]) || ((dR2_boosted < dR_upperbound[i+1]) && dR2_boosted>=dR_upperbound[i])) mp_dR_boosted_t[i] =vp_boosted.P()*cos(alpha_boosted);
        mp_dR_boosted_tracks[i]+=((1-fake)/eff)*mp_dR_boosted_t[i];
        mp_dR_boosted_tracks_uncorr[i]+=mp_dR_boosted_t[i];
	   }
   }
   //fill in the output tree

  }
 
   
  double trkrap=vtracks.Rapidity();
  
  nt_mpt_tracks->Fill(mpt_alpha_tracks);
  nt_mpt_boosted_tracks->Fill(mpt_alpha_boosted_tracks);
  nt_mp_tracks->Fill(mp_alpha_tracks);
  nt_mp_boosted_tracks->Fill(mp_alpha_boosted_tracks);
 
  nt_mpt_tracks_dR->Fill(mpt_dR_tracks);
  nt_mpt_boosted_tracks_dR->Fill(mpt_dR_boosted_tracks);
  nt_mp_tracks_dR->Fill(mp_dR_tracks);
  nt_mp_boosted_tracks_dR->Fill(mp_dR_boosted_tracks);
 
  nt_mpt_tracks_uncorr->Fill(mpt_alpha_tracks_uncorr);
  nt_mpt_boosted_tracks_uncorr->Fill(mpt_alpha_boosted_tracks_uncorr);
  nt_mp_tracks_uncorr->Fill(mp_alpha_tracks_uncorr);
  nt_mp_boosted_tracks_uncorr->Fill(mp_alpha_boosted_tracks_uncorr);
 
  nt_mpt_tracks_uncorr_dR->Fill(mpt_dR_tracks_uncorr);
  nt_mpt_boosted_tracks_uncorr_dR->Fill(mpt_dR_boosted_tracks_uncorr);
  nt_mp_tracks_uncorr_dR->Fill(mp_dR_tracks_uncorr);
  nt_mp_boosted_tracks_uncorr_dR->Fill(mp_dR_boosted_tracks_uncorr);
 
  float jtentry[]={pt1,eta1,phi1,pt1_boosted,eta1_boosted,phi1_boosted,pt2,eta2,phi2,pt2_boosted,eta2_boosted,phi2_boosted,pt3,eta3,phi3,pt3_boosted,eta3_boosted,phi3_boosted,dphi,ptratio,jetspz,jetspt,jetseta,jetsrap,jetsphi,jetspz_12,jetspt_12,jetseta_12,jetsrap_12,jetsphi_12,trkrap,mpt_tracks,mpt_boosted_tracks,mp_tracks,mp_boosted_tracks,mpt_tracks_uncorr,mpt_boosted_tracks_uncorr,mp_tracks_uncorr,mp_boosted_tracks_uncorr,hfp,hfm,vz};
  
  nt_jet->Fill(jtentry);
 }
 

 
 outf->cd();
 nt_track->Write();
 nt_jet->Write();

 nt_mpt_tracks->Write();
 nt_mpt_boosted_tracks->Write();
 nt_mp_tracks->Write();
 nt_mp_boosted_tracks->Write();
 
 nt_mpt_tracks_dR->Write();
 nt_mpt_boosted_tracks_dR->Write();
 nt_mp_tracks_dR->Write();
 nt_mp_boosted_tracks_dR->Write();
 
 nt_mpt_tracks_uncorr->Write();
 nt_mpt_boosted_tracks_uncorr->Write();
 nt_mp_tracks_uncorr->Write();
 nt_mp_boosted_tracks_uncorr->Write();
 
 nt_mpt_tracks_uncorr_dR->Write();
 nt_mpt_boosted_tracks_uncorr_dR->Write();
 nt_mp_tracks_uncorr_dR->Write();
 nt_mp_boosted_tracks_uncorr_dR->Write();
 
 outf->Close();
 }
