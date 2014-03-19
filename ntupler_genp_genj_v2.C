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
#include "../../trackEfficiency/v2/ntupler/trackTree.C"

void ntupler_genp_genj_v2(){ 
 TH1D::SetDefaultSumw2();
 double ptmin_trk=0.5; 
 double ptmax_trk=300;
 //input file
 cout<<"ptmin= "<<ptmin_trk<<" ptmax= "<<ptmax_trk<<endl;
 TString directory="/mnt/hadoop/cms/store/user/dgulhan/HIMC/Jet80/Track8_Jet4_STARTHI53_LV1/merged2/";
 TString infname="HiForest_Pythia_Hydjet_Jet80_Track8_Jet4_STARTHI53_LV1_merged_forest_0"; 
  
 trackTree * ftrk = new trackTree(Form("%s/%s.root",directory.Data(),infname.Data()));
 HiTree * fhi = new HiTree(Form("%s/%s.root",directory.Data(),infname.Data()));
 t * fjet = new t(Form("%s/%s.root",directory.Data(),infname.Data()),"akPu3PF");
 hi * fgen = new hi(Form("%s/%s.root",directory.Data(),infname.Data()));

 //pt bins for track efficiency correction
 int npt=5;
 double ptmin[]={0.4,1,3,8,100};
 double ptmax[]={1,3,8,100,300};
 double frac[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.01};
 
 //getting unfolding profile 
 TF1 * fit = new TF1("fit","[0]*x+[1]*x*x+[2]",0,50);
 fit->SetParameters(9.93930e-01,1.99078e-04, 4.18899e-03); 

 // getting histograms for track efficiency  

 TFile *f_eff[npt];
 TProfile *p_eff_cent[npt]; 
 TProfile2D *p_eff_accept[npt]; 
 TProfile *p_eff_pt[npt]; 
 TProfile *p_eff_rmin[npt]; 
 for(int ipt=0; ipt<npt;ipt++){
   // f_eff[ipt]= new TFile(Form("../../trackEfficiency/v2/final_hists/eff_pt%d_%d_step_cent3accept3pt3rmin2_akPu3PF_dogenjet1.root",(int)ptmin[ipt],(int)ptmax[ipt]));
   f_eff[ipt]= new TFile(Form("../../trackEfficiency/v2/final_hists/eff_pt%d_%d_step_cent3accept3pt3rmin2.root",(int)ptmin[ipt],(int)ptmax[ipt]));
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
   f_fake[ipt]= new TFile(Form("../../trackFake/final_hists/fake_pt%d_%d_step_cent2accept2pt2rmin1_akPu3PF_dogenjet1.root",(int)ptmin[ipt],(int)ptmax[ipt]));
   p_fake_cent[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_cent");
   p_fake_pt[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_pt");
   p_fake_accept[ipt]=(TProfile2D*)f_fake[ipt]->Get("p_fake_acceptance");
   p_fake_rmin[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_rmin");
 }
 
 //output file and tree
 TFile *outf= new TFile(Form("track_ntuple_%s_jetscut_trackercut_jetframe_v2_pt%d_%d.root",infname.Data(),(int)ptmin_trk,(int)ptmax_trk),"recreate");
 std::string partVars="pt:eta:phi:rmin:cent:pt_boosted:p_boosted:eta_boosted:phi_boosted:cent:pt1:eta1:phi1:pt1_boosted:eta1_boosted:phi1_boosted:pt2:eta2:phi2:pt2_boosted:eta2_boosted:phi2_boosted:pt3:eta3:phi3:pt3_boosted:eta3_boosted:phi3_boosted:dphi:ptratio:jetspz:jetspt:jetseta:jetsrap:jetspz_12:jetspt_12:jetsrap_12:alpha:alpha_boosted:mpt_track:mpt_boosted_track:mpt_track_s:mpt_boosted_track_s:mpt_track_b:mpt_boosted_track_b";
 std::string trackVars="trackselect:eff:fake:pt:pt_unfold:eta:phi:rmin:cent:pt_boosted:p_boosted:eta_boosted:phi_boosted:cent:pt1:eta1:phi1:pt1_boosted:eta1_boosted:phi1_boosted:pt2:eta2:phi2:pt2_boosted:eta2_boosted:phi2_boosted:pt3:eta3:phi3:pt3_boosted:eta3_boosted:phi3_boosted:dphi:ptratio:jetspz:jetspt:jetseta:jetsrap:jetspz_12:jetspt_12:jetsrap_12:alpha:alpha_boosted:mpt_track:mpt_boosted_track:mpt_track_s:mpt_boosted_track_s:mpt_track_b:mpt_boosted_track_b";
 
 TNtuple *nt_track = new TNtuple("nt_track","",trackVars.data());
 TNtuple *nt_particle = new TNtuple("nt_particle","",partVars.data());
 
 std:string jetVars="cent:pt1:eta1:phi1:pt1_boosted:eta1_boosted:phi1_boosted:pt2:eta2:phi2:pt2_boosted:eta2_boosted:phi2_boosted:pt3:eta3:phi3:pt3_boosted:eta3_boosted:phi3_boosted:dphi:ptratio:jetspz:jetspt:jetseta:jetsrap:jetsphi:jetspz_12:jetspt_12:jetseta_12:jetsrap_12:jetsphi_12:jetsrap_reco:partrap:trkrap:signalrap:backgroundrap:mpt:mpt_boosted:mp:mp_boosted";
 TNtuple *nt_jet = new TNtuple("nt_jet","",jetVars.data());

 std::string mptVars="alpha_0:alpha_1:alpha_2:alpha_3:alpha_4:alpha_5:alpha_6:alpha_7:alpha_8:alpha_9:alpha_10:alpha_11:alpha_12:alpha_13:alpha_14:alpha_15:alpha_16:alpha_17:alpha_18:alpha_19";
 
 TNtuple *nt_mpt=new TNtuple("nt_mpt","",mptVars.data());
 TNtuple *nt_mpt_boosted=new TNtuple("nt_mpt_boosted","",mptVars.data());
 TNtuple *nt_mp=new TNtuple("nt_mp","",mptVars.data());
 TNtuple *nt_mp_boosted=new TNtuple("nt_mp_boosted","",mptVars.data());
 
 TNtuple *nt_mpt_b=new TNtuple("nt_mpt_b","",mptVars.data());
 TNtuple *nt_mpt_boosted_b=new TNtuple("nt_mpt_boosted_b","",mptVars.data());
 TNtuple *nt_mp_b=new TNtuple("nt_mp_b","",mptVars.data());
 TNtuple *nt_mp_boosted_b=new TNtuple("nt_mp_boosted_b","",mptVars.data());
 
 TNtuple *nt_mpt_s=new TNtuple("nt_mpt_s","",mptVars.data());
 TNtuple *nt_mpt_boosted_s=new TNtuple("nt_mpt_boosted_s","",mptVars.data());
 TNtuple *nt_mp_s=new TNtuple("nt_mp_s","",mptVars.data());
 TNtuple *nt_mp_boosted_s=new TNtuple("nt_mp_boosted_s","",mptVars.data());
 
 TNtuple *nt_mpt_tracks=new TNtuple("nt_mpt_tracks","",mptVars.data());
 TNtuple *nt_mpt_boosted_tracks=new TNtuple("nt_mpt_boosted_tracks","",mptVars.data());
 TNtuple *nt_mp_tracks=new TNtuple("nt_mp_tracks","",mptVars.data());
 TNtuple *nt_mp_boosted_tracks=new TNtuple("nt_mp_boosted_tracks","",mptVars.data());
 
 cout<<"2"<<endl;
 //loop over events
 int nentries = ftrk->GetEntriesFast();
 for(int jentry=0;jentry<nentries;jentry++){
  if((jentry%1000)==0) std::cout<<jentry<<"/"<<nentries<<std::endl;

  ftrk->GetEntry(jentry);
  fhi->GetEntry(jentry);
  fjet->GetEntry(jentry);
  fgen->GetEntry(jentry);

  float cent=fhi->hiBin;
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
  for(int ijet=0;ijet<fjet->ngen;ijet++){
   if(fabs(fjet->geneta[ijet])>2 ||fjet->genpt[ijet]<30) continue;
   jets.push_back(std::make_pair(fjet->genpt[ijet],std::make_pair(fjet->geneta[ijet],fjet->genphi[ijet])));
  
   TLorentzVector vjet;
   vjet.SetPtEtaPhiM(fjet->genpt[ijet],fjet->geneta[ijet],fjet->genphi[ijet],0);
   vjets+=vjet;
   njet++;
  } 
  
  TLorentzVector vjets_reco;
  int njet_reco;
  for(int ijet=0;ijet<fjet->nref;ijet++){
   if(fabs(fjet->jteta[ijet])>2 ||fjet->jtpt[ijet]<30) continue;
    TLorentzVector vjet;
    vjet.SetPtEtaPhiM(fjet->jtpt[ijet],fjet->jteta[ijet],fjet->jtphi[ijet],0);
    vjets_reco+=vjet;
    njet_reco++;
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
   pt1=       jets[njet-1].first;
   eta1=      jets[njet-1].second.first;
   phi1=      jets[njet-1].second.second;
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
    v1_boosted.Boost(0,0,-tanh(jetsrap));
    pt1_boosted=v1_boosted.Pt();
    eta1_boosted=v1_boosted.Eta();
    phi1_boosted=v1_boosted.Phi();
  	v2_boosted.Boost(0,0,-tanh(jetsrap));
    pt2_boosted=v2_boosted.Pt();
    eta2_boosted=v2_boosted.Eta();
    phi2_boosted=v2_boosted.Phi();
	v3_boosted.Boost(0,0,-tanh(jetsrap));
    pt3_boosted=v3_boosted.Pt();
    eta3_boosted=v3_boosted.Eta();
    phi3_boosted=v3_boosted.Phi();
  }
  
  
  double jetspt_reco=-99;
  double jetspx_reco=-99;
  double jetspy_reco=-99;
  double jetspz_reco=-99;
  double jetseta_reco=-99;
  double jetsphi_reco=-99;
  double jetsrap_reco=-99;
  
    
  if(njet_reco>1){ 
    jetspt_reco=vjets_reco.Pt();
    jetspx_reco=vjets_reco.Px();
    jetspy_reco=vjets_reco.Py();
    jetspz_reco=vjets_reco.Pz();
    jetseta_reco=vjets_reco.Eta();
    jetsphi_reco=vjets_reco.Phi();
    jetsrap_reco=vjets_reco.Rapidity();
  }
  
  float mpt=0;
  float mp=0;
  float mpt_boosted=0;
  float mp_boosted=0;
  float mpt_tracks=0;
  float mp_tracks=0;
  float mpt_boosted_tracks=0;
  float mp_boosted_tracks=0;
  float mpt_s=0;
  float mp_s=0;
  float mpt_boosted_s=0;
  float mp_boosted_s=0;
  float mpt_b=0;
  float mp_b=0;
  float mpt_boosted_b=0;
  float mp_boosted_b=0;
  float mp_alpha_boosted_s[20];
  float mp_alpha_boosted_b[20];
  float mp_alpha_boosted[20];
  float mp_alpha_s[20];
  float mp_alpha_b[20];
  float mp_alpha[20];
  float mpt_alpha_boosted_s[20];
  float mpt_alpha_boosted_b[20];
  float mpt_alpha_boosted[20];
  float mpt_alpha_s[20];
  float mpt_alpha_b[20];
  float mpt_alpha[20];
  
  float mp_alpha_boosted_tracks[20];
  float mp_alpha_tracks[20];
  float mpt_alpha_boosted_tracks[20];
  float mpt_alpha_tracks[20];
  
  for(int i=0;i<20;i++){
   mp_alpha[i]=0;  
   mp_alpha_s[i]=0;  
   mp_alpha_b[i]=0;  
   mp_alpha_boosted[i]=0;  
   mp_alpha_boosted_s[i]=0;  
   mp_alpha_boosted_b[i]=0;  
   mpt_alpha[i]=0;  
   mpt_alpha_s[i]=0;  
   mpt_alpha_b[i]=0;  
   mpt_alpha_boosted[i]=0;  
   mpt_alpha_boosted_s[i]=0;  
   mpt_alpha_boosted_b[i]=0;  
   mp_alpha_tracks[i]=0;  
   mp_alpha_boosted_tracks[i]=0;  
   mpt_alpha_tracks[i]=0;  
   mpt_alpha_boosted_tracks[i]=0;  
  }
  
  TLorentzVector vparticles, vsignal, vbackground,vtracks;
  
  for(int itrk=0;itrk<ftrk->nTrk;itrk++){
   float eta=ftrk->trkEta[itrk];
   float pt=ftrk->trkPt[itrk];
   float phi=ftrk->trkPhi[itrk];
   float pt_unfold=pt;
   if(pt<50) pt_unfold=fit->Eval(pt);
   
   float eff_pt,eff_cent,eff_accept,eff_rmin;
   eff_pt=eff_cent=eff_accept=eff_rmin=1;
   float fake_pt,fake_cent,fake_accept,fake_rmin;
   fake_pt=fake_cent=fake_accept=fake_rmin=0;
   
   if(pt<ptmin_trk || pt>ptmax_trk) continue; 
   float trackselect=(ftrk->highPurity[itrk] && (ftrk->trkDxy1[itrk]/ftrk->trkDxyError1[itrk])<3.0 && (ftrk->trkDz1[itrk]/ftrk->trkDzError1[itrk])<3 && (ftrk->trkPtError[itrk]/ftrk->trkPt[itrk])<0.1);
   
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
      eff_cent=p_eff_cent[ipt]->GetBinContent(p_eff_cent[ipt]->FindBin(cent));
      eff_accept=p_eff_accept[ipt]->GetBinContent(p_eff_accept[ipt]->GetXaxis()->FindBin(phi),p_eff_accept[ipt]->GetYaxis()->FindBin(eta));
      eff_rmin=p_eff_rmin[ipt]->GetBinContent(p_eff_rmin[ipt]->FindBin(rmin));
      fake_pt=p_fake_pt[ipt]->GetBinContent(p_fake_pt[ipt]->FindBin(pt));
      fake_cent=p_fake_cent[ipt]->GetBinContent(p_fake_cent[ipt]->FindBin(cent));
      fake_accept=p_fake_accept[ipt]->GetBinContent(p_fake_accept[ipt]->GetXaxis()->FindBin(phi),p_fake_accept[ipt]->GetYaxis()->FindBin(eta));
      fake_rmin=p_fake_rmin[ipt]->GetBinContent(p_fake_rmin[ipt]->FindBin(rmin));//fakeiciency for rmin>3 is 1.
      // if(rmin<3)eff_rmin=p_eff_rmin[ipt]->GetBinContent(p_eff_rmin[ipt]->FindBin(rmin));//efficiency for rmin>3 is 1.
     }     
   }
   
   float eff=1;
   if(pt<100)eff=eff_accept*eff_cent*eff_pt*eff_rmin;
   float fake=0;
   if(pt<100)fake=fake_accept+fake_cent+fake_pt+fake_rmin;
   
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
   mpt_tracks+=(1/eff)*pt*cos(phi-phi1);
   alpha=vp.Angle(v1.Vect());
   mp_tracks+=(1/eff)*vp.P()*cos(alpha);

   if(njet>1){
    TLorentzVector vp_boosted;
    vp_boosted.SetPtEtaPhiM(pt,eta,phi,0.13957018);
    vp_boosted.Boost(0,0,-tanh(jetsrap)); 
	
   	pt_boosted=vp_boosted.Pt();
  	p_boosted=vp_boosted.P();
  	eta_boosted=vp_boosted.Eta();
  	phi_boosted=vp_boosted.Phi();
	
    if(abs(phi_boosted)>TMath::Pi())continue;
    if(abs(eta_boosted)>90)continue;
    
	mpt_boosted_track=pt_boosted*cos(phi_boosted-phi1_boosted);
    mpt_boosted_tracks+=pt_boosted*cos(phi_boosted-phi1_boosted);
    alpha_boosted=vp_boosted.Angle(v1_boosted.Vect());
    mp_boosted_tracks+=vp_boosted.P()*cos(alpha_boosted);
	  
    double mpt_alpha_t[20];
    double mpt_alpha_boosted_t[20];
    double mp_alpha_t[20];
    double mp_alpha_boosted_t[20];
    
   for(int i=0;i<20;i++){
      mpt_alpha_t[i]=0;
      mpt_alpha_boosted_t[i]=0;
	  if(alpha <(frac[i+1]*TMath::Pi()) && alpha >=(frac[i]*TMath::Pi())) mpt_alpha_t[i] =pt*cos(phi -phi1);
      mpt_alpha_tracks[i]+=(1/eff)*mpt_alpha_t[i];
	  if(alpha_boosted <(frac[i+1]*TMath::Pi()) && alpha_boosted >=(frac[i-1]*TMath::Pi())) mpt_alpha_boosted_t[i] =pt *cos(phi -phi1);
      mpt_alpha_boosted_tracks[i]+=(1/eff)*mpt_alpha_boosted_t[i];
    
       mp_alpha_t[i]=0;
       mp_alpha_boosted_t[i]=0;
	   if(alpha <(frac[i+1]*TMath::Pi()) && alpha >=(frac[i]*TMath::Pi())) mp_alpha_t[i] =vp.P()*cos(alpha);
       mp_alpha_tracks[i]+=(1/eff)*mp_alpha_t[i];
	   if(alpha_boosted <(frac[i+1]*TMath::Pi()) && alpha_boosted >=(frac[i]*TMath::Pi())) mp_alpha_boosted_t[i] =vp_boosted.P()*cos(alpha_boosted);
       mp_alpha_boosted_tracks[i]+=(1/eff)*mp_alpha_boosted_t[i];
	  }
   }
   //fill in the output tree
 
   float entry[]={trackselect,eff,fake,pt,pt_unfold,eta,phi,rmin,cent,pt_boosted,p_boosted,eta_boosted,phi_boosted,cent,pt1,eta1,phi1,pt1_boosted,eta1_boosted,phi1_boosted,pt2,eta2,phi2,pt2_boosted,eta2_boosted,phi2_boosted,pt3,eta3,phi3,pt3_boosted,eta3_boosted,phi3_boosted,dphi,ptratio,jetspz,jetspt,jetseta,jetsrap,jetspz_12,jetspt_12,jetsrap_12,alpha,alpha_boosted,mpt_track,mpt_boosted_track};
   nt_track->Fill(entry);
  }
  
  for(int itrk=0;itrk<fgen->mult;itrk++){
   float eta=fgen->eta[itrk];
   if(fabs(eta)>2.4) continue; //acceptance of the tracker
   float pt=fgen->pt[itrk];
   if(pt<ptmin_trk || pt>ptmax_trk) continue; 
   // if(pt>1) continue;
   float phi=fgen->phi[itrk];
   
   TLorentzVector vp;
   vp.SetPtEtaPhiM(pt,eta,phi,0.13957018);
   vparticles+=vp;
   float signal=0;
   if(fgen->sube[itrk]==0)signal=1;
   if(signal==1)vsignal+=vp;
   if(signal==0)vbackground+=vp;
   
   float p_boosted=-99;
   float pt_boosted=-99;
   float eta_boosted=-99;
   float phi_boosted=-99;
   float alpha=-99;
   float alpha_boosted=-99;
   float mpt_boosted_track=-99;
   float mpt_boosted_track_s=-99;
   float mpt_boosted_track_b=-99;
   float  mpt_track_s=-99;
   float  mpt_track_b=-99;
    
	 float mpt_track=pt*cos(phi-phi1);
	 if(signal==1)mpt_track_s=pt*cos(phi-phi1);
   else mpt_track_b=pt*cos(phi-phi1);
   mpt+=mpt_track;
	 if(signal==1)mpt_s+=pt*cos(phi-phi1);
   else mpt_b+=pt*cos(phi-phi1);
	 alpha=vp.Angle(v1.Vect());
   mp+=vp.P()*cos(alpha);
	 if(signal==1)mp_s+=vp.P()*cos(alpha);
	 else mp_b+=vp.P()*cos(alpha);

   if(njet>1){
    TLorentzVector vp_boosted;
    vp_boosted.SetPtEtaPhiM(pt,eta,phi,0.13957018);
    vp_boosted.Boost(0,0,-tanh(jetsrap));
	
   	pt_boosted=vp_boosted.Pt();
  	p_boosted=vp_boosted.P();
  	eta_boosted=vp_boosted.Eta();
  	phi_boosted=vp_boosted.Phi();
	
    if(abs(phi_boosted)>TMath::Pi())continue;
    if(abs(eta_boosted)>90)continue;
    
	  mpt_boosted_track=pt_boosted*cos(phi_boosted-phi1_boosted);
    if(signal==1) mpt_boosted_track_s=pt_boosted*cos(phi_boosted-phi1_boosted);
    else mpt_boosted_track_b=pt_boosted*cos(phi_boosted-phi1_boosted);
    mpt_boosted+=pt_boosted*cos(phi_boosted-phi1_boosted);
    alpha_boosted=vp_boosted.Angle(v1_boosted.Vect());
    mp_boosted+=vp_boosted.P()*cos(alpha_boosted);
    if(signal==1) mp_boosted_s+=vp_boosted.P()*cos(alpha_boosted);
    else mp_boosted_b+=pt_boosted*cos(phi_boosted-phi1_boosted);
    if(signal==1) mpt_boosted_s+=pt_boosted*cos(phi_boosted-phi1_boosted);
    else mpt_boosted_b+=pt_boosted*cos(phi_boosted-phi1_boosted);
	
    double mpt_alpha_t[20];
    double mpt_alpha_boosted_t[20];
    double mp_alpha_t[20];
    double mp_alpha_boosted_t[20];
    
   for(int i=0;i<20;i++){
     mpt_alpha_t[i]=0;
     mpt_alpha_boosted_t[i]=0;
	 if(alpha <(frac[i+1]*TMath::Pi()) && alpha >=(frac[i]*TMath::Pi())) mpt_alpha_t[i] =pt*cos(phi -phi1);
     mpt_alpha[i]+=mpt_alpha_t[i];
     if(signal==1) mpt_alpha_s[i]+=mpt_alpha_t[i];
     else mpt_alpha_b[i]+=mpt_alpha_t[i];
	 if(alpha_boosted <(frac[i+1]*TMath::Pi()) && alpha_boosted >=(frac[i-1]*TMath::Pi())) mpt_alpha_boosted_t[i] =pt *cos(phi -phi1);
     mpt_alpha_boosted[i]+=mpt_alpha_boosted_t[i];
     if(signal==1) mpt_alpha_boosted_s[i]+=mpt_alpha_boosted_t[i];
     else mpt_alpha_boosted_b[i]+=mpt_alpha_boosted_t[i];
    
     mp_alpha_t[i]=0;
     mp_alpha_boosted_t[i]=0;
	 if(alpha <(frac[i+1]*TMath::Pi()) && alpha >=(frac[i]*TMath::Pi())) mp_alpha_t[i] =vp.P()*cos(alpha);
     mp_alpha[i]+=mp_alpha_t[i];
     if(signal==1) mp_alpha_s[i]+=mp_alpha_t[i];
     else mp_alpha_b[i]+=mp_alpha_t[i];
	 if(alpha_boosted <(frac[i+1]*TMath::Pi()) && alpha_boosted >=(frac[i]*TMath::Pi())) mp_alpha_boosted_t[i] =vp_boosted.P()*cos(alpha_boosted);
     mp_alpha_boosted[i]+=mp_alpha_boosted_t[i];
     if(signal==1) mp_alpha_boosted_s[i]+=mp_alpha_boosted_t[i];
     else mp_alpha_boosted_b[i]+=mp_alpha_boosted_t[i];
	  }
	  
	     
   }
   
   
   float rmin=100;
 
   //find rmin; 
   for(int ijet=0;ijet<fjet->nref;ijet++){
    if(fabs(fjet->jteta[ijet])>2 || fjet->jtpt[ijet]<30) continue;
    float r_reco=sqrt(pow(eta-fjet->jteta[ijet],2)+pow(acos(cos(phi-fjet->jtphi[ijet])),2));
    if(r_reco<rmin)rmin=r_reco;
   }
   
   //fill in the output tree
 
   float entry[]={pt,eta,phi,rmin,cent,pt_boosted,p_boosted,eta_boosted,phi_boosted,cent,pt1,eta1,phi1,pt1_boosted,eta1_boosted,phi1_boosted,pt2,eta2,phi2,pt2_boosted,eta2_boosted,phi2_boosted,pt3,eta3,phi3,pt3_boosted,eta3_boosted,phi3_boosted,dphi,ptratio,jetspz,jetspt,jetseta,jetsrap,jetspz_12,jetspt_12,jetsrap_12,alpha,alpha_boosted,mpt_track,mpt_boosted_track,mpt_track_s,mpt_boosted_track_s,mpt_track_b,mpt_boosted_track_b};
   nt_particle->Fill(entry);
  }
   
  double partrap=vparticles.Rapidity();
  double trkrap=vtracks.Rapidity();
  double signalrap=vsignal.Rapidity();
  double backgroundrap=vbackground.Rapidity();
  
  
  nt_mpt->Fill(mpt_alpha);
  nt_mpt_boosted->Fill(mpt_alpha_boosted);
  nt_mp->Fill(mp_alpha);
  nt_mp_boosted->Fill(mp_alpha_boosted);
  
  nt_mpt_s->Fill(mpt_alpha_s);
  nt_mpt_boosted_s->Fill(mpt_alpha_boosted_s);
  nt_mp_s->Fill(mp_alpha_s);
  nt_mp_boosted_s->Fill(mp_alpha_boosted_s);
  
  nt_mpt_b->Fill(mpt_alpha_b);
  nt_mpt_boosted_b->Fill(mpt_alpha_boosted_b);
  nt_mp_b->Fill(mp_alpha_b);
  nt_mp_boosted_b->Fill(mp_alpha_boosted_b);
 
  nt_mpt_tracks->Fill(mpt_alpha_tracks);
  nt_mpt_boosted_tracks->Fill(mpt_alpha_boosted_tracks);
  nt_mp_tracks->Fill(mp_alpha_tracks);
  nt_mp_boosted_tracks->Fill(mp_alpha_boosted_tracks);
 
  float jtentry[]={cent,pt1,eta1,phi1,pt1_boosted,eta1_boosted,phi1_boosted,pt2,eta2,phi2,pt2_boosted,eta2_boosted,phi2_boosted,pt3,eta3,phi3,pt3_boosted,eta3_boosted,phi3_boosted,dphi,ptratio,jetspz,jetspt,jetseta,jetsrap,jetsphi,jetspz_12,jetspt_12,jetseta_12,jetsrap_12,jetsphi_12,jetsrap_reco,partrap,trkrap,signalrap,backgroundrap,mpt,mpt_boosted,mp,mp_boosted};
  
  nt_jet->Fill(jtentry);
 }
 

 
 outf->cd();
 // nt_track->Write();
 // nt_particle->Write();
 nt_jet->Write();
 
  nt_mpt->Write();
  nt_mpt_boosted->Write();
  nt_mp->Write();
  nt_mp_boosted->Write();
  
  nt_mpt_s->Write();
  nt_mpt_boosted_s->Write();
  nt_mp_s->Write();
  nt_mp_boosted_s->Write();
  
  nt_mpt_b->Write();
  nt_mpt_boosted_b->Write();
  nt_mp_b->Write();
  nt_mp_boosted_b->Write();
 
  nt_mpt_tracks->Write();
  nt_mpt_boosted_tracks->Write();
  nt_mp_tracks->Write();
  nt_mp_boosted_tracks->Write();
 // nt_track->Write();
 outf->Close();
 }