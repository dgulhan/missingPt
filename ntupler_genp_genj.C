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
#include "../trackEfficiency/ntupler/trackTree.C"

void ntupler_genp_genj(){
 TH1D::SetDefaultSumw2();
 
 //input file
 TString directory="root://eoscms//eos/cms/store/group/phys_heavyions/yjlee/PbPbMC2014/";
 TString infname="HiForest3_HydjetDrum_Pyquen_Dijet80_Embedded_d20140122_Track7_v1";
 
 trackTree * ftrk = new trackTree(Form("%s/%s.root",directory.Data(),infname.Data()));
 HiTree * fhi = new HiTree(Form("%s/%s.root",directory.Data(),infname.Data()));
 t * fjet = new t(Form("%s/%s.root",directory.Data(),infname.Data()));
 hi * fgen = new hi(Form("%s/%s.root",directory.Data(),infname.Data()));

 //pt bins for track efficiency correction
 int npt=5;
 double ptmin[]={0.4,1,3,8,100};
 double ptmax[]={1,3,8,100,300};

 //getting histograms for track efficiency  
 // TFile *f_eff[npt];
 // TProfile *p_eff_cent[npt]; 
 // TProfile2D *p_eff_accept[npt]; 
 // TProfile *p_eff_pt[npt]; 
 // for(int ipt=0; ipt<npt;ipt++){
   // f_eff[ipt]= new TFile(Form("../../trackEfficiency/final_hists/eff_pt%d_%d_step_cent2accept2pt2.root",(int)ptmin[ipt],(int)ptmax[ipt]));
   // p_eff_cent[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_cent");
   // p_eff_pt[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_pt");
   // p_eff_accept[ipt]=(TProfile2D*)f_eff[ipt]->Get("p_eff_acceptance");
 // }

 //output file and tree
 TFile *outf= new TFile(Form("/d01/dgulhan/missingPt/ntuples/track_ntuple_%s_jetscut_trackercut_jetframe_sb_info_pt1.root",infname.Data()),"recreate");
 std::string trackVars="pt:eta:phi:pt_boosted:p_boosted:eta_boosted:phi_boosted:cent:pt1:eta1:phi1:pt1_boosted:eta1_boosted:phi1_boosted:pt2:eta2:phi2:pt2_boosted:eta2_boosted:phi2_boosted:pt3:eta3:phi3:pt3_boosted:eta3_boosted:phi3_boosted:dphi:ptratio:jetspz:jetspt:jetseta:jetsrap:jetspz_12:jetspt_12:jetsrap_12:alpha:alpha_boosted:mpt_track:mpt_boosted_track:mpt_track_s:mpt_boosted_track_s:mpt_track_b:mpt_boosted_track_b";
 
 TNtuple *nt_track = new TNtuple("nt_track","",trackVars.data());
 std:string jetVars="cent:pt1:eta1:phi1:pt1_boosted:eta1_boosted:phi1_boosted:pt2:eta2:phi2:pt2_boosted:eta2_boosted:phi2_boosted:pt3:eta3:phi3:pt3_boosted:eta3_boosted:phi3_boosted:dphi:ptratio:jetspz:jetspt:jetseta:jetsrap:jetsphi:jetspz_12:jetspt_12:jetseta_12:jetsrap_12:jetsphi_12:jetsrap_reco:mpt:mpt_boosted:mp:mp_boosted:mpt_alpha1:mpt_alpha2:mpt_alpha3:mpt_alpha4:mpt_alpha5:mpt_alpha6:mpt_alpha7:mpt_alpha8:mpt_alpha9:mpt_alpha1_boosted:mpt_alpha2_boosted:mpt_alpha3_boosted:mpt_alpha4_boosted:mpt_alpha5_boosted:mpt_alpha6_boosted:mpt_alpha7_boosted:mpt_alpha8_boosted:mpt_alpha9_boosted:mpt_b:mpt_boosted_b:mpt_alpha1_b:mpt_alpha2_b:mpt_alpha3_b:mpt_alpha4_b:mpt_alpha5_b:mpt_alpha6_b:mpt_alpha7_b:mpt_alpha8_b:mpt_alpha9_b:mpt_alpha1_boosted_b:mpt_alpha2_boosted_b:mpt_alpha3_boosted_b:mpt_alpha4_boosted_b:mpt_alpha5_boosted_b:mpt_alpha6_boosted_b:mpt_alpha7_boosted_b:mpt_alpha8_boosted_b:mpt_alpha9_boosted_b:mpt_s:mpt_boosted_s:mpt_alpha1_s:mpt_alpha2_s:mpt_alpha3_s:mpt_alpha4_s:mpt_alpha5_s:mpt_alpha6_s:mpt_alpha7_s:mpt_alpha8_s:mpt_alpha9_s:mpt_alpha1_boosted_s:mpt_alpha2_boosted_s:mpt_alpha3_boosted_s:mpt_alpha4_boosted_s:mpt_alpha5_boosted_s:mpt_alpha6_boosted_s:mpt_alpha7_boosted_s:mpt_alpha8_boosted_s:mpt_alpha9_boosted_s:partrap:trkrap:signalrap:backgroundrap";
 
 TNtuple *nt_jet = new TNtuple("nt_jet","",jetVars.data());
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
  //loop over tracks
  // std::vector<std::pair<double, std::pair<double,std::pair<double,std::pair<double,std::pair<double,std::pair<double,std::pair<double,std::pair<double,double> > > > > > > > > jets;
  std::vector<std::pair<double, std::pair<double,double> >  >jets;
  
  int njet=0;
  TLorentzVector vjets;
  for(int ijet=0;ijet<fjet->ngen;ijet++){
   if(fabs(fjet->geneta[ijet])>2 ||fjet->genpt[ijet]<30) continue;
   // if(fjet->genmatchindex[ijet]<0){
    // jets.push_back(std::make_pair(fjet->genpt[ijet],std::make_pair(fjet->geneta[ijet], std::make_pair(fjet->genphi[ijet],std::make_pair(-100,std::make_pair(-100,std::make_pair(-100,std::make_pair(-100,std::make_pair(-100,-100)))))))));
   // }else{
    // jets.push_back(std::make_pair(fjet->genpt[ijet],std::make_pair(fjet->geneta[ijet], std::make_pair(fjet->genphi[ijet],std::make_pair(fjet->jtpt[fjet->genmatchindex[ijet]],std::make_pair(fjet->jteta[fjet->genmatchindex[ijet]],std::make_pair(fjet->jtphi[fjet->genmatchindex[ijet]],std::make_pair(fjet->neutralMax[fjet->genmatchindex[ijet]],std::make_pair(fjet->trackMax[fjet->genmatchindex[ijet]],fjet->trackN[fjet->genmatchindex[ijet]])))))))));
   // }
   jets.push_back(std::make_pair(fjet->genpt[ijet],std::make_pair(fjet->geneta[ijet],fjet->genphi[ijet])));
  
   TLorentzVector vjet;
   vjet.SetPtEtaPhiM(fjet->genpt[ijet],fjet->geneta[ijet],fjet->genphi[ijet],0);
   vjets+=vjet;
   njet++;
  } 
  
  // std::vector<std::pair<double, std::pair<double,std::pair<double,std::pair<double,std::pair<double,std::pair<double,std::pair<double,std::pair<double,double> > > > > > > > > jets_reco;
  TLorentzVector vjets_reco;
  int njet_reco;
  for(int ijet=0;ijet<fjet->nref;ijet++){
   if(fabs(fjet->jteta[ijet])>2 ||fjet->jtpt[ijet]<30) continue;
    // jets_reco.push_back(std::make_pair(fjet->jtpt[ijet],std::make_pair(fjet->jteta[ijet], std::make_pair(fjet->jtphi[ijet],std::make_pair(fjet->refpt[ijet],std::make_pair(fjet->refeta[ijet],std::make_pair(fjet->refphi[ijet],std::make_pair(fjet->neutralMax[ijet],std::make_pair(fjet->trackMax[ijet],fjet->trackN[ijet]))))))))); 
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
   // refpt1=    jets[njet-1].second.second.second.first;
   // refeta1=   jets[njet-1].second.second.second.second.first;
   // refphi1=   jets[njet-1].second.second.second.second.second.first;
   // neutralMax1=jets[njet-1].second.second.second.second.second.second.first;
   // trackMax1= jets[njet-1].second.second.second.second.second.second.second.first;
   // trackN1= jets[njet-1].second.second.second.second.second.second.second.second;
   v1.SetPtEtaPhiM(pt1,eta1,phi1,0);
   v1_boosted.SetPtEtaPhiM(pt1,eta1,phi1,0);
   vls+=v1;
   if(njet>1){
    pt2=jets[njet-2].first;
    eta2=jets[njet-2].second.first;
    phi2=jets[njet-2].second.second;
    // refpt2=jets[njet-2].second.second.second.first;
    // refeta2=jets[njet-2].second.second.second.second.first;
    // refphi2=jets[njet-2].second.second.second.second.second.first;
    // neutralMax2=jets[njet-2].second.second.second.second.second.second.first;
    // trackMax2=jets[njet-2].second.second.second.second.second.second.second.first;
    // trackN2=jets[njet-2].second.second.second.second.second.second.second.second;
    v2.SetPtEtaPhiM(pt2,eta2,phi2,0);
    v2_boosted.SetPtEtaPhiM(pt2,eta2,phi2,0);
    dphi=acos(cos(phi1-phi2));
    ptratio=pt2/pt1;
	  vls+=v2;

    if(njet>2){
     pt3=jets[njet-3].first;
     eta3=jets[njet-3].second.first;
     phi3=jets[njet-3].second.second;
     // refpt3=jets[njet-3].second.second.second.first;
     // refeta3=jets[njet-3].second.second.second.second.first;
     // refphi3=jets[njet-3].second.second.second.second.second.first;
     // neutralMax3=jets[njet-3].second.second.second.second.second.second.first;
     // trackMax3=jets[njet-3].second.second.second.second.second.second.second.first;
     // trackN3=jets[njet-3].second.second.second.second.second.second.second.second;
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
  
  double mpt=0;
  double mp=0;
  double mpt_boosted=0;
  double mp_boosted=0;
  double mpt_s=0;
  double mp_s=0;
  double mpt_boosted_s=0;
  double mp_boosted_s=0;
  double mpt_b=0;
  double mp_b=0;
  double mpt_boosted_b=0;
  double mp_boosted_b=0;
  double mpt_alpha1_boosted=0;
  double mpt_alpha2_boosted=0;
  double mpt_alpha3_boosted=0;
  double mpt_alpha4_boosted=0;
  double mpt_alpha5_boosted=0;
  double mpt_alpha6_boosted=0;
  double mpt_alpha7_boosted=0;
  double mpt_alpha8_boosted=0;
  double mpt_alpha9_boosted=0;
  double mpt_alpha1=0;
  double mpt_alpha2=0;
  double mpt_alpha3=0;
  double mpt_alpha4=0;
  double mpt_alpha5=0;
  double mpt_alpha6=0;
  double mpt_alpha7=0;
  double mpt_alpha8=0;
  double mpt_alpha9=0;
  double mpt_alpha1_boosted_s=0;
  double mpt_alpha2_boosted_s=0;
  double mpt_alpha3_boosted_s=0;
  double mpt_alpha4_boosted_s=0;
  double mpt_alpha5_boosted_s=0;
  double mpt_alpha6_boosted_s=0;
  double mpt_alpha7_boosted_s=0;
  double mpt_alpha8_boosted_s=0;
  double mpt_alpha9_boosted_s=0;
  double mpt_alpha1_s=0;
  double mpt_alpha2_s=0;
  double mpt_alpha3_s=0;
  double mpt_alpha4_s=0;
  double mpt_alpha5_s=0;
  double mpt_alpha6_s=0;
  double mpt_alpha7_s=0;
  double mpt_alpha8_s=0;
  double mpt_alpha9_s=0;
  double mpt_alpha1_boosted_b=0;
  double mpt_alpha2_boosted_b=0;
  double mpt_alpha3_boosted_b=0;
  double mpt_alpha4_boosted_b=0;
  double mpt_alpha5_boosted_b=0;
  double mpt_alpha6_boosted_b=0;
  double mpt_alpha7_boosted_b=0;
  double mpt_alpha8_boosted_b=0;
  double mpt_alpha9_boosted_b=0;
  double mpt_alpha1_b=0;
  double mpt_alpha2_b=0;
  double mpt_alpha3_b=0;
  double mpt_alpha4_b=0;
  double mpt_alpha5_b=0;
  double mpt_alpha6_b=0;
  double mpt_alpha7_b=0;
  double mpt_alpha8_b=0;
  double mpt_alpha9_b=0;
  TLorentzVector vparticles, vsignal, vbackground,vtracks;
  
  // for(int itrk=0;itrk<ftrk->nParticle;itrk++){
  for(int itrk=0;itrk<ftrk->nTrk;itrk++){
   float eta=ftrk->trkEta[itrk];
   float pt=ftrk->pPt[itrk];
   float phi=ftrk->pPhi[itrk];
   if(pt>1) continue;
   TLorentzVector vt;
   vt.SetPtEtaPhiM(pt,eta,phi,.13957018);
   vtracks+=vt;
  }
  
  for(int itrk=0;itrk<fgen->mult;itrk++){
   float eta=fgen->eta[itrk];
   if(fabs(eta)>2.4) continue; //acceptance of the tracker
   float pt=fgen->pt[itrk];
   if(pt>1) continue;
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
    mpt+=pt*cos(phi-phi1);
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
    else mpt_boosted_track_s=pt_boosted*cos(phi_boosted-phi1_boosted);
    mpt_boosted+=pt_boosted*cos(phi_boosted-phi1_boosted);
    alpha_boosted=vp_boosted.Angle(v1_boosted.Vect());
    mp_boosted+=vp_boosted.P()*cos(alpha_boosted);
	 
	double mpt_alpha1_t,mpt_alpha2_t,mpt_alpha3_t,mpt_alpha4_t,mpt_alpha5_t,mpt_alpha6_t,mpt_alpha7_t,mpt_alpha8_t,mpt_alpha9_t;
	double mpt_alpha1_boosted_t,mpt_alpha2_boosted_t,mpt_alpha3_boosted_t,mpt_alpha4_boosted_t,mpt_alpha5_boosted_t,mpt_alpha6_boosted_t,mpt_alpha7_boosted_t,mpt_alpha8_boosted_t,mpt_alpha9_boosted_t;
	mpt_alpha1_t=mpt_alpha2_t=mpt_alpha3_t=mpt_alpha4_t=mpt_alpha5_t=mpt_alpha6_t=mpt_alpha7_t=mpt_alpha8_t=mpt_alpha9_t=0;
	mpt_alpha1_boosted_t=mpt_alpha2_boosted_t=mpt_alpha3_boosted_t=mpt_alpha4_boosted_t=mpt_alpha5_boosted_t=mpt_alpha6_boosted_t=mpt_alpha7_boosted_t=mpt_alpha8_boosted_t=mpt_alpha9_boosted_t=0;
	  if(alpha <(0.1*TMath::Pi()/2) ||alpha >((1-0.1*0.5)*TMath::Pi())) mpt_alpha1_t =pt *cos(phi -phi1);
	  if(alpha <(0.2*TMath::Pi()/2) ||alpha >((1-0.2*0.5)*TMath::Pi())) mpt_alpha2_t =pt *cos(phi -phi1 );
	  if(alpha <(0.3*TMath::Pi()/2) ||alpha >((1-0.3*0.5)*TMath::Pi())) mpt_alpha3_t =pt *cos(phi -phi1 );
	  if(alpha <(0.4*TMath::Pi()/2) ||alpha >((1-0.4*0.5)*TMath::Pi())) mpt_alpha4_t =pt *cos(phi -phi1 );
	  if(alpha <(0.5*TMath::Pi()/2) ||alpha >((1-0.5*0.5)*TMath::Pi())) mpt_alpha5_t =pt *cos(phi -phi1 );
	  if(alpha <(0.6*TMath::Pi()/2) ||alpha >((1-0.6*0.5)*TMath::Pi())) mpt_alpha6_t =pt *cos(phi -phi1 );
	  if(alpha <(0.7*TMath::Pi()/2) ||alpha >((1-0.7*0.5)*TMath::Pi())) mpt_alpha7_t =pt *cos(phi -phi1 );
	  if(alpha <(0.8*TMath::Pi()/2) ||alpha >((1-0.8*0.5)*TMath::Pi())) mpt_alpha8_t =pt *cos(phi -phi1 );
	  if(alpha <(0.9*TMath::Pi()/2) ||alpha >((1-0.9*0.5)*TMath::Pi())) mpt_alpha9_t =pt *cos(phi -phi1 );
	
	  if(alpha_boosted<(0.1*TMath::Pi()/2) ||alpha_boosted>((1-0.1*0.5)*TMath::Pi())) mpt_alpha1_boosted_t=pt_boosted*cos(phi_boosted-phi1_boosted);
	  if(alpha_boosted<(0.2*TMath::Pi()/2) ||alpha_boosted>((1-0.2*0.5)*TMath::Pi())) mpt_alpha2_boosted_t=pt_boosted*cos(phi_boosted-phi1_boosted);
     if(alpha_boosted<(0.3*TMath::Pi()/2) ||alpha_boosted>((1-0.3*0.5)*TMath::Pi())) mpt_alpha3_boosted_t=pt_boosted*cos(phi_boosted-phi1_boosted);
	  if(alpha_boosted<(0.4*TMath::Pi()/2) ||alpha_boosted>((1-0.4*0.5)*TMath::Pi())) mpt_alpha4_boosted_t=pt_boosted*cos(phi_boosted-phi1_boosted);
	  if(alpha_boosted<(0.5*TMath::Pi()/2) ||alpha_boosted>((1-0.5*0.5)*TMath::Pi())) mpt_alpha5_boosted_t=pt_boosted*cos(phi_boosted-phi1_boosted);
	  if(alpha_boosted<(0.6*TMath::Pi()/2) ||alpha_boosted>((1-0.6*0.5)*TMath::Pi())) mpt_alpha6_boosted_t=pt_boosted*cos(phi_boosted-phi1_boosted);
	  if(alpha_boosted<(0.7*TMath::Pi()/2) ||alpha_boosted>((1-0.7*0.5)*TMath::Pi())) mpt_alpha7_boosted_t=pt_boosted*cos(phi_boosted-phi1_boosted);
   	 if(alpha_boosted<(0.8*TMath::Pi()/2) ||alpha_boosted>((1-0.8*0.5)*TMath::Pi()))  mpt_alpha8_boosted_t=pt_boosted*cos(phi_boosted-phi1_boosted);
	  if(alpha_boosted<(0.9*TMath::Pi()/2) ||alpha_boosted>((1-0.9*0.5)*TMath::Pi())) mpt_alpha9_boosted_t=pt_boosted*cos(phi_boosted-phi1_boosted);
       mpt_alpha1_boosted+=mpt_alpha1_boosted_t;
       mpt_alpha2_boosted+=mpt_alpha2_boosted_t;
       mpt_alpha3_boosted+=mpt_alpha3_boosted_t;
       mpt_alpha4_boosted+=mpt_alpha4_boosted_t;
       mpt_alpha5_boosted+=mpt_alpha5_boosted_t;
       mpt_alpha6_boosted+=mpt_alpha6_boosted_t;
       mpt_alpha7_boosted+=mpt_alpha7_boosted_t;
       mpt_alpha8_boosted+=mpt_alpha8_boosted_t;
       mpt_alpha9_boosted+=mpt_alpha9_boosted_t;
       mpt_alpha1+=mpt_alpha1_t;
       mpt_alpha2+=mpt_alpha2_t;
       mpt_alpha3+=mpt_alpha3_t;
       mpt_alpha4+=mpt_alpha4_t;
       mpt_alpha5+=mpt_alpha5_t;
       mpt_alpha6+=mpt_alpha6_t;
       mpt_alpha7+=mpt_alpha7_t;
       mpt_alpha8+=mpt_alpha8_t;
       mpt_alpha9+=mpt_alpha9_t;
       if(signal==1){
	       mpt_alpha1_boosted_s+=mpt_alpha1_boosted_t;
         mpt_alpha2_boosted_s+=mpt_alpha2_boosted_t;
         mpt_alpha3_boosted_s+=mpt_alpha3_boosted_t;
         mpt_alpha4_boosted_s+=mpt_alpha4_boosted_t;
         mpt_alpha5_boosted_s+=mpt_alpha5_boosted_t;
         mpt_alpha6_boosted_s+=mpt_alpha6_boosted_t;
         mpt_alpha7_boosted_s+=mpt_alpha7_boosted_t;
         mpt_alpha8_boosted_s+=mpt_alpha8_boosted_t;
         mpt_alpha9_boosted_s+=mpt_alpha9_boosted_t;
         mpt_alpha1_s+=mpt_alpha1_t;
         mpt_alpha2_s+=mpt_alpha2_t;
         mpt_alpha3_s+=mpt_alpha3_t;
         mpt_alpha4_s+=mpt_alpha4_t;
         mpt_alpha5_s+=mpt_alpha5_t;
         mpt_alpha6_s+=mpt_alpha6_t;
         mpt_alpha7_s+=mpt_alpha7_t;
         mpt_alpha8_s+=mpt_alpha8_t;
         mpt_alpha9_s+=mpt_alpha9_t;
       }else{
	   	 mpt_alpha1_boosted_b+=mpt_alpha1_boosted_t;
         mpt_alpha2_boosted_b+=mpt_alpha2_boosted_t;
         mpt_alpha3_boosted_b+=mpt_alpha3_boosted_t;
         mpt_alpha4_boosted_b+=mpt_alpha4_boosted_t;
         mpt_alpha5_boosted_b+=mpt_alpha5_boosted_t;
         mpt_alpha6_boosted_b+=mpt_alpha6_boosted_t;
         mpt_alpha7_boosted_b+=mpt_alpha7_boosted_t;
         mpt_alpha8_boosted_b+=mpt_alpha8_boosted_t;
         mpt_alpha9_boosted_b+=mpt_alpha9_boosted_t;
         mpt_alpha1_b+=mpt_alpha1_t;
         mpt_alpha2_b+=mpt_alpha2_t;
         mpt_alpha3_b+=mpt_alpha3_t;
         mpt_alpha4_b+=mpt_alpha4_t;
         mpt_alpha5_b+=mpt_alpha5_t;
         mpt_alpha6_b+=mpt_alpha6_t;
         mpt_alpha7_b+=mpt_alpha7_t;
         mpt_alpha8_b+=mpt_alpha8_t;
         mpt_alpha9_b+=mpt_alpha9_t;
	   }
   }
   //fill in the output tree
 
   // float entry[]={pt,eta,phi,pt_boosted,p_boosted,eta_boosted,phi_boosted,cent,pt1,phi1,eta1,jetspx,jetspy,jetspz,jetspt,jetseta,jetsrap,jetsphi,jetspx_12,jetspy_12,jetspz_12,jetspt_12,jetseta_12,jetsrap_12,jetsphi_12,refpt1,refeta1,refphi1,neutralMax1,trackMax1,trackN1, pt2,phi2, eta2,refpt2,refphi2,refeta2,neutralMax2,trackMax2,trackN2,pt3,phi3,eta3, refpt3,refeta3,refphi3,neutralMax3,trackMax3,trackN3, dphi, ptratio};
   float entry[]={pt,eta,phi,pt_boosted,p_boosted,eta_boosted,phi_boosted,cent,pt1,eta1,phi1,pt1_boosted,eta1_boosted,phi1_boosted,pt2,eta2,phi2,pt2_boosted,eta2_boosted,phi2_boosted,pt3,eta3,phi3,pt3_boosted,eta3_boosted,phi3_boosted,dphi,ptratio,jetspz,jetspt,jetseta,jetsrap,jetspz_12,jetspt_12,jetsrap_12,alpha,alpha_boosted,mpt_track,mpt_boosted_track,mpt_track_s,mpt_boosted_track_s,mpt_track_b,mpt_boosted_track_b};
   nt_track->Fill(entry);
  }
   
  double partrap=vparticles.Rapidity();
  double trkrap=vtracks.Rapidity();
  double signalrap=vsignal.Rapidity();
  double backgroundrap=vbackground.Rapidity();
  
  float jtentry[]={cent,pt1,eta1,phi1,pt1_boosted,eta1_boosted,phi1_boosted,pt2,eta2,phi2,pt2_boosted,eta2_boosted,phi2_boosted,pt3,eta3,phi3,pt3_boosted,eta3_boosted,phi3_boosted,dphi,ptratio,jetspz,jetspt,jetseta,jetsrap,jetsphi,jetspz_12,jetspt_12,jetseta_12,jetsrap_12,jetsphi_12,jetsrap_reco,mpt,mpt_boosted,mp,mp_boosted,mpt_alpha1,mpt_alpha2,mpt_alpha3,mpt_alpha4,mpt_alpha5,mpt_alpha6,mpt_alpha7,mpt_alpha8,mpt_alpha9,mpt_alpha1_boosted,mpt_alpha2_boosted,mpt_alpha3_boosted,mpt_alpha4_boosted,mpt_alpha5_boosted,mpt_alpha6_boosted,mpt_alpha7_boosted,mpt_alpha8_boosted,mpt_alpha9_boosted,mpt_b,mpt_boosted_b,mpt_alpha1_b,mpt_alpha2_b,mpt_alpha3_b,mpt_alpha4_b,mpt_alpha5_b,mpt_alpha6_b,mpt_alpha7_b,mpt_alpha8_b,mpt_alpha9_b,mpt_alpha1_boosted_b,mpt_alpha2_boosted_b,mpt_alpha3_boosted_b,mpt_alpha4_boosted_b,mpt_alpha5_boosted_b,mpt_alpha6_boosted_b,mpt_alpha7_boosted_b,mpt_alpha8_boosted_b,mpt_alpha9_boosted_b,mpt_s,mpt_boosted_s,mpt_alpha1_s,mpt_alpha2_s,mpt_alpha3_s,mpt_alpha4_s,mpt_alpha5_s,mpt_alpha6_s,mpt_alpha7_s,mpt_alpha8_s,mpt_alpha9_s,mpt_alpha1_boosted_s,mpt_alpha2_boosted_s,mpt_alpha3_boosted_s,mpt_alpha4_boosted_s,mpt_alpha5_boosted_s,mpt_alpha6_boosted_s,mpt_alpha7_boosted_s,mpt_alpha8_boosted_s,mpt_alpha9_boosted_s,partrap,trkrap,signalrap,backgroundrap};
  
  
  nt_jet->Fill(jtentry);
 }
 
 outf->cd();
 nt_jet->Write();
 nt_track->Write();
 outf->Close();
 }