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
#include "../../trackEfficiency/v2/ntupler/trackTree.C"

void ntupler_genp_genj(){ 
 TH1D::SetDefaultSumw2();
 
 //input file
 TString directory="/mnt/hadoop/cms/store/user/dgulhan/HIMC/Jet80/Track8_Jet4_STARTHI53_LV1/merged2/";
 TString infname="HiForest_Pythia_Hydjet_Jet80_Track8_Jet4_STARTHI53_LV1_merged_forest_0"; 
  
 trackTree * ftrk = new trackTree(Form("%s/%s.root",directory.Data(),infname.Data()));
 HiTree * fhi = new HiTree(Form("%s/%s.root",directory.Data(),infname.Data()));
 t * fjet = new t(Form("%s/%s.root",directory.Data(),infname.Data()));
 hi * fgen = new hi(Form("%s/%s.root",directory.Data(),infname.Data()));

 //pt bins for track efficiency correction
 int npt=5;
 double ptmin[]={0.4,1,3,8,100};
 double ptmax[]={1,3,8,100,300};
 double frac[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.01};
 
 // getting histograms for track efficiency  

 TFile *f_eff[npt];
 TProfile *p_eff_cent[npt]; 
 TProfile2D *p_eff_accept[npt]; 
 TProfile *p_eff_pt[npt]; 
 TProfile *p_eff_rmin[npt]; 
 for(int ipt=0; ipt<npt;ipt++){
   f_eff[ipt]= new TFile(Form("../../trackEfficiency/v2/final_hists/eff_pt%d_%d_step_cent3accept3pt3rmin2.root",(int)ptmin[ipt],(int)ptmax[ipt]));
   p_eff_cent[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_cent");
   p_eff_pt[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_pt");
   p_eff_accept[ipt]=(TProfile2D*)f_eff[ipt]->Get("p_eff_acceptance");
   p_eff_rmin[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_rmin");
 }
 
 //output file and tree
 TFile *outf= new TFile(Form("track_ntuple_dR_%s_jetscut_trackercut_jetframe_v2_noptcut.root",infname.Data()),"recreate");
 std::string partVars="pt:eta:phi:rmin:cent:pt_boosted:p_boosted:eta_boosted:phi_boosted:cent:pt1:eta1:phi1:pt1_boosted:eta1_boosted:phi1_boosted:pt2:eta2:phi2:pt2_boosted:eta2_boosted:phi2_boosted:pt3:eta3:phi3:pt3_boosted:eta3_boosted:phi3_boosted:dphi:ptratio:jetspz:jetspt:jetseta:jetsrap:jetspz_12:jetspt_12:jetsrap_12:deltaR:deltaR_boosted:mpt_track:mpt_boosted_track:mpt_track_s:mpt_boosted_track_s:mpt_track_b:mpt_boosted_track_b";
 std::string trackVars="trackselect:eff:pt:eta:phi:rmin:cent:pt_boosted:p_boosted:eta_boosted:phi_boosted:cent:pt1:eta1:phi1:pt1_boosted:eta1_boosted:phi1_boosted:pt2:eta2:phi2:pt2_boosted:eta2_boosted:phi2_boosted:pt3:eta3:phi3:pt3_boosted:eta3_boosted:phi3_boosted:dphi:ptratio:jetspz:jetspt:jetseta:jetsrap:jetspz_12:jetspt_12:jetsrap_12:deltaR:deltaR_boosted:mpt_track:mpt_boosted_track:mpt_track_s:mpt_boosted_track_s:mpt_track_b:mpt_boosted_track_b";
 
 TNtuple *nt_track = new TNtuple("nt_track","",trackVars.data());
  TNtuple *nt_particle = new TNtuple("nt_particle","",partVars.data());
 
 std:string jetVars="cent:pt1:eta1:phi1:pt1_boosted:eta1_boosted:phi1_boosted:pt2:eta2:phi2:pt2_boosted:eta2_boosted:phi2_boosted:pt3:eta3:phi3:pt3_boosted:eta3_boosted:phi3_boosted:dphi:ptratio:jetspz:jetspt:jetseta:jetsrap:jetsphi:jetspz_12:jetspt_12:jetseta_12:jetsrap_12:jetsphi_12:jetsrap_reco:mpt:mpt_boosted:mp:mp_boosted:mp_deltaR_0:mp_deltaR_1:mp_deltaR_2:mp_deltaR_3:mp_deltaR_4:mp_deltaR_5:mp_deltaR_6:mp_deltaR_7:mp_deltaR_8:mp_deltaR_9:mp_deltaR_10:mp_deltaR_11:mp_deltaR_12:mp_deltaR_13:mp_deltaR_14:mp_deltaR_15:mp_deltaR_16:mp_deltaR_17:mp_deltaR_18:mp_deltaR_19:mp_deltaR_boosted_0:mp_deltaR_boosted_1:mp_deltaR_boosted_2:mp_deltaR_boosted_3:mp_deltaR_boosted_4:mp_deltaR_boosted_5:mp_deltaR_boosted_6:mp_deltaR_boosted_7:mp_deltaR_boosted_8:mp_deltaR_boosted_9:mp_deltaR_boosted_10:mp_deltaR_boosted_11:mp_deltaR_boosted_12:mp_deltaR_boosted_13:mp_deltaR_boosted_14:mp_deltaR_boosted_15:mp_deltaR_boosted_16:mp_deltaR_boosted_17:mp_deltaR_boosted_18:mp_deltaR_boosted_19:mpt_deltaR_0:mpt_deltaR_1:mpt_deltaR_2:mpt_deltaR_3:mpt_deltaR_4:mpt_deltaR_5:mpt_deltaR_6:mpt_deltaR_7:mpt_deltaR_8:mpt_deltaR_9:mpt_deltaR_10:mpt_deltaR_11:mpt_deltaR_12:mpt_deltaR_13:mpt_deltaR_14:mpt_deltaR_15:mpt_deltaR_16:mpt_deltaR_17:mpt_deltaR_18:mpt_deltaR_19:mpt_deltaR_boosted_1:mpt_deltaR_boosted_2:mpt_deltaR_boosted_3:mpt_deltaR_boosted_4:mpt_deltaR_boosted_5:mpt_deltaR_boosted_6:mpt_deltaR_boosted_7:mpt_deltaR_boosted_8:mpt_deltaR_boosted_9:mpt_deltaR_boosted_10:mpt_deltaR_boosted_11:mpt_deltaR_boosted_12:mpt_deltaR_boosted_13:mpt_deltaR_boosted_14:mpt_deltaR_boosted_15:mpt_deltaR_boosted_16:mpt_deltaR_boosted_17:mpt_deltaR_boosted_18:mpt_deltaR_boosted_19:mpt_tracks:mpt_boosted_tracks:mp_tracks:mp_boosted_tracks:mp_deltaR_tracks_0:mp_deltaR_tracks_1:mp_deltaR_tracks_2:mp_deltaR_tracks_3:mp_deltaR_tracks_4:mp_deltaR_tracks_5:mp_deltaR_tracks_6:mp_deltaR_tracks_7:mp_deltaR_tracks_8:mp_deltaR_tracks_9:mp_deltaR_tracks_10:mp_deltaR_tracks_11:mp_deltaR_tracks_12:mp_deltaR_tracks_13:mp_deltaR_tracks_14:mp_deltaR_tracks_15:mp_deltaR_tracks_16:mp_deltaR_tracks_17:mp_deltaR_tracks_18:mp_deltaR_tracks_19:mp_deltaR_boosted_tracks_0:mp_deltaR_boosted_tracks_1:mp_deltaR_boosted_tracks_2:mp_deltaR_boosted_tracks_3:mp_deltaR_boosted_tracks_4:mp_deltaR_boosted_tracks_5:mp_deltaR_boosted_tracks_6:mp_deltaR_boosted_tracks_7:mp_deltaR_boosted_tracks_8:mp_deltaR_boosted_tracks_9:mp_deltaR_boosted_tracks_10:mp_deltaR_boosted_tracks_11:mp_deltaR_boosted_tracks_12:mp_deltaR_boosted_tracks_13:mp_deltaR_boosted_tracks_14:mp_deltaR_boosted_tracks_15:mp_deltaR_boosted_tracks_16:mp_deltaR_boosted_tracks_17:mp_deltaR_boosted_tracks_18:mp_deltaR_boosted_tracks_19:mpt_deltaR_tracks_0:mpt_deltaR_tracks_1:mpt_deltaR_tracks_2:mpt_deltaR_tracks_3:mpt_deltaR_tracks_4:mpt_deltaR_tracks_5:mpt_deltaR_tracks_6:mpt_deltaR_tracks_7:mpt_deltaR_tracks_8:mpt_deltaR_tracks_9:mpt_deltaR_tracks_10:mpt_deltaR_tracks_11:mpt_deltaR_tracks_12:mpt_deltaR_tracks_13:mpt_deltaR_tracks_14:mpt_deltaR_tracks_15:mpt_deltaR_tracks_16:mpt_deltaR_tracks_17:mpt_deltaR_tracks_18:mpt_deltaR_tracks_19:mpt_deltaR_boosted_tracks_0:mpt_deltaR_boosted_tracks_1:mpt_deltaR_boosted_tracks_2:mpt_deltaR_boosted_tracks_3:mpt_deltaR_boosted_tracks_4:mpt_deltaR_boosted_tracks_5:mpt_deltaR_boosted_tracks_6:mpt_deltaR_boosted_tracks_7:mpt_deltaR_boosted_tracks_8:mpt_deltaR_boosted_tracks_9:mpt_deltaR_boosted_tracks_10:mpt_deltaR_boosted_tracks_11:mpt_deltaR_boosted_tracks_12:mpt_deltaR_boosted_tracks_13:mpt_deltaR_boosted_tracks_14:mpt_deltaR_boosted_tracks_15:mpt_deltaR_boosted_tracks_16:mpt_deltaR_boosted_tracks_17:mpt_deltaR_boosted_tracks_18:mpt_deltaR_boosted_tracks_19:partrap:trkrap:signalrap:backgroundrap";
 cout<<"2"<<endl;
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
  double mpt_tracks=0;
  double mp_tracks=0;
  double mpt_boosted_tracks=0;
  double mp_boosted_tracks=0;
  double mpt_s=0;
  double mp_s=0;
  double mpt_boosted_s=0;
  double mp_boosted_s=0;
  double mpt_b=0;
  double mp_b=0;
  double mpt_boosted_b=0;
  double mp_boosted_b=0;
  double mp_deltaR_boosted_s[20];
  double mp_deltaR_boosted_b[20];
  double mp_deltaR_boosted[20];
  double mp_deltaR_s[20];
  double mp_deltaR_b[20];
  double mp_deltaR[20];
  double mpt_deltaR_boosted_s[20];
  double mpt_deltaR_boosted_b[20];
  double mpt_deltaR_boosted[20];
  double mpt_deltaR_s[20];
  double mpt_deltaR_b[20];
  double mpt_deltaR[20];
  
  double mp_deltaR_boosted_tracks[20];
  double mp_deltaR_tracks[20];
  double mpt_deltaR_boosted_tracks[20];
  double mpt_deltaR_tracks[20];
  
  for(int i=0;i<20;i++){
   mp_deltaR[i]=0;  
   mp_deltaR_s[i]=0;  
   mp_deltaR_b[i]=0;  
   mp_deltaR_boosted[i]=0;  
   mp_deltaR_boosted_s[i]=0;  
   mp_deltaR_boosted_b[i]=0;  
   mpt_deltaR[i]=0;  
   mpt_deltaR_s[i]=0;  
   mpt_deltaR_b[i]=0;  
   mpt_deltaR_boosted[i]=0;  
   mpt_deltaR_boosted_s[i]=0;  
   mpt_deltaR_boosted_b[i]=0;  
   mp_deltaR_tracks[i]=0;  
   mp_deltaR_boosted_tracks[i]=0;  
   mpt_deltaR_tracks[i]=0;  
   mpt_deltaR_boosted_tracks[i]=0;  
  }
  
  TLorentzVector vparticles, vsignal, vbackground,vtracks;
  
  for(int itrk=0;itrk<ftrk->nTrk;itrk++){
   float eta=ftrk->trkEta[itrk];
   float pt=ftrk->pPt[itrk];
   float phi=ftrk->pPhi[itrk];
   
   float eff_pt,eff_cent,eff_accept,eff_rmin;
   eff_pt=eff_cent=eff_accept=eff_rmin=1;
   if(pt<0.4) continue;
   float trackselect=(ftrk->mtrkQual[itrk] && (ftrk->mtrkDxy1[itrk]/ftrk->mtrkDxyError1[itrk])<3.0 && (ftrk->mtrkDz1[itrk]/ftrk->mtrkDzError1[itrk])<3 && (ftrk->mtrkPtError[itrk]/ftrk->mtrkPt[itrk])<0.1);
   
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
      eff_rmin=p_eff_rmin[ipt]->GetBinContent(p_eff_rmin[ipt]->FindBin(rmin));//efficiency for rmin>3 is 1.
      // if(rmin<3)eff_rmin=p_eff_rmin[ipt]->GetBinContent(p_eff_rmin[ipt]->FindBin(rmin));//efficiency for rmin>3 is 1.
     }     
   }
  
   float eff=eff_accept*eff_cent*eff_pt*eff_rmin;
   
   TLorentzVector vp;
   vp.SetPtEtaPhiM(pt,eta,phi,0.13957018);
   vtracks+=vp;
   
   float p_boosted=-99;
   float pt_boosted=-99;
   float eta_boosted=-99;
   float phi_boosted=-99;
   float deltaR=-99;
   float deltaR_boosted=-99;
   float mpt_boosted_track=-99;
    
	 float mpt_track=pt*cos(phi-phi1);
   mpt+=pt*cos(phi-phi1);
	 deltaR=sqrt(pow(eta1-eta,2)+pow(acos(cos(phi1-phi)),2));
   mp+=vp.P()*cos(deltaR);

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
    deltaR_boosted=sqrt(pow(eta1_boosted-eta_boosted,2)+pow(acos(cos(phi1_boosted-phi_boosted)),2));
    mp_boosted_tracks+=vp_boosted.P()*cos(deltaR_boosted);
	  
    double mpt_deltaR_t[20];
    double mpt_deltaR_boosted_t[20];
    double mp_deltaR_t[20];
    double mp_deltaR_boosted_t[20];
    
   for(int i=0;i<20;i++){
     mpt_deltaR_t[i]=0;
     mpt_deltaR_boosted_t[i]=0;
	   if(deltaR <(frac[i+1]) && deltaR >=(frac[i])) mpt_deltaR_t[i] =pt*cos(phi -phi1);
     mpt_deltaR_tracks[i]+=mpt_deltaR_t[i];
	   if(deltaR_boosted <(frac[i+1]) && deltaR_boosted >=(frac[i-1])) mpt_deltaR_boosted_t[i] =pt *cos(phi -phi1);
     mpt_deltaR_boosted_tracks[i]+=mpt_deltaR_boosted_t[i];
    
     mp_deltaR_t[i]=0;
     mp_deltaR_boosted_t[i]=0;
	   if(deltaR <(frac[i+1]) && deltaR >=(frac[i])) mp_deltaR_t[i] =vp.P()*cos(deltaR);
     mp_deltaR_tracks[i]+=mp_deltaR_t[i];
	   if(deltaR_boosted <(frac[i+1]) && deltaR_boosted >=(frac[i])) mp_deltaR_boosted_t[i] =vp_boosted.P()*cos(deltaR_boosted);
     mp_deltaR_boosted_tracks[i]+=mp_deltaR_boosted_t[i];
	  }
   }
   //fill in the output tree
 
   float entry[]={trackselect,eff,pt,eta,phi,rmin,cent,pt_boosted,p_boosted,eta_boosted,phi_boosted,cent,pt1,eta1,phi1,pt1_boosted,eta1_boosted,phi1_boosted,pt2,eta2,phi2,pt2_boosted,eta2_boosted,phi2_boosted,pt3,eta3,phi3,pt3_boosted,eta3_boosted,phi3_boosted,dphi,ptratio,jetspz,jetspt,jetseta,jetsrap,jetspz_12,jetspt_12,jetsrap_12,deltaR,deltaR_boosted,mpt_track,mpt_boosted_track};
   nt_track->Fill(entry);
  }
  
  for(int itrk=0;itrk<fgen->mult;itrk++){
   float eta=fgen->eta[itrk];
   if(fabs(eta)>2.4) continue; //acceptance of the tracker
   float pt=fgen->pt[itrk];
   if(pt<0.4) continue; //acceptance of the tracker
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
   float deltaR=-99;
   float deltaR_boosted=-99;
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
	 deltaR=sqrt(pow(eta-eta1,2)+pow(acos(cos(phi-phi1)),2));
   mp+=vp.P()*cos(deltaR);
	 if(signal==1)mp_s+=vp.P()*cos(deltaR);
	 else mp_b+=vp.P()*cos(deltaR);

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
    deltaR_boosted=sqrt(pow(eta_boosted-eta1_boosted,2)+pow(acos(cos(phi_boosted-phi1_boosted)),2));
    mp_boosted+=vp_boosted.P()*cos(deltaR_boosted);
	  
    double mpt_deltaR_t[20];
    double mpt_deltaR_boosted_t[20];
    double mp_deltaR_t[20];
    double mp_deltaR_boosted_t[20];
    
   for(int i=0;i<20;i++){
     mpt_deltaR_t[i]=0;
     mpt_deltaR_boosted_t[i]=0;
	 if(deltaR <(frac[i+1]) && deltaR >=(frac[i])) mpt_deltaR_t[i] =pt*cos(phi -phi1);
     mpt_deltaR[i]+=mpt_deltaR_t[i];
     if(signal==1) mpt_deltaR_s[i]+=mpt_deltaR_t[i];
     else mpt_deltaR_b[i]+=mpt_deltaR_t[i];
	 if(deltaR_boosted <(frac[i+1]) && deltaR_boosted >=(frac[i-1])) mpt_deltaR_boosted_t[i] =pt *cos(phi -phi1);
     mpt_deltaR_boosted[i]+=mpt_deltaR_boosted_t[i];
     if(signal==1) mpt_deltaR_boosted_s[i]+=mpt_deltaR_boosted_t[i];
     else mpt_deltaR_boosted_b[i]+=mpt_deltaR_boosted_t[i];
    
     mp_deltaR_t[i]=0;
     mp_deltaR_boosted_t[i]=0;
	 if(deltaR <(frac[i+1]) && deltaR >=(frac[i])) mp_deltaR_t[i] =vp.P()*cos(deltaR);
     mp_deltaR[i]+=mp_deltaR_t[i];
     if(signal==1) mp_deltaR_s[i]+=mp_deltaR_t[i];
     else mp_deltaR_b[i]+=mp_deltaR_t[i];
	 if(deltaR_boosted <(frac[i+1]) && deltaR_boosted >=(frac[i])) mp_deltaR_boosted_t[i] =vp_boosted.P()*cos(deltaR_boosted);
     mp_deltaR_boosted[i]+=mp_deltaR_boosted_t[i];
     if(signal==1) mp_deltaR_boosted_s[i]+=mp_deltaR_boosted_t[i];
     else mp_deltaR_boosted_b[i]+=mp_deltaR_boosted_t[i];
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
 
   float entry[]={pt,eta,phi,rmin,cent,pt_boosted,p_boosted,eta_boosted,phi_boosted,cent,pt1,eta1,phi1,pt1_boosted,eta1_boosted,phi1_boosted,pt2,eta2,phi2,pt2_boosted,eta2_boosted,phi2_boosted,pt3,eta3,phi3,pt3_boosted,eta3_boosted,phi3_boosted,dphi,ptratio,jetspz,jetspt,jetseta,jetsrap,jetspz_12,jetspt_12,jetsrap_12,deltaR,deltaR_boosted,mpt_track,mpt_boosted_track,mpt_track_s,mpt_boosted_track_s,mpt_track_b,mpt_boosted_track_b};
   nt_particle->Fill(entry);
  }
   
  double partrap=vparticles.Rapidity();
  double trkrap=vtracks.Rapidity();
  double signalrap=vsignal.Rapidity();
  double backgroundrap=vbackground.Rapidity();
  
  float jtentry[]={cent,pt1,eta1,phi1,pt1_boosted,eta1_boosted,phi1_boosted,pt2,eta2,phi2,pt2_boosted,eta2_boosted,phi2_boosted,pt3,eta3,phi3,pt3_boosted,eta3_boosted,phi3_boosted,dphi,ptratio,jetspz,jetspt,jetseta,jetsrap,jetsphi,jetspz_12,jetspt_12,jetseta_12,jetsrap_12,jetsphi_12,jetsrap_reco,mpt,mpt_boosted,mp,mp_boosted,mp_deltaR[0],mp_deltaR[1],mp_deltaR[2],mp_deltaR[3],mp_deltaR[4],mp_deltaR[5],mp_deltaR[6],mp_deltaR[7],mp_deltaR[8],mp_deltaR[9],mp_deltaR[10],mp_deltaR[11],mp_deltaR[12],mp_deltaR[13],mp_deltaR[14],mp_deltaR[15],mp_deltaR[16],mp_deltaR[17],mp_deltaR[18],mp_deltaR[19],mp_deltaR_boosted[0],mp_deltaR_boosted[1],mp_deltaR_boosted[2],mp_deltaR_boosted[3],mp_deltaR_boosted[4],mp_deltaR_boosted[5],mp_deltaR_boosted[6],mp_deltaR_boosted[7],mp_deltaR_boosted[8],mp_deltaR_boosted[9],mp_deltaR_boosted[10],mp_deltaR_boosted[11],mp_deltaR_boosted[12],mp_deltaR_boosted[13],mp_deltaR_boosted[14],mp_deltaR_boosted[15],mp_deltaR_boosted[16],mp_deltaR_boosted[17],mp_deltaR_boosted[18],mp_deltaR_boosted[19],mpt_deltaR[0],mpt_deltaR[1],mpt_deltaR[2],mpt_deltaR[3],mpt_deltaR[4],mpt_deltaR[5],mpt_deltaR[6],mpt_deltaR[7],mpt_deltaR[8],mpt_deltaR[9],mpt_deltaR[10],mpt_deltaR[11],mpt_deltaR[12],mpt_deltaR[13],mpt_deltaR[14],mpt_deltaR[15],mpt_deltaR[16],mpt_deltaR[17],mpt_deltaR[18],mpt_deltaR[19],mpt_deltaR_boosted[0],mpt_deltaR_boosted[1],mpt_deltaR_boosted[2],mpt_deltaR_boosted[3],mpt_deltaR_boosted[4],mpt_deltaR_boosted[5],mpt_deltaR_boosted[6],mpt_deltaR_boosted[7],mpt_deltaR_boosted[8],mpt_deltaR_boosted[9],mpt_deltaR_boosted[10],mpt_deltaR_boosted[11],mpt_deltaR_boosted[12],mpt_deltaR_boosted[13],mpt_deltaR_boosted[14],mpt_deltaR_boosted[15],mpt_deltaR_boosted[16],mpt_deltaR_boosted[17],mpt_deltaR_boosted[18],mpt_deltaR_boosted[19],mpt_tracks,mpt_boosted_tracks,mp_tracks,mp_boosted_tracks,mp_deltaR_tracks[0],mp_deltaR_tracks[1],mp_deltaR_tracks[2],mp_deltaR_tracks[3],mp_deltaR_tracks[4],mp_deltaR_tracks[5],mp_deltaR_tracks[6],mp_deltaR_tracks[7],mp_deltaR_tracks[8],mp_deltaR_tracks[9],mp_deltaR_tracks[10],mp_deltaR_tracks[11],mp_deltaR_tracks[12],mp_deltaR_tracks[13],mp_deltaR_tracks[14],mp_deltaR_tracks[15],mp_deltaR_tracks[16],mp_deltaR_tracks[17],mp_deltaR_tracks[18],mp_deltaR_tracks[19],mp_deltaR_boosted_tracks[0],mp_deltaR_boosted_tracks[1],mp_deltaR_boosted_tracks[2],mp_deltaR_boosted_tracks[3],mp_deltaR_boosted_tracks[4],mp_deltaR_boosted_tracks[5],mp_deltaR_boosted_tracks[6],mp_deltaR_boosted_tracks[7],mp_deltaR_boosted_tracks[8],mp_deltaR_boosted_tracks[9],mp_deltaR_boosted_tracks[10],mp_deltaR_boosted_tracks[11],mp_deltaR_boosted_tracks[12],mp_deltaR_boosted_tracks[13],mp_deltaR_boosted_tracks[14],mp_deltaR_boosted_tracks[15],mp_deltaR_boosted_tracks[16],mp_deltaR_boosted_tracks[17],mp_deltaR_boosted_tracks[18],mp_deltaR_boosted_tracks[19],mpt_deltaR_tracks[0],mpt_deltaR_tracks[1],mpt_deltaR_tracks[2],mpt_deltaR_tracks[3],mpt_deltaR_tracks[4],mpt_deltaR_tracks[5],mpt_deltaR_tracks[6],mpt_deltaR_tracks[7],mpt_deltaR_tracks[8],mpt_deltaR_tracks[9],mpt_deltaR_tracks[10],mpt_deltaR_tracks[11],mpt_deltaR_tracks[12],mpt_deltaR_tracks[13],mpt_deltaR_tracks[14],mpt_deltaR_tracks[15],mpt_deltaR_tracks[16],mpt_deltaR_tracks[17],mpt_deltaR_tracks[18],mpt_deltaR_tracks[19],mpt_deltaR_boosted_tracks[0],mpt_deltaR_boosted_tracks[1],mpt_deltaR_boosted_tracks[2],mpt_deltaR_boosted_tracks[3],mpt_deltaR_boosted_tracks[4],mpt_deltaR_boosted_tracks[5],mpt_deltaR_boosted_tracks[6],mpt_deltaR_boosted_tracks[7],mpt_deltaR_boosted_tracks[8],mpt_deltaR_boosted_tracks[9],mpt_deltaR_boosted_tracks[10],mpt_deltaR_boosted_tracks[11],mpt_deltaR_boosted_tracks[12],mpt_deltaR_boosted_tracks[13],mpt_deltaR_boosted_tracks[14],mpt_deltaR_boosted_tracks[15],mpt_deltaR_boosted_tracks[16],mpt_deltaR_boosted_tracks[17],mpt_deltaR_boosted_tracks[18],mpt_deltaR_boosted_tracks[19],partrap,trkrap,signalrap,backgroundrap};
  
  nt_jet->Fill(jtentry);
 }
 

 
 outf->cd();
 nt_track->Write();
 nt_particle->Write();
 nt_jet->Write();
 // nt_track->Write();
 outf->Close();
 }