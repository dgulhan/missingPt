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

 
void ntupler_pPb_FF(double ptmin_trk=0.5,double ptmax_trk=300){ 
TH1D::SetDefaultSumw2();
 TString algo="ak3Calo"; 
 bool doSaveTrackInfo=true;

 //input file  
 cout<<"ptmin= "<<ptmin_trk<<" ptmax= "<<ptmax_trk<<endl;
 // TString directory="root://eoscms//eos/cms/store/group/phys_heavyions/dgulhan/pp2013/";
 // TString infname="HiForest_pp_Jet80_v8_PP2013_HiForest_PromptReco_JsonPP_Jet80_PPReco_merged_forest_0"; 
    
 // TString directory="root://eoscms//eos/cms/store/group/phys_heavyions/yjlee/pPb2013/promptReco/";
 // TString infname="PA2013_HiForest_PromptReco_JSonPPb_forestv77"; 
 // TString directory="root://eoscms//eos/cms/store/caf/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt50/HiForest_v77_merged01/";
 // TString infname="pt50_HP04_prod16_v77_merged_forest_0"; 
 
 TString directory="root://eoscms//eos/cms/store/group/phys_heavyions/yjlee/pPb2013/promptReco/";
 TString infname="PA2013_HiForest_PromptReco_JSonPPb_forestv72_HLT40_HLT60"; 
 ppTrack * ftrk = new ppTrack(Form("%s/%s.root",directory.Data(),infname.Data()));
 HiTree * fhi = new HiTree(Form("%s/%s.root",directory.Data(),infname.Data()));
 skimTree_pp * fskim = new skimTree_pp(Form("%s/%s.root",directory.Data(),infname.Data()));
 // t * fjet = new t(Form("%s/%s.root",directory.Data(),infname.Data()),algo.Data());
 t * fjet = new t(Form("%s/%s.root",directory.Data(),infname.Data()),"ak3PF");
 hi * fgen = new hi(Form("%s/%s.root",directory.Data(),infname.Data()));
 
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
 double dR_upperbound[]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,3,3.5,4,4.4};
  
 //getting histograms for track efficiency correction 
 TFile *f_eff[npt];
 TProfile2D *p_eff_accept[npt]; 
 TProfile *p_eff_pt[npt]; 
 TProfile *p_eff_rmin[npt]; 
 for(int ipt=0; ipt<npt;ipt++){
   if(ipt<npt-1)f_eff[ipt]= new TFile(Form("/afs/cern.ch/user/d/dgulhan/workDir/trackEfficiency/final_hists_%s_20140414/eff_pt%d_%d_accept4pt4rmin3_%s_dogenjet0.root",algo.Data(),(int)ptmin[ipt],(int)ptmax[ipt],algo.Data()));
   else f_eff[ipt]= new TFile(Form("/afs/cern.ch/user/d/dgulhan/workDir/trackEfficiency/final_hists_%s_20140414/eff_pt%d_%d_accept4pt4rmin4_%s_dogenjet0.root",algo.Data(),(int)ptmin[ipt],(int)ptmax[ipt],algo.Data()));
   p_eff_pt[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_pt");
   p_eff_accept[ipt]=(TProfile2D*)f_eff[ipt]->Get("p_eff_acceptance");
   p_eff_rmin[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_rmin");
 }
 TFile *f_fake[npt];
 TProfile2D *p_fake_accept[npt]; 
 TProfile *p_fake_pt[npt]; 
 TProfile *p_fake_rmin[npt]; 
 for(int ipt=0; ipt<npt;ipt++){
   if(ipt<npt-1)f_fake[ipt]= new TFile(Form("/afs/cern.ch/user/d/dgulhan/workDir/trackFake_pp/final_hists_%s_20140414/fake_pt%d_%d_step_accept4pt4rmin3_%s_dogenjet0.root",algo.Data(),(int)ptmin[ipt],(int)ptmax[ipt],algo.Data()));
   else f_fake[ipt]= new TFile(Form("/afs/cern.ch/user/d/dgulhan/workDir/trackFake_pp/final_hists_%s_20140414/fake_pt%d_%d_step_accept4pt4rmin4_%s_dogenjet0.root",algo.Data(),(int)ptmin[ipt],(int)ptmax[ipt],algo.Data()));
   p_fake_pt[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_pt");
   p_fake_accept[ipt]=(TProfile2D*)f_fake[ipt]->Get("p_fake_acceptance");
   p_fake_rmin[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_rmin");
 }
 
 
 //output file and tree
 TFile *outf= new TFile(Form("ntuples_pPb_FF/small_ntuple_%s_pt%d_%d_%s.root",infname.Data(),(int)ptmin_trk,(int)ptmax_trk,algo.Data()),"recreate");
 std::string trackVars="trackselect:eff:fake:trkfake:trkstatus:weight_unfold:pt:eta:phi:rmin:pt1:eta1:phi1:pt2:eta2:phi2:pt3:eta3:phi3:dphi:ptratio:vz:pileupfilter:secondary:correction:assojet_pt:assojet_rawpt"; 
 
 TNtuple *nt_track = new TNtuple("nt_track","",trackVars.data());
  
 std::string FFVars="trackselect:eff:fake:trkfake:trkstatus:weight_unfold:pt:eta:phi:rmin:vz:pileupfilter:secondary:correction:assojet_pt:assojet_rawpt:jtpt:jtphi:jteta:rawpt:trackN:alpha"; 
 
 TNtuple *nt_FF = new TNtuple("nt_FF","",FFVars.data());
 
 std::string dijetVars="pt1:eta1:phi1:pt2:eta2:phi2:pt3:eta3:phi3:dphi:ptratio:hfp:hfm:vz:pileupfilter";
 TNtuple *nt_dijet = new TNtuple("nt_dijet","",dijetVars.data());

 std::string jetVars="jtpt:jtphi:jteta:rawpt:trackN";
 TNtuple *nt_jet = new TNtuple("nt_jet","",jetVars.data());
  
 TFile *f_multrec;
 
 f_multrec= new TFile("../JetTrack/trackMultipleRec/multreco_pp.root");
 TFile *f_secondary;
 f_secondary= new TFile("../JetTrack/trackSecondary/secondary_2d.root");
 
 TH2D * hsecondary = (TH2D*)f_secondary->Get("hpt_eta");
 TH2D * hmultrec = (TH2D*)f_multrec->Get("heta_pt");
 
 //loop over events
 int nentries = ftrk->GetEntriesFast();
 // for(int jentry=0;jentry<nentries;jentry++){
 for(int jentry=0;jentry<2000000;jentry++){
  if((jentry%1000)==0) std::cout<<jentry<<"/"<<nentries<<std::endl;
  fskim->GetEntry(jentry);
  ftrk->GetEntry(jentry);
  fhi->GetEntry(jentry);
  fjet->GetEntry(jentry);
  fgen->GetEntry(jentry);
  
  // cout<<"pcoll "<<fskim->pcollisionEventSelection<<endl;
  // cout<<"noisefilter "<<fskim->pHBHENoiseFilter<<endl;
  if(!(fskim->pPAcollisionEventSelectionPA && fskim->pHBHENoiseFilter))continue;
  
  // cout<<"passed selection"<<endl;
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
   if(fabs(fjet->jteta[ijet])>2 ||fjet->jtpt[ijet]<40) continue;
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
  TLorentzVector v1,v2,v3,vave;

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
   vave=(v1-v2);
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
  float n_sublead=0;
  float n_lead=0;
  
  TLorentzVector vtracks;

  
  int Npassingtrk=0;
  for(int itrk=0;itrk<ftrk->nTrk;itrk++){
   float eta=ftrk->trkEta[itrk];
   float pt=ftrk->trkPt[itrk];
   float phi=ftrk->trkPhi[itrk];
   
   float trkfake=ftrk->trkFake[itrk];
   float trkstatus=ftrk->trkStatus[itrk];
   float eff_pt,eff_accept,eff_rmin;
   eff_pt=eff_accept=eff_rmin=1;
   float fake_pt,fake_accept,fake_rmin;
   fake_pt=fake_accept=fake_rmin=0;
   
   if(fabs(eta)>2.4) continue;
   if(pt<ptmin_trk || pt>ptmax_trk) continue; 
   float trackselect=(ftrk->highPurity[itrk] && fabs(ftrk->trkDxy1[itrk]/ftrk->trkDxyError1[itrk])<3.0 && fabs(ftrk->trkDz1[itrk]/ftrk->trkDzError1[itrk])<3 && (ftrk->trkPtError[itrk]/ftrk->trkPt[itrk])<0.1);
   
   float weight_unfold=hweight->GetBinContent(hweight->FindBin(pt));
   
   float rmin=100;
   float assojet_pt=-99;
   float assojet_rawpt=-99;
 
   //find rmin; 
   for(int ijet=0;ijet<fjet->nref;ijet++){
    if(fabs(fjet->jteta[ijet])>2 || fjet->jtpt[ijet]<30) continue;
    float r_reco=sqrt(pow(eta-fjet->jteta[ijet],2)+pow(acos(cos(phi-fjet->jtphi[ijet])),2));
    if(r_reco<rmin){ 
     rmin=r_reco;
     assojet_pt=fjet->jtpt[ijet];
     assojet_rawpt=fjet->rawpt[ijet];
    }
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
   float fake=0;
   if(pt<300)fake=fake_accept+fake_pt+fake_rmin;
   if(fake<0) fake=0;
   
   if(eff==0){
    cout<<"zero efficiency"<<" eta="<<eta<<" pt="<<pt<<" phi="<<phi<<endl;
	  if(pt>100)eff=0.8;
	  else eff=1;
   }
   
   float multrec=0;
   multrec=hmultrec->GetBinContent(hmultrec->FindBin(pt,eta));
   float secondary=0;
   secondary=hsecondary->GetBinContent(hsecondary->FindBin(pt,eta));
   
   
   float correction = ((1-fake)*(1-secondary)/(eff));
   TLorentzVector vp;
   vp.SetPtEtaPhiM(pt,eta,phi,0.13957018);
   vtracks+=vp;
   
   
    float entry[]={trackselect,eff,fake,trkfake,trkstatus,weight_unfold,pt,eta,phi,rmin,pt1,eta1,phi1,pt2,eta2,phi2,pt3,eta3,phi3,dphi,ptratio,vz,fskim->pVertexFilterCutGplus,secondary,correction,assojet_pt};
    if(doSaveTrackInfo)nt_track->Fill(entry);
     
     // std::string FFVars="trackselect:eff:fake:trkfake:trkstatus:weight_unfold:pt:eta:phi:rmin:vz:pileupfilter:secondary:correction:assojet_pt:jtpt"; 

	 for(int ijet=0;ijet<fjet->nref;ijet++){
     if(!(fjet->jtpt[ijet]>60  && fjet->jtpt[ijet]<80)  || fabs(fjet->jteta[ijet])>2 ) continue;     
     float r_reco=sqrt(pow(eta-fjet->jteta[ijet],2)+pow(acos(cos(phi-fjet->jtphi[ijet])),2));
     TLorentzVector vjetincl;
     vjetincl.SetPtEtaPhiM(fjet->jtpt[ijet],fjet->jteta[ijet],fjet->jtphi[ijet],0);
     float alpha=vp.Angle(vjetincl.Vect());
     float FF_entry[]={trackselect,eff,fake,trkfake,trkstatus,weight_unfold,pt,eta,phi,rmin,vz,fskim->pVertexFilterCutGplus,secondary,correction,assojet_pt,assojet_rawpt,fjet->jtpt[ijet],fjet->jtphi[ijet],fjet->jteta[ijet],fjet->rawpt[ijet],fjet->trackN[ijet],alpha};
     if(r_reco<0.3) nt_FF->Fill(FF_entry);
   } 

  }
  for(int ijet=0;ijet<fjet->nref;ijet++){
    if(!(fjet->jtpt[ijet]>60  && fjet->jtpt[ijet]<80) || fabs(fjet->jteta[ijet])>2 ) continue;     
    float jtentry[]={fjet->jtpt[ijet],fjet->jtphi[ijet],fjet->jteta[ijet],fjet->rawpt[ijet],fjet->trackN[ijet]};
    nt_jet->Fill(jtentry);
   }
  float jtentry[]={pt1,eta1,phi1,pt2,eta2,phi2,pt3,eta3,phi3,dphi,ptratio,hfp,hfm,vz,fskim->pVertexFilterCutGplus};
  
 
  nt_dijet->Fill(jtentry);
 }
 

 
 outf->cd();
 if(doSaveTrackInfo){
  nt_track->Write();
 }
 nt_jet->Write();
 nt_dijet->Write();
 nt_FF->Write();
 outf->Close();
 }
