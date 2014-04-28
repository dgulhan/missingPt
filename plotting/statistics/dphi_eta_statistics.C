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

void dphi_eta_statistics(){

   TFile *f = TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_PbPb_mptonly_20140323_v3/full_ntuple_hiForest_Jet80or95_GR_R_53_LV6_12Mar2014_0000CET_Track8_Jet21_0_pt0_300_akVs3Calo.root")); 
   TFile * f_ref= TFile::Open(Form("/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntuples_pp_mptonly_20140323_v3/full_ntuple_pt80_pp2013_P01_prod22_v81_merged_forest_0_pt0_300_ak3Calo.root")); 
  
  cout<<1<<endl;
  TTree * t_jet = (TTree*)f->Get("nt_jet");
  TTree * t_jet_ref = (TTree*)f_ref->Get("nt_jet");
   
   
  double   nentries =   t_jet->GetEntries("pt1>120 && pt2>50 && abs(vz)<15");
  double   nentries_ref =   t_jet_ref->GetEntries("pt1>120 && pt2>50 && abs(vz)<15");
  
  TH2D * prof1 = new TH2D("prof1",";#Delta#phi;#eta_{1}",40,2,TMath::Pi(),20,0,2);
  TH2D * prof2 = new TH2D("prof2",";#Delta#phi;#eta_{2}",40,2,TMath::Pi(),20,0,2);
  
  TH2D * prof1_ref = new TH2D("prof1_ref",";#Delta#phi;#eta_{1}",40,2,TMath::Pi(),20,0,2);
  TH2D * prof2_ref = new TH2D("prof2_ref",";#Delta#phi;#eta_{2}",40,2,TMath::Pi(),20,0,2);
  
  t_jet->Draw("eta1:dphi>>prof1","pt1>120 && pt2>50 && abs(vz)<15");
  t_jet->Draw("eta2:dphi>>prof2","pt1>120 && pt2>50 && abs(vz)<15");
  
  t_jet_ref->Draw("eta1:dphi>>prof1_ref","pt1>120 && pt2>50 && abs(vz)<15");
  t_jet_ref->Draw("eta2:dphi>>prof2_ref","pt1>120 && pt2>50 && abs(vz)<15");
  
  prof1->Scale(1/nentries);
  prof2->Scale(1/nentries);
  prof1_ref->Scale(1/nentries_ref);
  prof2_ref->Scale(1/nentries_ref);
  
  TCanvas * c1 = new TCanvas("c1","",600,600);
  c1->SetRightMargin(0.2);
  prof1->Draw("colz");
  c1->SaveAs("eta1_dphi.png");
  
  TCanvas * c2 = new TCanvas("c2","",600,600);
  c2->SetRightMargin(0.2);
  prof2->Draw("colz");
  c2->SaveAs("eta2_dphi.png");
  
  TCanvas * c3 = new TCanvas("c3","",600,600);
  c3->SetRightMargin(0.2);

  prof1_ref->Draw("colz");
  c3->SaveAs("eta1_dphi_ref.png");
  
  TCanvas * c4 = new TCanvas("c4","",600,600);
  c4->SetRightMargin(0.2);
  prof2_ref->Draw("colz");
  c4->SaveAs("eta2_dphi_ref.png");
  
  cout<<"|eta1,2|<0.5 "<<"dphi>9*TMath::Pi()/10 "<<t_jet->GetEntries("dphi>9*TMath::Pi()/10 && pt1>120 && pt2>50 && abs(vz)<15 && abs(eta1)<0.5 && abs(eta2)<0.5")<<"/"<<nentries<<endl;
  cout<<"|eta1,2|<0.5 "<<"dphi>9*TMath::Pi()/10 "<<t_jet_ref->GetEntries("dphi>9*TMath::Pi()/10 && pt1>120 && pt2>50 && abs(vz)<15 && abs(eta1)<0.5 && abs(eta2)<0.5")<<"/"<<nentries_ref<<endl;
  
  cout<<"|eta1,2|<0.6 "<<"dphi>9*TMath::Pi()/10 "<<t_jet->GetEntries("dphi>9*TMath::Pi()/10 && pt1>120 && pt2>50 && abs(vz)<15 && abs(eta1)<0.6 && abs(eta2)<0.6")<<"/"<<nentries<<endl;
  cout<<"|eta1,2|<0.6 "<<"dphi>9*TMath::Pi()/10 "<<t_jet_ref->GetEntries("dphi>9*TMath::Pi()/10 && pt1>120 && pt2>50 && abs(vz)<15 && abs(eta1)<0.6 && abs(eta2)<0.6")<<"/"<<nentries_ref<<endl;
  
  cout<<"|eta1,2|<0.7 "<<"dphi>9*TMath::Pi()/10 "<<t_jet->GetEntries("dphi>9*TMath::Pi()/10 && pt1>120 && pt2>50 && abs(vz)<15 && abs(eta1)<0.7 && abs(eta2)<0.7")<<"/"<<nentries<<endl;
  cout<<"|eta1,2|<0.7 "<<"dphi>9*TMath::Pi()/10 "<<t_jet_ref->GetEntries("dphi>9*TMath::Pi()/10 && pt1>120 && pt2>50 && abs(vz)<15 && abs(eta1)<0.7 && abs(eta2)<0.7")<<"/"<<nentries_ref<<endl;
  
  cout<<"|eta1,2|<0.5 "<<"dphi>7*TMath::Pi()/8 "<<t_jet->GetEntries("dphi>7*TMath::Pi()/8 && pt1>120 && pt2>50 && abs(vz)<15 && abs(eta1)<0.5 && abs(eta2)<0.5")<<"/"<<nentries<<endl;
  cout<<"|eta1,2|<0.5 "<<"dphi>7*TMath::Pi()/8 "<<t_jet_ref->GetEntries("dphi>7*TMath::Pi()/8 && pt1>120 && pt2>50 && abs(vz)<15 && abs(eta1)<0.5 && abs(eta2)<0.5")<<"/"<<nentries_ref<<endl;
  
  cout<<"|eta1,2|<0.6 "<<"dphi>7*TMath::Pi()/8 "<<t_jet->GetEntries("dphi>7*TMath::Pi()/8 && pt1>120 && pt2>50 && abs(vz)<15 && abs(eta1)<0.6 && abs(eta2)<0.6")<<"/"<<nentries<<endl;
  cout<<"|eta1,2|<0.6 "<<"dphi>7*TMath::Pi()/8 "<<t_jet_ref->GetEntries("dphi>7*TMath::Pi()/8 && pt1>120 && pt2>50 && abs(vz)<15 && abs(eta1)<0.6 && abs(eta2)<0.6")<<"/"<<nentries_ref<<endl;
  
  cout<<"|eta1,2|<0.7 "<<"dphi>7*TMath::Pi()/8 "<<t_jet->GetEntries("dphi>7*TMath::Pi()/8 && pt1>120 && pt2>50 && abs(vz)<15 && abs(eta1)<0.7 && abs(eta2)<0.7")<<"/"<<nentries<<endl;
  cout<<"|eta1,2|<0.7 "<<"dphi>7*TMath::Pi()/8 "<<t_jet_ref->GetEntries("dphi>7*TMath::Pi()/8 && pt1>120 && pt2>50 && abs(vz)<15 && abs(eta1)<0.7 && abs(eta2)<0.7")<<"/"<<nentries_ref<<endl;
}