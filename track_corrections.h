#include "utilities.h"

class Track_corr{
 private:
  int                 table;
  int                 npt; 
  static const int    npt_pp = 4;
  static const int    npt_PbPb = 29;
  
  double              ptmin[29];
  double              ptmax[29];

  static const double ptmin_PbPb[29];
  static const double ptmax_PbPb[29];

  static const double ptmin_pp[4];
  static const double ptmax_pp[4];
  
  static const int   cent_min[29];
  static const int   cent_max[29];
 
  TFile              *f_eff[29];
  TProfile           *p_eff_cent[29]; 
  TProfile2D         *p_eff_accept[29]; 
  TProfile           *p_eff_pt[29]; 
  TProfile           *p_eff_rmin[29]; 

  TFile              *f_fake[29];
  TProfile           *p_fake_cent[29]; 
  TProfile2D         *p_fake_accept[29]; 
  TProfile           *p_fake_pt[29]; 
  TProfile           *p_fake_rmin[29]; 
 
  TFile              *f_multrec;
  TFile              *f_secondary;
  TH2D               *hsecondary; 
  TH2D               *hmultrec ;
  
 public:
  Track_corr(int table){ //constructor opens the files and gets the histograms 
  //0 for PbPb
  //1 for pp with pp tracking
  //2 for pp with HI tracking 
   this->table = table;
   
   switch(table){
    case 0:
     npt=npt_PbPb;
     break;
    case 1: case 2:
     npt=npt_pp;
     break;
    default:
     cout << "Invalid option for tracking table" << endl;
     cout << "0 for PbPb, "
     << "1 for pp with pp tracking, "
     << "2 for pp with HI tracking" << endl;
     break;
   }
 
   for(int ipt=0; ipt<npt;ipt++){
   
    switch(table){
     case 0:
      ptmin[ipt]=ptmin_PbPb[ipt];
      ptmax[ipt]=ptmax_PbPb[ipt];

      f_eff[ipt]= new TFile(Form("TrackCorrectionsPbPb/eff/eff_pt%d_%d_cent%d_%d.root",
      (int)(ptmin[ipt]*100),(int)(ptmax[ipt]*100),
      (int)(0.5*cent_min[ipt]),(int)(0.5*cent_max[ipt])));
      f_fake[ipt]= new TFile(Form("TrackCorrectionsPbPb/fake/fake_pt%d_%d_cent%d_%d.root",
      (int)(ptmin[ipt]*100),(int)(ptmax[ipt]*100),
      (int)(0.5*cent_min[ipt]),(int)(0.5*cent_max[ipt])));    
      
      f_multrec=NULL;
      f_secondary=NULL;
      hsecondary=NULL;
      hmultrec=NULL;
      
      break;
     case 1:
      ptmin[ipt]=ptmin_pp[ipt];
      ptmax[ipt]=ptmax_pp[ipt];

      f_eff[ipt]= new TFile(Form("TrackCorrections_pp_ppTracking/eff/eff_pt%d_%d_accept4pt4rmin3_ak3Calo_dogenjet0.root",
      (int)ptmin[ipt],(int)ptmax[ipt]));
      f_fake[ipt]= new TFile(Form("TrackCorrections_pp_ppTracking/fake/fake_pt%d_%d_step_accept4pt4rmin3_ak3Calo_dogenjet0.root",
      (int)ptmin[ipt],(int)ptmax[ipt]));


      f_multrec= new TFile("TrackCorrections_pp_ppTracking/multreco/multreco_pp.root");
      f_secondary= new TFile("TrackCorrections_pp_ppTracking/secondary/secondary_pp.root");
      hsecondary = (TH2D*)f_secondary->Get("hpt_eta");
      hmultrec = (TH2D*)f_multrec->Get("hpt_eta");
      
      break;
     case 2:
      ptmin[ipt]=ptmin_pp[ipt];
      ptmax[ipt]=ptmax_pp[ipt];

      f_eff[ipt]= new TFile(Form("TrackCorrections_pp_HITracking/eff/eff_pt%d_%d_accept4pt4rmin3_ak3Calo_dogenjet0.root",
      (int)ptmin[ipt],(int)ptmax[ipt]));
      f_fake[ipt]= new TFile(Form("TrackCorrections_pp_HITracking/fake/fake_pt%d_%d_step_accept4pt4rmin3_ak3Calo_dogenjet0.root",
      (int)ptmin[ipt],(int)ptmax[ipt]));  
      
      f_multrec= new TFile("TrackCorrections_pp_HITracking/multreco/multreco_pp_hireco.root");
      f_secondary= new TFile("TrackCorrections_pp_HITracking/secondary/secondary_pp_HIReco.root");
      hsecondary = (TH2D*)f_secondary->Get("hpt_eta");
      hmultrec = (TH2D*)f_multrec->Get("hpt_eta");
      
      break;
     default: 
      cout << "Invalid option for tracking table" << endl;
      cout << "0 for PbPb, "
      << "1 for pp with pp tracking, "
      << "2 for pp with HI tracking" << endl;
      break;
    }
   }
  
   for(int ipt=0; ipt<npt;ipt++){
    if(table==0) p_eff_cent[ipt] = (TProfile*)f_eff[ipt]->Get("p_eff_cent");
    else p_eff_cent[ipt] = NULL;
    p_eff_pt[ipt] = (TProfile*)f_eff[ipt]->Get("p_eff_pt");
    p_eff_accept[ipt] = (TProfile2D*)f_eff[ipt]->Get("p_eff_acceptance");    
    p_eff_rmin[ipt] = (TProfile*)f_eff[ipt]->Get("p_eff_rmin");
    if(table==0) p_fake_cent[ipt] = (TProfile*)f_fake[ipt]->Get("p_fake_cent");
    else p_fake_cent[ipt] = NULL;
    p_fake_pt[ipt] = (TProfile*)f_fake[ipt]->Get("p_fake_pt");
    p_fake_accept[ipt] = (TProfile2D*)f_fake[ipt]->Get("p_fake_acceptance");
    p_fake_rmin[ipt] = (TProfile*)f_fake[ipt]->Get("p_fake_rmin");
   }
  }
  
  double get_correction(double pt, double eta, double phi, double rmin, double cent=0){ //cent is not used if table is 1 or 2
   
   double correction = 1;
   
   double eff_pt,eff_cent,eff_accept,eff_rmin;
   eff_pt = eff_cent = eff_accept = eff_rmin = 1;
   double fake_pt,fake_cent,fake_accept,fake_rmin;
   fake_pt = fake_cent = fake_accept = fake_rmin = 0;
  
   for(int ipt=0;ipt<npt;ipt++){
    if(pt>=ptmin[ipt] && pt<ptmax[ipt] && cent>=cent_min[ipt] && cent<cent_max[ipt]){
      eff_pt = p_eff_pt[ipt]->GetBinContent(p_eff_pt[ipt]->FindBin(pt));
      if(table==0) eff_cent = p_eff_cent[ipt]->GetBinContent(p_eff_cent[ipt]->FindBin(cent));
      eff_accept = p_eff_accept[ipt]->GetBinContent(p_eff_accept[ipt]->GetXaxis()->FindBin(phi),p_eff_accept[ipt]->GetYaxis()->FindBin(eta));
      if(rmin<5)eff_rmin = p_eff_rmin[ipt]->GetBinContent(p_eff_rmin[ipt]->FindBin(rmin));
      else eff_rmin = p_eff_rmin[ipt]->GetBinContent(p_eff_rmin[ipt]->FindBin(5));
     }     
   } 
   
   double eff=1;
   eff=eff_accept*eff_cent*eff_pt*eff_rmin;
   if(eff==0){
    if(pt>100)eff=0.8;
	  else eff=1;
   }
      
   for(int ipt=0;ipt<npt;ipt++){
    if(pt>=ptmin[ipt] && pt<ptmax[ipt] && cent>=cent_min[ipt] && cent<cent_max[ipt]){
      fake_pt = p_fake_pt[ipt]->GetBinContent(p_fake_pt[ipt]->FindBin(pt));
      if(table==0) fake_cent = p_fake_cent[ipt]->GetBinContent(p_fake_cent[ipt]->FindBin(cent));
      fake_accept = p_fake_accept[ipt]->GetBinContent(p_fake_accept[ipt]->GetXaxis()->FindBin(phi),p_fake_accept[ipt]->GetYaxis()->FindBin(eta));
 	    if(rmin<5) fake_rmin = p_fake_rmin[ipt]->GetBinContent(p_fake_rmin[ipt]->FindBin(rmin));
 	    else fake_rmin = p_fake_rmin[ipt]->GetBinContent(p_fake_rmin[ipt]->FindBin(5));
     }     
   }
   
   float fake=0;
   if(pt<100)fake = fake_accept+fake_cent+fake_pt+fake_rmin;
   
   if(table==0){
    correction = (1-fake)/eff;
   }else{
    double multrec = hmultrec->GetBinContent(hmultrec->FindBin(pt,eta));
    double secondary = hsecondary->GetBinContent(hsecondary->FindBin(pt,eta));
    correction = ((1-secondary)*(1-fake))/(eff*(1+multrec));
   }
   
   return correction;
  }
  
  void delete_table(){
   delete[] f_eff;
   delete[] p_eff_cent; 
   delete[] p_eff_accept; 
   delete[] p_eff_pt; 
   delete[] p_eff_rmin; 

   delete[] f_fake;
   delete[] p_fake_cent; 
   delete[] p_fake_accept; 
   delete[] p_fake_pt; 
   delete[] p_fake_rmin; 
 
   delete f_multrec;
   delete f_secondary;
   delete hsecondary; 
   delete hmultrec ;
  }
};

 const int Track_corr::cent_min[29] =      {  0,   20, 40,   60, 100,   0,  20,  40,  60, 100,   0,  20,  40,  60, 100,  0, 20, 40, 60,100, 0,20,40, 60,100, 0,20, 40,  0};
 const int Track_corr::cent_max[29] =      { 20,   40, 60,  100, 200,  20,  40,  60, 100, 200,  20,  40,  60, 100, 200, 20, 40, 60,100,200,20,40,60,100,200,20,40,200,200};
 
 const double Track_corr::ptmin_PbPb[29] = { 0.5, 0.5, 0.5, 0.5, 0.5,0.55,0.55,0.55,0.55,0.55,0.65,0.65,0.65,0.65,0.65,0.8,0.8,0.8,0.8,0.8, 1, 1, 1,  1,  1, 3, 3,  3,  8};
 const double Track_corr::ptmax_PbPb[29] = {0.55,0.55,0.55,0.55,0.55,0.65,0.65,0.65,0.65,0.65, 0.8, 0.8, 0.8, 0.8, 0.8,  1,  1,  1,  1,  1, 3, 3, 3,  3,  3, 8, 8,  8,300};

 const double Track_corr::ptmin_pp[4] = {0.5, 1, 3,   8};
 const double Track_corr::ptmax_pp[4] = {  1, 3,  8,300};