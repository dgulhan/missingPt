#include <iostream>
#include "TCanvas.h"
#include "TError.h"
#include "TPad.h"
#include "TString.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TF1.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TLatex.h"
#include "TString.h"  
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCut.h"
#include "TNtuple.h"
#include "TProfile2D.h"
#include "TMath.h"
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
void drawPatch(float x1, float y1, float x2, float y2){
  TLegend *t1=new TLegend(x1,y1,x2,y2);
  t1->SetFillColor(kWhite);
  t1->SetBorderSize(0);
  t1->SetFillStyle(1001);
  t1->Draw("");
}

void drawText(const char *text, float xp, float yp, int size=22){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(43);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw("same");
}


void compare_Calo_PF(){
TH1D::SetDefaultSumw2();

TFile *f[2];
f[0] =new TFile("/d01/dgulhan/ntuple_HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_0_50k_akPuPF.root");
f[1] =new TFile("/d01/dgulhan/track_ntuple_HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_0_jetscut_trackercut_jetframe_v2_pt0_100_akPu3Calo.root");
 
int ncent=6;
int cent_min[]={0,20,40,60,100,140};
int cent_max[]={20,40,60,100,140,200};

TTree *t[2];
TTree *t_jet[2];
TH1D * h[2];
TH1D * h_ref[2];
TH1D * h_binned[2][ncent];
TH1D * h_binned_ref[2][ncent];
TH1D * h_wrong[2];
TH1D * h_wrong_ref[2];
double n_dijet[2];
double n_dijet_binned[2][ncent];

for(int ialgo=0;ialgo<2;ialgo++){
 t[ialgo]=(TTree*)f[ialgo]->Get("nt_track");
 t_jet[ialgo]=(TTree*)f[ialgo]->Get("nt_jet");
 
 n_dijet[ialgo]=t_jet[ialgo]->GetEntries("pt1>120 && pt2>30 && dphi>2.094");
 
 h[ialgo]=new TH1D(Form("h_%d",ialgo),";leading jet #xi;(1/N_{dijet})dN/d#xi",20,0,5.5);
 h_ref[ialgo]=new TH1D(Form("h_ref_%d",ialgo),";leading jet #xi;(1/N_{dijet})dN/d#xi",20,0,5.5);
 t[ialgo]->Draw(Form("-log((pt/(pt1*cosh(-eta1)*cosh(-eta1)))*(cos(phi)*cos(phi1)+sin(phi)*sin(phi1)+sinh(-eta1)*sinh(eta)))>>h_ref_%d",ialgo),"abs(eta)<2.4 && pt1>120 && pt2>30 && dphi>2.094  && sqrt(pow(acos(cos(phi-phi1)),2)+pow(eta+eta1,2))<0.3 && pt>1 && abs(eta1)>0.3");
 t[ialgo]->Draw(Form("-log((pt/(pt1*cosh(eta1)*cosh(eta1)))*(cos(phi)*cos(phi1)+sin(phi)*sin(phi1)+sinh(eta1)*sinh(eta)))>>h_%d",ialgo),"abs(eta)<2.4 && pt1>120 && pt2>30 && dphi>2.094 && sqrt(pow(acos(cos(phi-phi1)),2)+pow(eta-eta1,2))<0.3 && pt>1 && abs(eta1)>0.3");
 
 h_wrong[ialgo]=new TH1D(Form("h_wrong_%d",ialgo),";leading jet #xi;(1/N_{dijet})dN/d#xi",20,0,5.5);
 h_wrong_ref[ialgo]=new TH1D(Form("h_wrong_ref_%d",ialgo),";leading jet #xi;(1/N_{dijet})dN/d#xi",20,0,5.5);
 t[ialgo]->Draw(Form("-log(pt*cos(sqrt(pow(acos(cos(phi-phi1)),2)+pow(eta+eta1,2)))*cosh(eta)/cosh(eta1)/pt1)>>h_wrong_ref_%d",ialgo),"abs(eta)<2.4 && pt1>120 && pt2>30 && dphi>2.094  && sqrt(pow(acos(cos(phi-phi1)),2)+pow(eta+eta1,2))<0.3 && pt>1 && abs(eta1)>0.3");
 t[ialgo]->Draw(Form("-log(pt*cos(sqrt(pow(acos(cos(phi-phi1)),2)+pow(eta-eta1,2)))*cosh(eta)/cosh(eta1)/pt1)>>h_wrong_%d",ialgo),"abs(eta)<2.4 && pt1>120 && pt2>30 && dphi>2.094 && sqrt(pow(acos(cos(phi-phi1)),2)+pow(eta-eta1,2))<0.3 && pt>1 && abs(eta1)>0.3");
 
 for(int icent=0;icent<ncent;icent++){
  h_binned[ialgo][icent]= new TH1D(Form("h_binned_%d_%d",ialgo,icent),";leading jet #xi;(1/N_{dijet})dN/d#xi",10,0,5.5);
  h_binned_ref[ialgo][icent]= new TH1D(Form("h_binned_ref_%d_%d",ialgo,icent),";leading jet #xi;(1/N_{dijet})dN/d#xi",10,0,5.5);
  
  n_dijet_binned[ialgo][icent]=t_jet[ialgo]->GetEntries(Form("pt1>120 && pt2>30 && dphi>2.094 && cent>=%d && cent<%d",cent_min[icent],cent_max[icent]));
  t[ialgo]->Draw(Form("-log((pt/(pt1*cosh(-eta1)*cosh(-eta1)))*(cos(phi)*cos(phi1)+sin(phi)*sin(phi1)+sinh(-eta1)*sinh(eta)))>>h_binned_ref_%d_%d",ialgo,icent),Form("abs(eta)<2.4 && pt1>120 && pt2>30 && dphi>2.094  && sqrt(pow(acos(cos(phi-phi1)),2)+pow(eta+eta1,2))<0.3 && pt>1 && abs(eta1)>0.3 && cent>=%d && cent<%d",cent_min[icent],cent_max[icent]));
  t[ialgo]->Draw(Form("-log((pt/(pt1*cosh(eta1)*cosh(eta1)))*(cos(phi)*cos(phi1)+sin(phi)*sin(phi1)+sinh(eta1)*sinh(eta)))>>h_binned_%d_%d",ialgo,icent),Form("abs(eta)<2.4 && pt1>120 && pt2>30 && dphi>2.094 && sqrt(pow(acos(cos(phi-phi1)),2)+pow(eta-eta1,2))<0.3 && pt>1 && abs(eta1)>0.3 && cent>=%d && cent<%d",cent_min[icent],cent_max[icent]));
  
   h_binned[ialgo][icent]->Add(h_binned_ref[ialgo][icent],-1);
   h_binned[ialgo][icent]->Scale(1/n_dijet_binned[ialgo][icent]);
   h_binned[ialgo][icent]->Scale(1/h_binned[ialgo][icent]->GetBinWidth(1));
 }
 
 h[ialgo]->Add(h_ref[ialgo],-1);
 cout<<n_dijet[ialgo]<<endl;
 h[ialgo]->Scale(1/n_dijet[ialgo]);
 h[ialgo]->Scale(1/h[ialgo]->GetBinWidth(1));
 h_wrong[ialgo]->Add(h_wrong_ref[ialgo],-1);
 h_wrong[ialgo]->Scale(1/n_dijet[ialgo]);
 h_wrong[ialgo]->Scale(1/h_wrong[ialgo]->GetBinWidth(1));

}

 TLegend *t4=new TLegend(0.65,0.8,0.95,0.95); 
 t4->SetFillColor(0);
 t4->SetBorderSize(0);
 t4->SetFillStyle(0); 
 t4->SetTextFont(63);
 t4->SetTextSize(22);

 TLegend *t2=new TLegend(0.65,0.8,0.95,0.95); 
 t2->SetFillColor(0);
 t2->SetBorderSize(0);
 t2->SetFillStyle(0); 
 t2->SetTextFont(63);
 t2->SetTextSize(22);

 TLegend *t1=new TLegend(0.65,0.8,0.95,0.95); 
 t1->SetFillColor(0);
 t1->SetBorderSize(0);
 t1->SetFillStyle(0); 
 t1->SetTextFont(63);
 t1->SetTextSize(22);

TCanvas *c0 = new TCanvas("c0","",600,600);
c0->SetLogy();
h_ref[0]->SetMarkerStyle(24);
h_ref[0]->SetMarkerColor(kRed);
h_ref[0]->SetLineColor(kRed);
h[0]->SetMarkerColor(kRed);
h[0]->SetLineColor(kRed);
h_ref[1]->SetMarkerStyle(24);
h_ref[0]->SetMaximum(500000);
h_ref[0]->SetMinimum(0.001);
h_ref[0]->Draw();
h_ref[1]->Draw("same");
h[0]->Draw("same");
h[1]->Draw("same");
t2->AddEntry(h[0],"signal","p");
t2->AddEntry(h_ref[0],"background","p");
t2->Draw("same");
 drawText("p_{T,1}>120, p_{T,2}>30 GeV/c",0.25,0.9);
 drawText("p_{T}>1 GeV/c, no track correct.",0.25,0.85);
 drawText("0.3<|#eta_{1}|<2, #Delta#phi>2#pi/3",0.25,0.8);
 drawText("0-100%, PYTHIA+HYDJET",0.25,0.75);
c0->SaveAs("sig_bckg.png");

TCanvas *c1 = new TCanvas("c1","",600,600);
 h[0]->SetMaximum(2.5);

h[0]->Draw();
h[1]->Draw("same");
 t4->AddEntry(h[0],"Pu, PF","p");
 t4->AddEntry(h[1],"Pu, Calo","p");
 t4->Draw("same");
 drawText("p_{T,1}>120, p_{T,2}>30 GeV/c",0.25,0.9);
 drawText("p_{T}>1 GeV/c, no track correct.",0.25,0.85);
 drawText("0.3<|#eta_{1}|<2, #Delta#phi>2#pi/3",0.25,0.8);
 drawText("0-100%, PYTHIA+HYDJET",0.25,0.75);
c1->SaveAs("FF_algo.png");

TCanvas *c2 = new TCanvas("c2","",600,600);
TH1D *hrat= (TH1D*)h[0]->Clone("hrat");
hrat->Divide(h[1]);
hrat->SetMaximum(1.5);
hrat->SetMinimum(0.8);
hrat->GetYaxis()->SetTitle("akPu3PF/akPu3Calo");
hrat->Draw();
 drawText("p_{T,1}>120, p_{T,2}>30 GeV/c",0.25,0.9);
 drawText("p_{T}>1 GeV/c, no track correct.",0.25,0.85);
 drawText("0.3<|#eta_{1}|<2, #Delta#phi>2#pi/3",0.25,0.8);
 drawText("0-100%, PYTHIA+HYDJET",0.25,0.75);
c2->SaveAs("FF_algo_rat.png");
 
TCanvas *c3 = new TCanvas("c3","",600,600);
h_wrong[0]->SetMarkerStyle(25);
h_wrong[0]->SetMaximum(2.5);

h_wrong[0]->Draw();
h[0]->Draw("same");
t1->AddEntry(h[0],"correct","p");
t1->AddEntry(h_wrong[0],"wrong","p");
t1->Draw("same");
 drawText("p_{T,1}>120, p_{T,2}>30 GeV/c",0.25,0.9);
 drawText("p_{T}>1 GeV/c, no track correct.",0.25,0.85);
 drawText("0.3<|#eta_{1}|<2, #Delta#phi>2#pi/3",0.25,0.8);
 drawText("0-100%, PYTHIA+HYDJET",0.25,0.75);
c3->SaveAs("wrong_FF.png"); 

TCanvas *c4 = new TCanvas("c4","",600,600);
TH1D * h_wrong_rat = (TH1D*)h_wrong[0]->Clone("h_wrong_rat");
 h_wrong_rat->Divide(h[0]);
 h_wrong_rat->GetYaxis()->SetTitle("angle proj/dR proj");
 h_wrong_rat->Draw();
 drawText("p_{T,1}>120, p_{T,2}>30 GeV/c",0.25,0.9);
 drawText("p_{T}>1 GeV/c, no track correct.",0.25,0.85);
 drawText("0.3<|#eta_{1}|<2, #Delta#phi>2#pi/3",0.25,0.8);
 drawText("0-100%, PYTHIA+HYDJET",0.25,0.75);
 drawText("akPu3PF",0.25,0.7);
c4->SaveAs("wrong_FF_rat.png"); 

  TCanvas * c5= new TCanvas("c5","",1050,700);
  makeMultiPanelCanvas(c5,3,2,0.0,0.0,0.28,0.2,0.02);
  TH1D *h_binned_rat[ncent];
  for(int icent=0;icent<ncent;icent++){
   h_binned_rat[icent]=(TH1D*)h_binned[0][icent]->Clone(Form("h_binned_rat_%d",icent));
   h_binned_rat[icent]->Divide(h_binned[1][icent]);
   c5->cd(icent+1);
   h_binned_rat[icent]->SetMaximum(1.5);
   h_binned_rat[icent]->SetMinimum(0.8);
   h_binned_rat[icent]->Draw();
   drawText(Form("%d-%d %%",(int)(0.5*cent_min[icent]),(int)(0.5*cent_max[icent])),0.65,0.8);
  }
  c5->SaveAs("centrality_binned.pdf");
  c5->SaveAs("centrality_binned.png");
}