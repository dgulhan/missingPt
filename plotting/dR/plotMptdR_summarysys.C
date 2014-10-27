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
#include "TF1.h"
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
void plot(){
TH1D::SetDefaultSumw2();

TFile * f1 = new TFile("PbPb-pp_systematic_errors_recoreco_gengen.root");
TFile * f2 = new TFile("PbPb-pp_systematic_errors_genreco_gengen.root");
TFile * filejes = new TFile("PbPb-pp_systematic_JES_LeadJetUp.root");

int ncent=2;
int nAj=3;
TH1D *hrecoreco[ncent][nAj];
TH1D *hrecojet[ncent][nAj];
TH1D *hrecotrk[ncent][nAj];
TH1D *hjes[ncent][nAj];
TString scent[]={"peripheral","central"};
TString sAj[]={"Aj01","Aj022","Aj221"};
int centmin[]={0,30};
int centmax[]={30,100};
double Ajmin[]={0.,0.,0.22};
double Ajmax[]={1,0.22,1};
TF1 *frecoreco[ncent][nAj];
TF1 *frecojet[ncent][nAj];
TF1 *frecotrk[ncent][nAj];
TF1 *frecoreco2[ncent][nAj];
TF1 *frecojet2[ncent][nAj];
TF1 *frecotrk2[ncent][nAj];
TF1 *frecoreco3[ncent][nAj];
TF1 *frecojet3[ncent][nAj];
TF1 *frecotrk3[ncent][nAj];
TF1 *fjes[ncent][nAj];
TF1 *fjes2[ncent][nAj];
TF1 *ftrk[ncent][nAj];
TF1 *ftrk2[ncent][nAj];

 TLegend *leg;
 leg= new TLegend(0.5,0.5,0.9,0.8);
 leg->SetBorderSize(0); 
 leg->SetFillStyle(0);
 leg->SetTextFont(43);
 leg->SetTextSize(18);
 
TFile *filedata[nAj];
filedata[0] = new TFile("results_Aj2.root");
filedata[1] = new TFile("results_Aj0.root");
filedata[2] = new TFile("results_Aj1.root");

TH1D *hdata[ncent][nAj];
TH1D *htrk[ncent][nAj];
for(int icent=0;icent<ncent;icent++){
 for(int iAj=0;iAj<nAj;iAj++){
  // frecoreco[icent][iAj]= new TF1(Form("frecoreco%d%d",icent,iAj),"[0]*(x-[1])*(x-[2])*(x-3)",0.3,2);
  // frecojet[icent][iAj]= new TF1(Form("frecoreco%d%d",icent,iAj),"[0]*(x-[1])*(x-[2])*(x-3)",0.3,2);
  // frecotrk[icent][iAj]= new TF1(Form("frecoreco%d%d",icent,iAj),"[0]*(x-[1])*(x-[2])*(x-3)",0.3,2);
  
  if(iAj==2){
   frecoreco2[icent][iAj]= new TF1(Form("frecoreco2%d%d",icent,iAj),"[0]*(x-[1])",0.2,1.2);
   frecojet2[icent][iAj]= new TF1(Form("frecojet2%d%d",icent,iAj),"[0]*(x-[1])",0.2,1.2);
   frecoreco[icent][iAj]= new TF1(Form("frecoreco%d%d",icent,iAj),"[0]*(x-[1])",1.2,1.8);
   frecojet[icent][iAj]= new TF1(Form("frecojet%d%d",icent,iAj),"[0]*(x-[1])",1.2,1.8);
   frecotrk[icent][iAj]= new TF1(Form("frecotrk%d%d",icent,iAj),"[0]*(x-[1])",1.2,1.8);
   frecotrk2[icent][iAj]= new TF1(Form("frecotrk2%d%d",icent,iAj),"[0]*(x-[1])",0.2,1.2);
   fjes2[icent][iAj]= new TF1(Form("fjes2%d%d",icent,iAj),"[0]*(x-[1])",0.2,1.2);
   fjes[icent][iAj]= new TF1(Form("fjes%d%d",icent,iAj),"[0]*(x-[1])",1.2,1.8);
   ftrk2[icent][iAj]= new TF1(Form("ftrk2%d%d",icent,iAj),"[0]*(x-[1])",0.2,1.2);
   ftrk[icent][iAj]= new TF1(Form("ftrk%d%d",icent,iAj),"[0]*(x-[1])",1.2,1.8);
  }
  else{
   frecoreco2[icent][iAj]= new TF1(Form("frecoreco2%d%d",icent,iAj),"[0]*(x-[1])",0.2,0.8);
   frecojet2[icent][iAj]= new TF1(Form("frecojet2%d%d",icent,iAj),"[0]*(x-[1])",0.2,0.8);
   frecoreco[icent][iAj]= new TF1(Form("frecoreco%d%d",icent,iAj),"[0]*(x-[1])",0.8,1.8);
   frecojet[icent][iAj]= new TF1(Form("frecojet%d%d",icent,iAj),"[0]*(x-[1])",0.8,1.8);
   frecotrk[icent][iAj]= new TF1(Form("frecotrk%d%d",icent,iAj),"[0]*(x-[1])",0.8,1.8);
   frecotrk2[icent][iAj]= new TF1(Form("frecotrk2%d%d",icent,iAj),"[0]*(x-[1])",0.2,0.8);
   fjes2[icent][iAj]= new TF1(Form("fjes2%d%d",icent,iAj),"[0]*(x-[1])",0.2,0.8);
   fjes[icent][iAj]= new TF1(Form("fjes%d%d",icent,iAj),"[0]*(x-[1])",0.8,1.8);
   ftrk2[icent][iAj]= new TF1(Form("ftrk2%d%d",icent,iAj),"[0]*(x-[1])",0.2,0.8);
   ftrk[icent][iAj]= new TF1(Form("ftrk%d%d",icent,iAj),"[0]*(x-[1])",0.8,1.8);
  }
  frecoreco3[icent][iAj]= new TF1(Form("frecoreco3%d%d",icent,iAj),"[0]*(x-[1])",0.,0.2);
  frecojet3[icent][iAj]= new TF1(Form("frecojet3%d%d",icent,iAj),"[0]*(x-[1])",0.,0.2);
  frecotrk3[icent][iAj]= new TF1(Form("frecotrk3%d%d",icent,iAj),"[0]*(x-[1])",0.,0.2);
  
  // frecoreco[icent][iAj]->SetParameters(2,0.5,1,1.5);
  // frecojet[icent][iAj]->SetParameters(2,0.5,1,1.5);
  // frecotrk[icent][iAj]->SetParameters(2,0.5,1,1.5);
  
  frecoreco[icent][iAj]->SetParameters(2,-0.5);
  frecojet[icent][iAj]->SetParameters(2,-0.5);
  frecotrk[icent][iAj]->SetParameters(2,-0.5);
  frecoreco2[icent][iAj]->SetParameters(2,-0.5);
  frecojet2[icent][iAj]->SetParameters(2,-0.5);
  frecotrk2[icent][iAj]->SetParameters(2,-0.5);
  frecoreco3[icent][iAj]->SetParameters(2,-0.5);
  frecojet3[icent][iAj]->SetParameters(2,-0.5);
  frecotrk3[icent][iAj]->SetParameters(2,-0.5);
  fjes2[icent][iAj]->SetParameters(2,-0.5);
  fjes[icent][iAj]->SetParameters(2,-0.5);
  ftrk2[icent][iAj]->SetParameters(2,-0.5);
  ftrk[icent][iAj]->SetParameters(2,-0.5);
  
  hrecoreco[icent][iAj]=(TH1D*)f1->Get(Form("%s_%s",scent[icent].Data(),sAj[iAj].Data()));
  hrecojet[icent][iAj]=(TH1D*)f2->Get(Form("%s_%s",scent[icent].Data(),sAj[iAj].Data()));
  hjes[icent][iAj]=(TH1D*)filejes->Get(Form("%s_%s",scent[icent].Data(),sAj[iAj].Data()));
  
  hrecotrk[icent][iAj]=(TH1D*)hrecoreco[icent][iAj]->Clone(Form("hrecotrk%d%d",icent,iAj));
  
  hdata[icent][iAj]=(TH1D*)filedata[iAj]->Get(Form("hmpt_diff_5_%d",icent));
  
  htrk[icent][iAj]=(TH1D*)hdata[icent][iAj]->Clone(Form("htrk%d%d",icent,iAj));
  htrk[icent][iAj]->Scale(0.05); 
 
  hrecotrk[icent][iAj]->Add(hrecojet[icent][iAj],-1);
  hrecotrk[icent][iAj]->SetMarkerColor(kBlue);
  hrecotrk[icent][iAj]->SetLineColor(kBlue);
  hrecojet[icent][iAj]->SetMarkerColor(kBlack);
  hrecojet[icent][iAj]->SetLineColor(kBlack);
  hrecoreco[icent][iAj]->SetMarkerColor(kRed);
  hrecoreco[icent][iAj]->SetLineColor(kRed);
  hjes[icent][iAj]->SetMarkerColor(kGreen+1);
  hjes[icent][iAj]->SetLineColor(kGreen+1);
  htrk[icent][iAj]->SetMarkerColor(kViolet);
  htrk[icent][iAj]->SetLineColor(kViolet);
  for(int ibin=0;ibin<10;ibin++){
   hrecoreco[icent][iAj]->SetBinContent(ibin+1,fabs(hrecoreco[icent][iAj]->GetBinContent(ibin+1)));
   hrecojet[icent][iAj]->SetBinContent(ibin+1,fabs(hrecojet[icent][iAj]->GetBinContent(ibin+1)));
   hrecotrk[icent][iAj]->SetBinContent(ibin+1,fabs(hrecotrk[icent][iAj]->GetBinContent(ibin+1)));
   hjes[icent][iAj]->SetBinContent(ibin+1,fabs(hjes[icent][iAj]->GetBinContent(ibin+1)));
   htrk[icent][iAj]->SetBinContent(ibin+1,fabs(htrk[icent][iAj]->GetBinContent(ibin+1)));
  }
  hrecoreco[icent][iAj]->Fit(frecoreco[icent][iAj],"R");
  hrecojet[icent][iAj]->Fit(frecojet[icent][iAj],"R");
  hrecotrk[icent][iAj]->Fit(frecotrk[icent][iAj],"R");
  frecoreco[icent][iAj]->SetLineColor(kRed);
  frecojet[icent][iAj]->SetLineColor(kBlack);
  frecotrk[icent][iAj]->SetLineColor(kBlue);
  hrecoreco[icent][iAj]->Fit(frecoreco2[icent][iAj],"R");
  hrecojet[icent][iAj]->Fit(frecojet2[icent][iAj],"R");
  hrecotrk[icent][iAj]->Fit(frecotrk2[icent][iAj],"R");
  frecoreco2[icent][iAj]->SetLineColor(kRed);
  frecojet2[icent][iAj]->SetLineColor(kBlack);
  frecotrk2[icent][iAj]->SetLineColor(kBlue);
  hrecoreco[icent][iAj]->Fit(frecoreco3[icent][iAj],"R");
  hrecojet[icent][iAj]->Fit(frecojet3[icent][iAj],"R");
  hrecotrk[icent][iAj]->Fit(frecotrk3[icent][iAj],"R");
  hjes[icent][iAj]->Fit(fjes2[icent][iAj],"R");
  hjes[icent][iAj]->Fit(fjes[icent][iAj],"R");
  htrk[icent][iAj]->Fit(ftrk2[icent][iAj],"R");
  htrk[icent][iAj]->Fit(ftrk[icent][iAj],"R");
  fjes2[icent][iAj]->SetLineColor(kGreen+1);
  fjes[icent][iAj]->SetLineColor(kGreen+1);
  ftrk2[icent][iAj]->SetLineColor(kViolet);
  ftrk[icent][iAj]->SetLineColor(kViolet);
 }
}

leg->AddEntry(hrecotrk[0][0],"trk reco","l");
leg->AddEntry(hrecojet[0][0],"jet reco","l");
leg->AddEntry(hrecoreco[0][0],"trk+jet reco","l");
leg->AddEntry(hjes[0][0],"residual jes","l");
leg->AddEntry(htrk[0][0],"residual trk","l");

TH1D * empty = new TH1D("empty",";#DeltaR;PbPb-pp systematics",20,0.,2);
empty->SetMaximum(3.499);
empty->SetMinimum(0);

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
 empty->GetXaxis()->SetNdivisions(505);
 empty->GetYaxis()->SetNdivisions(505);
 TCanvas *c1 = new TCanvas("c1","",(ncent+1)*300,700);
 makeMultiPanelCanvas(c1,nAj,ncent,0.0,0.0,0.2,0.2,0.02);
 for(int icent=0;icent<ncent;icent++){
  for(int iAj=0;iAj<nAj;iAj++){
   c1->cd(icent*3+iAj+1);
   empty->Draw();
   hrecojet[icent][iAj]->Draw("same hist");
   hrecoreco[icent][iAj]->Draw("same hist");
   hrecotrk[icent][iAj]->Draw("same hist");
   hjes[icent][iAj]->Draw("same hist");
   htrk[icent][iAj]->Draw("same hist");
   frecojet[icent][iAj]->Draw("same hist");
   frecoreco[icent][iAj]->Draw("same hist");
   frecotrk[icent][iAj]->Draw("same hist");
   frecojet2[icent][iAj]->Draw("same hist");
   frecoreco2[icent][iAj]->Draw("same hist");
   frecotrk2[icent][iAj]->Draw("same hist");
   fjes[icent][iAj]->Draw("same hist");
   fjes2[icent][iAj]->Draw("same hist");
   ftrk[icent][iAj]->Draw("same hist");
   ftrk2[icent][iAj]->Draw("same hist");
   drawText(Form("%d - %d %%",centmin[icent],centmax[icent]),0.5,0.9);
   drawText(Form("%.2f < A_{J} < %.2f",Ajmin[iAj],Ajmax[iAj]),0.5,0.8);
   if(icent==0 && iAj==nAj-1) leg->Draw("same");
  }
 }
 c1->SaveAs("systematics/PbPbmPP_systematics.png");
 c1->SaveAs("systematics/PbPbmPP_systematics.pdf");

}