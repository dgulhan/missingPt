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
		 
		 if(i==0 && j==1){
          pad[i][j] = new TPad(padName.Data(),padName.Data(),Xlow[i],Ylow[j],Xup[i]*0.837,Yup[j]*0.935);
		 }
		 else if(i==0 && j==0){
		  pad[i][j] = new TPad(padName.Data(),padName.Data(),Xlow[i],Ylow[j]*0.935,Xup[i],Yup[j]);
		 }
		 else if(i==1 && j==1){
		  pad[i][j] = new TPad(padName.Data(),padName.Data(),0.837*Xlow[i],Ylow[j],Xup[i],Yup[j]);
		 }
         else pad[i][j] = new TPad(padName.Data(),padName.Data(),
            Xlow[i],Ylow[j],Xup[i],Yup[j]);
         if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
		 else if(i==1 && j==1) pad[i][j]->SetLeftMargin(0.57*PadWidth);
         else pad[i][j]->SetLeftMargin(0);

         if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
		   else pad[i][j]->SetRightMargin(0);

         if(j==0){
          if(i==0)pad[i][j]->SetTopMargin(edge);
          else pad[i][j]->SetTopMargin(edge*0.9);
         }else pad[i][j]->SetTopMargin(0);

         if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
		    else if(i==0 && j==0) pad[i][j]->SetBottomMargin(0.17*PadHeight);
         else pad[i][j]->SetBottomMargin(0);

         pad[i][j]->Draw();
         pad[i][j]->cd();
         pad[i][j]->SetNumber(columns*j+i+1);
      }
   }
   pad[0][0]->cd();

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


void cumulative(int gensigbkgreco=3, int doMC=0, int index_var=2, bool doIntegrate=true,  int iAj=0, int doAve=1){
 TH1D::SetDefaultSumw2(); 
 int ncent=2;  
 int npt=6;
 int  cent_min[]={0,60};
 int cent_max[]={60,200};
 double ptmin[]={  8,4,2,1,0.5,0.5};
 double ptmax[]={300,8,4,2,  1,300};

 int col[]={kRed+1,kGreen+3,kOrange+1,kYellow-9,kBlue-9,1};
 
 TH1D * hmpt_vs_alpha[npt][ncent]; //PbPb
 TH1D * hmpt_vs_alpha_diff[npt][ncent]; //PbPb-pp
 TH1D * hmpt_vs_alpha_ref[npt];//pp
    
 TLegend *leg,*leg2, *leg3;
 leg= new TLegend(0.45,0.1,0.99,0.4);
 leg->SetBorderSize(0); 
 leg->SetFillStyle(0);
 leg->SetTextFont(43);
 leg->SetTextSize(18);
 leg3 = new TLegend(0.5,0.1,0.95,0.4);
 leg3->SetBorderSize(0); 
 leg3->SetFillStyle(0); 
 leg3->SetTextFont(43);
 leg3->SetTextSize(18);
 leg2 = new TLegend(0.05,0.1,0.6,0.4); 
 leg2->SetBorderSize(0); 
 leg2->SetFillStyle(0); 
 leg2->SetTextFont(43);
 leg2->SetTextSize(18);
 
 for(int ipt=0;ipt<npt;ipt++){ 
  if(ipt<npt-1){
   hmpt_vs_alpha_ref[ipt]->SetLineColor(1);
   hmpt_vs_alpha_ref[ipt]->SetMarkerColor(1);
   hmpt_vs_alpha_ref[ipt]->SetMarkerSize(0);
   hmpt_vs_alpha_ref[ipt]->SetFillColor(col[ipt]);
  }else{ 
     hmpt_vs_alpha_ref[ipt]->SetMarkerStyle(25);
  }
  for(int icent=0;icent<ncent;icent++){
   if(ipt<npt-1){
	hmpt_vs_alpha[ipt][icent]->SetLineColor(1);
    hmpt_vs_alpha[ipt][icent]->SetMarkerColor(1);
    hmpt_vs_alpha[ipt][icent]->SetMarkerSize(0);
    hmpt_vs_alpha[ipt][icent]->SetFillColor(col[ipt]); 
	if(icent==0){
     if(ipt>2) leg2->AddEntry(hmpt_vs_alpha[ipt][icent],Form("%.1f - %.1f",ptmin[ipt],ptmax[ipt]),"f");
     else leg->AddEntry(hmpt_vs_alpha[ipt][icent],Form("%.1f - %.1f",ptmin[ipt],ptmax[ipt]),"f");
    }
   }else{
     hmpt_vs_alpha[ipt][icent]->SetMarkerStyle(5);
	 if(icent==0)leg2->AddEntry(hmpt_vs_alpha[npt-1][icent],">0.5","p");
   }   
    
   if(ipt<npt-1){
    hmpt_vs_alpha_diff[ipt][icent]->SetLineColor(1);
    hmpt_vs_alpha_diff[ipt][icent]->SetMarkerColor(1);
    hmpt_vs_alpha_diff[ipt][icent]->SetMarkerSize(0);
    hmpt_vs_alpha_diff[ipt][icent]->SetFillColor(col[ipt]);
   }else{ 
     hmpt_vs_alpha_diff[ipt][icent]->SetMarkerStyle(24);
	 if(icent==0)leg2->AddEntry(hmpt_vs_alpha[npt-1][icent],">0.5","p");

   }
  }
 }
 
 TFile *syst_PbPb = new TFile(Form("systematics/outf_Aj%d.root",iAj));
 TFile *syst_pp = new TFile(Form("systematics_pp/outf_Aj%d.root",iAj));

 TH1D * hsyst[ncent];
 TH1D * hsyst_PbPb_pp[ncent];
 TH1D * hsyst_pp=(TH1D *)syst_pp->Get(Form("hmpt_diff_%d_%d",npt-1,0));
 
 for(int icent=0;icent<ncent;icent++){
  hsyst[icent]=(TH1D *)syst_PbPb->Get(Form("hmpt_diff_%d_%d",npt-1,icent));
 }
 
 TLine* zeroLine_p = new TLine(0., 0., 2., 0.);
 zeroLine_p->SetLineColor(1);
 zeroLine_p->SetLineStyle(2);
  
 TH1D * empty=new TH1D("empty",Form(";%s;<#slash{p}_{T}^{#parallel}> (GeV/c)",axistitle[index_var].Data()),nalpha/2+1,frac);
 empty->Fill(0.5,1000); 
 empty->SetMaximum(20); 
 empty->SetMinimum(-40);  
 empty->GetXaxis()->CenterTitle();
 empty->GetYaxis()->CenterTitle();
 empty->GetXaxis()->SetNdivisions(505);
 empty->GetXaxis()->SetNdivisions(505);
 empty->GetYaxis()->SetNdivisions(505);
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
  
 TH1D * empty2=new TH1D("empty2",Form(";%s;<#slash{p}_{T}^{#parallel}> (GeV/c)",axistitle[index_var].Data()),nalpha/2+1,frac);
 empty2->Fill(0.5,1000); 
 empty2->SetMaximum(10); 
 empty2->SetMinimum(-10);
 empty2->GetXaxis()->CenterTitle();
 empty2->GetYaxis()->CenterTitle();
 empty2->GetXaxis()->SetNdivisions(505);
 empty2->GetYaxis()->SetNdivisions(505);
 empty2->GetXaxis()->SetTitleSize(28);
 empty2->GetXaxis()->SetTitleFont(43); 
 empty2->GetXaxis()->SetTitleOffset(2.2);
 empty2->GetXaxis()->SetLabelSize(22);
 empty2->GetXaxis()->SetLabelFont(43);
 empty2->GetYaxis()->SetTitleSize(28);
 empty2->GetYaxis()->SetTitleFont(43); 
 empty2->GetYaxis()->SetTitleOffset(1.8);
 empty2->GetYaxis()->SetLabelSize(22);
 empty2->GetYaxis()->SetLabelFont(43);
   
   
 TCanvas *c1 = new TCanvas("c1","",(ncent+1)*300,700);
 makeMultiPanelCanvas(c1,ncent+1,2,0.0,0.0,0.2,0.2,0.02);
 
  //###############panel with text and legend####################################

 c1->cd(ncent+2); 
 leg2->Draw("same");
 leg->Draw("same");
 if(Ajmin[iAj]==0 && Ajmax[iAj]==1) drawText(Form("|#eta_{1}|,|#eta_{2}|<%.2f",etadijet),0.05,0.54);
 else if(Ajmin[iAj]==0) drawText(Form("|#eta_{1}|,|#eta_{2}|<%.2f, A_{J} < %.2f",etadijet,Ajmax[iAj]),0.05,0.54);
 else if(Ajmax[iAj]==1)drawText(Form("|#eta_{1}|,|#eta_{2}|<%.2f,  A_{J} > %.2f",etadijet,Ajmin[iAj]),0.05,0.54);
 else drawText(Form("|#eta_{1}|,|#eta_{2}|<%.2f, %.2f < A_{J} < %.0f",etadijet,Ajmin[iAj],Ajmax[iAj]),0.05,0.54);
 drawText(Form("p_{T,1}>%d, p_{T,2}>%d GeV/c",jetpt1,jetpt2),0.05,0.9);
 drawText(Form("|#eta_{jet}|<2, #Delta#phi_{1,2}>5#pi/6"),0.05,0.78);
 drawText(Form("anti-kt Vs Calo R=0.3"),0.05,0.66);
 drawText(Form("|#eta_{trk}|<2.4"),0.05,0.42);
 
 //###############PbPb-pp panels####################################
 for(int icent=0;icent<ncent;icent++){ 
  c1->cd((2*(ncent+1))-icent);
  else empty2->Draw();
  for(int ipt=npt-2;ipt>=0;ipt--){ 
   hmpt_vs_alpha_diff[ipt][icent]->Draw("same");
   hmpt_vs_alpha_diff[ipt][icent]->Draw("same hist");
  }
 
  for(int ialpha = 0 ;ialpha<hmpt_vs_alpha_diff[npt-1][icent]->GetNbinsX();ialpha++){
   double  y = hmpt_vs_alpha_diff[npt-1][icent]->GetBinContent(ialpha+1);
   double yerr;
   if(ialpha==0) yerr = 3;
   else  yerr = systematics[icent];
   double x=hmpt_vs_alpha_diff[npt-1][icent]->GetBinCenter(ialpha+1);
  
   TLine * l1= new TLine(x-0.035,y-yerr,x+0.035,y-yerr);
   TLine * l2= new TLine(x-0.035,y-yerr,x+0.035,y-yerr);
   TLine * l3= new TLine(x-0.035,y+yerr,x+0.035,y+yerr);
   TLine * l4= new TLine(x-0.035,y+yerr,x+0.035,y+yerr);
   TLine * l5= new TLine(x-0.035,y-yerr,x-0.035,y-yerr+0.25);
   TLine * l6= new TLine(x-0.035,y+yerr,x-0.035,y+yerr-0.25);
   TLine * l7= new TLine(x+0.035,y-yerr,x+0.035,y-yerr+0.25);
   TLine * l8= new TLine(x+0.035,y+yerr,x+0.035,y+yerr-0.25);
  
   l1->Draw("same");
   l2->Draw("same");
   l3->Draw("same");
   l4->Draw("same");
   l5->Draw("same");
   l6->Draw("same");
   l7->Draw("same");
   l8->Draw("same");
  }

  hmpt_vs_alpha_diff[npt-1][icent]->Draw("same");

  zeroLine_p->Draw("same");
  c1->cd((2*(ncent+1))-icent)->RedrawAxis();

  if(icent==1)drawText(Form("%d-%d %%",(int)(0.5*cent_min[icent]),(int)(0.5*cent_max[icent])),0.2,0.85);
  else drawText(Form("%d-%d %%",(int)(0.5*cent_min[icent]),(int)(0.5*cent_max[icent])),0.05,0.85);
 
  if(icent==1)  drawText("PbPb-pp",0.2,0.9);
  else  drawText("PbPb-pp",0.05,0.9);
  
  if(icent!=ncent-1)drawPatch(0.0,0.0,0.05,0.19);
  if(icent!=0)drawPatch(0.95,0.0,1,0.19);
  if(icent!=0)drawPatch(0.0,0.95,0.16,1);
 }

 //###############pp panel####################################
 c1->cd(1);
 empty->Draw();
 
 leg3->AddEntry(hmpt_vs_alpha_ref[npt-1],"pp","p");
 leg3->AddEntry(hmpt_vs_alpha[npt-1][0],"PbPb","p");
 leg3->AddEntry(hmpt_vs_alpha_diff[npt-1][0],"PbPb-pp","p");

 leg3->Draw("same");
 
 for(int ipt=npt-2;ipt>=0;ipt--){
  hmpt_vs_alpha_ref[ipt]->Draw("same");
  hmpt_vs_alpha_ref[ipt]->Draw("same hist");
 }
 
 hmpt_vs_alpha_ref[npt-1]->Draw("same");

 for(int ialpha = 0 ;ialpha<hmpt_vs_alpha_ref[npt-1]->GetNbinsX();ialpha++){
  double  y = hmpt_vs_alpha_ref[npt-1]->GetBinContent(ialpha+1);
  double yerr = fabs(hsyst_pp->GetBinContent(ialpha+1));
  double x=hmpt_vs_alpha_ref[npt-1]->GetBinCenter(ialpha+1);
  
   TLine * l1= new TLine(x-0.035,y-yerr,x+0.035,y-yerr);
   TLine * l2= new TLine(x-0.035,y-yerr,x+0.035,y-yerr);
   TLine * l3= new TLine(x-0.035,y+yerr,x+0.035,y+yerr);
   TLine * l4= new TLine(x-0.035,y+yerr,x+0.035,y+yerr);
   TLine * l5= new TLine(x-0.035,y-yerr,x-0.035,y-yerr+0.25);
   TLine * l6= new TLine(x-0.035,y+yerr,x-0.035,y+yerr-0.25);
   TLine * l7= new TLine(x+0.035,y-yerr,x+0.035,y-yerr+0.25);
   TLine * l8= new TLine(x+0.035,y+yerr,x+0.035,y+yerr-0.25);
 
   l1->Draw("same");
   l2->Draw("same");
   l3->Draw("same");
   l4->Draw("same");
   l5->Draw("same");
   l6->Draw("same");
   l7->Draw("same");
   l8->Draw("same");
  }
  
 
 zeroLine_p->Draw("same");
 c1->cd(1)->RedrawAxis();
 
 drawText("pp #sqrt{s_{NN}}=2.76 TeV 5.3 pb^{-1}",0.3,0.9);
 drawText("CMS Preliminary",0.5,0.84);
 
 //###############PbPb panels################################### 

 for(int icent=0;icent<ncent;icent++){
 c1->cd(ncent+1-icent);
 empty->Draw();
  for(int ipt=npt-2;ipt>=0;ipt--){ 
   hmpt_vs_alpha[ipt][icent]->Draw("same");
   hmpt_vs_alpha[ipt][icent]->Draw("same hist");
  }
  hmpt_vs_alpha[npt-1][icent]->Draw("same");

  for(int ialpha = 0 ;ialpha<hmpt_vs_alpha[npt-1][icent]->GetNbinsX();ialpha++){
   double  y = hmpt_vs_alpha[npt-1][icent]->GetBinContent(ialpha+1);
   double yerr = fabs(hsyst[icent]->GetBinContent(ialpha+1));
   double x=hmpt_vs_alpha[npt-1][icent]->GetBinCenter(ialpha+1);
  
   TLine * l1= new TLine(x-0.035,y-yerr,x+0.035,y-yerr);
   TLine * l2= new TLine(x-0.035,y-yerr,x+0.035,y-yerr);
   TLine * l3= new TLine(x-0.035,y+yerr,x+0.035,y+yerr);
   TLine * l4= new TLine(x-0.035,y+yerr,x+0.035,y+yerr);
   TLine * l5= new TLine(x-0.035,y-yerr,x-0.035,y-yerr+TMath::Min(0.25,yerr));
   TLine * l6= new TLine(x-0.035,y+yerr,x-0.035,y+yerr-TMath::Min(0.25,yerr));
   TLine * l7= new TLine(x+0.035,y-yerr,x+0.035,y-yerr+TMath::Min(0.25,yerr));
   TLine * l8= new TLine(x+0.035,y+yerr,x+0.035,y+yerr-TMath::Min(0.25,yerr));
   
   l1->Draw("same");
   l2->Draw("same");
   l3->Draw("same");
   l4->Draw("same");
   l5->Draw("same");
   l6->Draw("same");
   l7->Draw("same");
   l8->Draw("same");
  }
   
  drawText(Form("%d-%d %%",(int)(0.5*cent_min[icent]),(int)(0.5*cent_max[icent])),0.05,0.84);
 
  zeroLine_p->Draw("same");
  c1->cd(ncent+1-icent)->RedrawAxis();
  drawText("PbPb #sqrt{s_{NN}}=2.76 TeV 150 #mub^{-1}",0.05,0.9);
  if(icent==1)drawText(" ",0.4,0.82);

 } 
  
 c1->SaveAs(Form("cumulative_eta%d_doMC%d_%s_doIntegrate%d_doGenJet%d_kind%d_Aj%d_%d_ncent%d.png",(int)(etadijet*10),doMC,svariable[index_var].Data(),doIntegrate,doGenJet,gensigbkgreco,(int)(Ajmin[iAj]*10),(int)(Ajmax[iAj]*10),ncent));
 c1->SaveAs(Form("cumulative_eta%d_doMC%d_%s_doIntegrate%d_doGenJet%d_kind%d_Aj%d_%d_ncent%d.pdf",(int)(etadijet*10),doMC,svariable[index_var].Data(),doIntegrate,doGenJet,gensigbkgreco,(int)(Ajmin[iAj]*10),(int)(Ajmax[iAj]*10),ncent));
 
}