//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Mar 15 17:19:57 2014 by ROOT version 5.32/00
// from TTree HltTree/
// found on file: root://eoscms//eos/cms/store/caf/user/dgulhan/PYTHIA/prod22_ppTracking/pt80_pp2013_P01_prod22_v81_merged_forest_0.root
//////////////////////////////////////////////////////////

#ifndef skimTree_pp_h
#define skimTree_pp_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class skimTree_pp {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           HLT_PAZeroBiasPixel_SingleTrack_v1;
   Int_t           reco_extra;
   Int_t           reco_extra_jet;
   Int_t           gen_step;
   Int_t           pat_step;
   Int_t           ana_step;
   Int_t           pcollisionEventSelection;
   Int_t           pHBHENoiseFilter;
   Int_t           phiEcalRecHitSpikeFilter;
   Int_t           pPAcollisionEventSelectionPA;
   Int_t           phfPosFilter3;
   Int_t           phfNegFilter3;
   Int_t           phfPosFilter2;
   Int_t           phfNegFilter2;
   Int_t           phfPosFilter1;
   Int_t           phfNegFilter1;
   Int_t           phltPixelClusterShapeFilter;
   Int_t           pprimaryvertexFilter;
   Int_t           pBeamScrapingFilter;
   Int_t           pVertexFilterCutG;
   Int_t           pVertexFilterCutGloose;
   Int_t           pVertexFilterCutGtight;
   Int_t           pVertexFilterCutGplus;
   Int_t           pVertexFilterCutE;
   Int_t           pVertexFilterCutEandG;

   // List of branches
   TBranch        *b_HLT_PAZeroBiasPixel_SingleTrack_v1;   //!
   TBranch        *b_reco_extra;   //!
   TBranch        *b_reco_extra_jet;   //!
   TBranch        *b_gen_step;   //!
   TBranch        *b_pat_step;   //!
   TBranch        *b_ana_step;   //!
   TBranch        *b_pcollisionEventSelection;   //!
   TBranch        *b_pHBHENoiseFilter;   //!
   TBranch        *b_phiEcalRecHitSpikeFilter;   //!
   TBranch        *b_pPAcollisionEventSelectionPA;   //!
   TBranch        *b_phfPosFilter3;   //!
   TBranch        *b_phfNegFilter3;   //!
   TBranch        *b_phfPosFilter2;   //!
   TBranch        *b_phfNegFilter2;   //!
   TBranch        *b_phfPosFilter1;   //!
   TBranch        *b_phfNegFilter1;   //!
   TBranch        *b_phltPixelClusterShapeFilter;   //!
   TBranch        *b_pprimaryvertexFilter;   //!
   TBranch        *b_pBeamScrapingFilter;   //!
   TBranch        *b_pVertexFilterCutG;   //!
   TBranch        *b_pVertexFilterCutGloose;   //!
   TBranch        *b_pVertexFilterCutGtight;   //!
   TBranch        *b_pVertexFilterCutGplus;   //!
   TBranch        *b_pVertexFilterCutE;   //!
   TBranch        *b_pVertexFilterCutEandG;   //!

   skimTree_pp(TString infile,TTree *tree=0);
   virtual ~skimTree_pp();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     Close();
   TFile *f;
};

#endif

#ifdef skimTree_pp_cxx
skimTree_pp::skimTree_pp(TString infile,TTree *tree) 
{
    f = TFile::Open(infile.Data());
    tree=(TTree*)f->Get("skimanalysis/HltTree");
	  Init(tree);
	  // f->Close();
}


skimTree_pp::~skimTree_pp()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t skimTree_pp::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t skimTree_pp::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void skimTree_pp::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1", &HLT_PAZeroBiasPixel_SingleTrack_v1, &b_HLT_PAZeroBiasPixel_SingleTrack_v1);
   fChain->SetBranchAddress("reco_extra", &reco_extra, &b_reco_extra);
   fChain->SetBranchAddress("reco_extra_jet", &reco_extra_jet, &b_reco_extra_jet);
   fChain->SetBranchAddress("gen_step", &gen_step, &b_gen_step);
   fChain->SetBranchAddress("pat_step", &pat_step, &b_pat_step);
   fChain->SetBranchAddress("ana_step", &ana_step, &b_ana_step);
   fChain->SetBranchAddress("pcollisionEventSelection", &pcollisionEventSelection, &b_pcollisionEventSelection);
   fChain->SetBranchAddress("pHBHENoiseFilter", &pHBHENoiseFilter, &b_pHBHENoiseFilter);
   fChain->SetBranchAddress("phiEcalRecHitSpikeFilter", &phiEcalRecHitSpikeFilter, &b_phiEcalRecHitSpikeFilter);
   fChain->SetBranchAddress("pPAcollisionEventSelectionPA", &pPAcollisionEventSelectionPA, &b_pPAcollisionEventSelectionPA);
   fChain->SetBranchAddress("phfPosFilter3", &phfPosFilter3, &b_phfPosFilter3);
   fChain->SetBranchAddress("phfNegFilter3", &phfNegFilter3, &b_phfNegFilter3);
   fChain->SetBranchAddress("phfPosFilter2", &phfPosFilter2, &b_phfPosFilter2);
   fChain->SetBranchAddress("phfNegFilter2", &phfNegFilter2, &b_phfNegFilter2);
   fChain->SetBranchAddress("phfPosFilter1", &phfPosFilter1, &b_phfPosFilter1);
   fChain->SetBranchAddress("phfNegFilter1", &phfNegFilter1, &b_phfNegFilter1);
   fChain->SetBranchAddress("phltPixelClusterShapeFilter", &phltPixelClusterShapeFilter, &b_phltPixelClusterShapeFilter);
   fChain->SetBranchAddress("pprimaryvertexFilter", &pprimaryvertexFilter, &b_pprimaryvertexFilter);
   fChain->SetBranchAddress("pBeamScrapingFilter", &pBeamScrapingFilter, &b_pBeamScrapingFilter);
   fChain->SetBranchAddress("pVertexFilterCutG", &pVertexFilterCutG, &b_pVertexFilterCutG);
   fChain->SetBranchAddress("pVertexFilterCutGloose", &pVertexFilterCutGloose, &b_pVertexFilterCutGloose);
   fChain->SetBranchAddress("pVertexFilterCutGtight", &pVertexFilterCutGtight, &b_pVertexFilterCutGtight);
   fChain->SetBranchAddress("pVertexFilterCutGplus", &pVertexFilterCutGplus, &b_pVertexFilterCutGplus);
   fChain->SetBranchAddress("pVertexFilterCutE", &pVertexFilterCutE, &b_pVertexFilterCutE);
   fChain->SetBranchAddress("pVertexFilterCutEandG", &pVertexFilterCutEandG, &b_pVertexFilterCutEandG);
   Notify();
}

Bool_t skimTree_pp::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void skimTree_pp::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t skimTree_pp::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void skimTree_pp::Close(){
  f->Close();
}
#endif // #ifdef skimTree_pp_cxx
