//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar 12 21:23:22 2014 by ROOT version 5.32/00
// from TTree trackTree/v1
// found on file: root://eoscms//eos/cms/store/caf/user/dgulhan/PYTHIA/prod22_ppTracking/pt80_pp2013_P01_prod22_v81_merged_forest_0.root
//////////////////////////////////////////////////////////

#ifndef ppTrack_h
#define ppTrack_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "HltTree_pp.C"
#include "skimTree_pp.C"
#include "HiTree.C"
#include "t.C"
#include "genPart.C"
#include "hi.C"

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class ppTrack {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           nEv;
   Int_t           nLumi;
   Int_t           nBX;
   Int_t           nRun;
   Int_t           N;
   Int_t           nv;
   Float_t         vx[2];   //[nv]
   Float_t         vy[2];   //[nv]
   Float_t         vz[2];   //[nv]
   Float_t         vxErr[2];   //[nv]
   Float_t         vyErr[2];   //[nv]
   Float_t         vzErr[2];   //[nv]
   Int_t           nDaugher[2];   //[nv]
   Int_t           nVtx;
   Int_t           maxVtx;
   Int_t           nTrkVtx[5];   //[nVtx]
   Float_t         normChi2Vtx[5];   //[nVtx]
   Float_t         sumPtVtx[5];   //[nVtx]
   Float_t         xVtx[5];   //[nVtx]
   Float_t         yVtx[5];   //[nVtx]
   Float_t         zVtx[5];   //[nVtx]
   Float_t         xVtxErr[5];   //[nVtx]
   Float_t         yVtxErr[5];   //[nVtx]
   Float_t         zVtxErr[5];   //[nVtx]
   Float_t         vtxDist2D[5];   //[nVtx]
   Float_t         vtxDist2DErr[5];   //[nVtx]
   Float_t         vtxDist2DSig[5];   //[nVtx]
    Float_t         vtxDist3DErr[5];   //[nVtx]
   Float_t         vtxDist3DSig[5];   //[nVtx]
   Int_t           nTrk;
   Float_t         trkPt[1000];   //[nTrk]
   Float_t         trkPtError[1000];   //[nTrk]
   Int_t           trkNHit[1000];   //[nTrk]
   Int_t           trkNlayer[1000];   //[nTrk]
   Float_t         trkEta[1000];   //[nTrk]
   Float_t         trkPhi[1000];   //[nTrk]
   Int_t           trkCharge[1000];   //[nTrk]
   Int_t           trkVtxIndex[1000];   //[nTrk]
   Bool_t          highPurity[1000];   //[nTrk]
   Bool_t          highPuritySetWithPV[1000];   //[nTrk]
   Float_t         trkChi2[1000];   //[nTrk]
   Float_t         trkNdof[1000];   //[nTrk]
   Float_t         trkDxy1[1000];   //[nTrk]
   Float_t         trkDxyError1[1000];   //[nTrk]
   Float_t         trkDz1[1000];   //[nTrk]
   Float_t         trkDzError1[1000];   //[nTrk]
   Bool_t          trkFake[1000];   //[nTrk]
   Float_t         trkAlgo[1000];   //[nTrk]
   Int_t           pfType[1000];   //[nTrk]
   Float_t         pfCandPt[1000];   //[nTrk]
   Float_t         pfSumEcal[1000];   //[nTrk]
   Float_t         pfSumHcal[1000];   //[nTrk]
   Float_t         trkStatus[1000];   //[nTrk]
   Int_t           nParticle;
   Float_t         pStatus[1000];   //[nParticle]
   Float_t         pPId[162];   //[nParticle]
   Float_t         pEta[162];   //[nParticle]
   Float_t         pPhi[162];   //[nParticle]
   Float_t         pPt[162];   //[nParticle]
   Float_t         pAcc[162];   //[nParticle]
   Float_t         pAccPair[162];   //[nParticle]
   Float_t         pNRec[162];   //[nParticle]
   Int_t           pNHit[162];   //[nParticle]
   Float_t         mtrkPt[162];   //[nParticle]
   Float_t         mtrkPtError[162];   //[nParticle]
   Int_t           mtrkNHit[162];   //[nParticle]
   Int_t           mtrkNlayer[162];   //[nParticle]
   Int_t           mtrkNlayer3D[162];   //[nParticle]
   Int_t           mtrkQual[162];   //[nParticle]
   Float_t         mtrkChi2[162];   //[nParticle]
   Float_t         mtrkNdof[162];   //[nParticle]
   Float_t         mtrkDz1[162];   //[nParticle]
   Float_t         mtrkDzError1[162];   //[nParticle]
   Float_t         mtrkDxy1[162];   //[nParticle]
   Float_t         mtrkDxyError1[162];   //[nParticle]
   Float_t         mtrkAlgo[162];   //[nParticle]
   Int_t           mtrkPfType[162];   //[nParticle]
   Float_t         mtrkPfCandPt[162];   //[nParticle]
   Float_t         mtrkPfSumEcal[162];   //[nParticle]
   Float_t         mtrkPfSumHcal[162];   //[nParticle]

   // List of branches
   TBranch        *b_nEv;   //!
   TBranch        *b_nLumi;   //!
   TBranch        *b_nBX;   //!
   TBranch        *b_nRun;   //!
   TBranch        *b_N;   //!
   TBranch        *b_nv;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_vxErr;   //!
   TBranch        *b_vyErr;   //!
   TBranch        *b_vzErr;   //!
   TBranch        *b_nDaugher;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_maxVtx;   //!
   TBranch        *b_nTrkVtx;   //!
   TBranch        *b_normChi2Vtx;   //!
   TBranch        *b_sumPtVtx;   //!
   TBranch        *b_xVtx;   //!
   TBranch        *b_yVtx;   //!
   TBranch        *b_zVtx;   //!
   TBranch        *b_xVtxErr;   //!
   TBranch        *b_yVtxErr;   //!
   TBranch        *b_zVtxErr;   //!
   TBranch        *b_vtxDist2D;   //!
   TBranch        *b_vtxDist2DErr;   //!
   TBranch        *b_vtxDist2DSig;   //!
   TBranch        *b_vtxDist3DErr;   //!
   TBranch        *b_vtxDist3DSig;   //!
   TBranch        *b_nTrk;   //!
   TBranch        *b_trkPt;   //!
   TBranch        *b_trkPtError;   //!
   TBranch        *b_trkNHit;   //!
   TBranch        *b_trkNlayer;   //!
   TBranch        *b_trkEta;   //!
   TBranch        *b_trkPhi;   //!
   TBranch        *b_trkCharge;   //!
   TBranch        *b_trkVtxIndex;   //!
   TBranch        *b_highPurity;   //!
   TBranch        *b_highPuritySetWithPV;   //!
   TBranch        *b_trkChi2;   //!
   TBranch        *b_trkNdof;   //!
   TBranch        *b_trkDxy1;   //!
   TBranch        *b_trkDxyError1;   //!
   TBranch        *b_trkDz1;   //!
   TBranch        *b_trkDzError1;   //!
   TBranch        *b_trkFake;   //!
   TBranch        *b_trkAlgo;   //!
   TBranch        *b_pfType;   //!
   TBranch        *b_pfCandPt;   //!
   TBranch        *b_pfSumEcal;   //!
   TBranch        *b_pfSumHcal;   //!
   TBranch        *b_trkStatus;   //!
   TBranch        *b_nParticle;   //!
   TBranch        *b_pStatus;   //!
   TBranch        *b_pPId;   //!
   TBranch        *b_pEta;   //!
   TBranch        *b_pPhi;   //!
   TBranch        *b_pPt;   //!
   TBranch        *b_pAcc;   //!
   TBranch        *b_pAccPair;   //!
   TBranch        *b_pNRec;   //!
   TBranch        *b_pNHit;   //!
   TBranch        *b_mtrkPt;   //!
   TBranch        *b_mtrkPtError;   //!
   TBranch        *b_mtrkNHit;   //!
   TBranch        *b_mtrkNlayer;   //!
   TBranch        *b_mtrkNlayer3D;   //!
   TBranch        *b_mtrkQual;   //!
   TBranch        *b_mtrkChi2;   //!
   TBranch        *b_mtrkNdof;   //!
   TBranch        *b_mtrkDz1;   //!
   TBranch        *b_mtrkDzError1;   //!
   TBranch        *b_mtrkDxy1;   //!
   TBranch        *b_mtrkDxyError1;   //!
   TBranch        *b_mtrkAlgo;   //!
   TBranch        *b_mtrkPfType;   //!
   TBranch        *b_mtrkPfCandPt;   //!
   TBranch        *b_mtrkPfSumEcal;   //!
   TBranch        *b_mtrkPfSumHcal;   //!

   ppTrack(TString infile,TTree *tree=0);
   virtual ~ppTrack();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Int_t    GetEntriesFast();
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
      virtual void             Close();
	  TFile *f;

};

#endif

#ifdef ppTrack_cxx
ppTrack::ppTrack(TString infile,TTree *tree)
{
   f = TFile::Open(infile);
   cout<<infile.Data()<<endl;
   tree = (TTree*) f->Get("ppTrack/trackTree");
   // fhi = new HiTree(infile);
   // fhlt = new HltTree(infile);
   // fskim = new skimTree(infile);
   // fjet = new t(infile);
   // fgen = new genPart(infile);
   Init(tree);
}

ppTrack::~ppTrack()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ppTrack::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Int_t ppTrack::GetEntriesFast()
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntriesFast();
}

Long64_t ppTrack::LoadTree(Long64_t entry)
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

void ppTrack::Init(TTree *tree)
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

   fChain->SetBranchAddress("nEv", &nEv, &b_nEv);
   fChain->SetBranchAddress("nLumi", &nLumi, &b_nLumi);
   fChain->SetBranchAddress("nBX", &nBX, &b_nBX);
   fChain->SetBranchAddress("nRun", &nRun, &b_nRun);
   fChain->SetBranchAddress("N", &N, &b_N);
   fChain->SetBranchAddress("nv", &nv, &b_nv);
   fChain->SetBranchAddress("vx", vx, &b_vx);
   fChain->SetBranchAddress("vy", vy, &b_vy);
   fChain->SetBranchAddress("vz", vz, &b_vz);
   fChain->SetBranchAddress("vxErr", vxErr, &b_vxErr);
   fChain->SetBranchAddress("vyErr", vyErr, &b_vyErr);
   fChain->SetBranchAddress("vzErr", vzErr, &b_vzErr);
   fChain->SetBranchAddress("nDaugher", nDaugher, &b_nDaugher);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("maxVtx", &maxVtx, &b_maxVtx);
   fChain->SetBranchAddress("nTrkVtx", nTrkVtx, &b_nTrkVtx);
   fChain->SetBranchAddress("normChi2Vtx", normChi2Vtx, &b_normChi2Vtx);
   fChain->SetBranchAddress("sumPtVtx", sumPtVtx, &b_sumPtVtx);
   fChain->SetBranchAddress("xVtx", xVtx, &b_xVtx);
   fChain->SetBranchAddress("yVtx", yVtx, &b_yVtx);
   fChain->SetBranchAddress("zVtx", zVtx, &b_zVtx);
   fChain->SetBranchAddress("xVtxErr", xVtxErr, &b_xVtxErr);
   fChain->SetBranchAddress("yVtxErr", yVtxErr, &b_yVtxErr);
   fChain->SetBranchAddress("zVtxErr", zVtxErr, &b_zVtxErr);
   fChain->SetBranchAddress("vtxDist2D", vtxDist2D, &b_vtxDist2D);
   fChain->SetBranchAddress("vtxDist2DErr", vtxDist2DErr, &b_vtxDist2DErr);
   fChain->SetBranchAddress("vtxDist2DSig", vtxDist2DSig, &b_vtxDist2DSig);
   fChain->SetBranchAddress("vtxDist3DErr", vtxDist3DErr, &b_vtxDist3DErr);
   fChain->SetBranchAddress("vtxDist3DSig", vtxDist3DSig, &b_vtxDist3DSig);
   fChain->SetBranchAddress("nTrk", &nTrk, &b_nTrk);
   fChain->SetBranchAddress("trkPt", trkPt, &b_trkPt);
   fChain->SetBranchAddress("trkPtError", trkPtError, &b_trkPtError);
   fChain->SetBranchAddress("trkNHit", trkNHit, &b_trkNHit);
   fChain->SetBranchAddress("trkNlayer", trkNlayer, &b_trkNlayer);
   fChain->SetBranchAddress("trkEta", trkEta, &b_trkEta);
   fChain->SetBranchAddress("trkPhi", trkPhi, &b_trkPhi);
   fChain->SetBranchAddress("trkCharge", trkCharge, &b_trkCharge);
   fChain->SetBranchAddress("trkVtxIndex", trkVtxIndex, &b_trkVtxIndex);
   fChain->SetBranchAddress("highPurity", highPurity, &b_highPurity);
   fChain->SetBranchAddress("highPuritySetWithPV", highPuritySetWithPV, &b_highPuritySetWithPV);
   fChain->SetBranchAddress("trkChi2", trkChi2, &b_trkChi2);
   fChain->SetBranchAddress("trkNdof", trkNdof, &b_trkNdof);
   fChain->SetBranchAddress("trkDxy1", trkDxy1, &b_trkDxy1);
   fChain->SetBranchAddress("trkDxyError1", trkDxyError1, &b_trkDxyError1);
   fChain->SetBranchAddress("trkDz1", trkDz1, &b_trkDz1);
   fChain->SetBranchAddress("trkDzError1", trkDzError1, &b_trkDzError1);
   fChain->SetBranchAddress("trkFake", trkFake, &b_trkFake);
   fChain->SetBranchAddress("trkAlgo", trkAlgo, &b_trkAlgo);
   fChain->SetBranchAddress("pfType", pfType, &b_pfType);
   fChain->SetBranchAddress("pfCandPt", pfCandPt, &b_pfCandPt);
   fChain->SetBranchAddress("pfSumEcal", pfSumEcal, &b_pfSumEcal);
   fChain->SetBranchAddress("pfSumHcal", pfSumHcal, &b_pfSumHcal);
   fChain->SetBranchAddress("trkStatus", trkStatus, &b_trkStatus);
   fChain->SetBranchAddress("nParticle", &nParticle, &b_nParticle);
   fChain->SetBranchAddress("pStatus", pStatus, &b_pStatus);
   fChain->SetBranchAddress("pPId", pPId, &b_pPId);
   fChain->SetBranchAddress("pEta", pEta, &b_pEta);
   fChain->SetBranchAddress("pPhi", pPhi, &b_pPhi);
   fChain->SetBranchAddress("pPt", pPt, &b_pPt);
   fChain->SetBranchAddress("pAcc", pAcc, &b_pAcc);
   fChain->SetBranchAddress("pAccPair", pAccPair, &b_pAccPair);
   fChain->SetBranchAddress("pNRec", pNRec, &b_pNRec);
   fChain->SetBranchAddress("pNHit", pNHit, &b_pNHit);
   fChain->SetBranchAddress("mtrkPt", mtrkPt, &b_mtrkPt);
   fChain->SetBranchAddress("mtrkPtError", mtrkPtError, &b_mtrkPtError);
   fChain->SetBranchAddress("mtrkNHit", mtrkNHit, &b_mtrkNHit);
   fChain->SetBranchAddress("mtrkNlayer", mtrkNlayer, &b_mtrkNlayer);
   fChain->SetBranchAddress("mtrkNlayer3D", mtrkNlayer3D, &b_mtrkNlayer3D);
   fChain->SetBranchAddress("mtrkQual", mtrkQual, &b_mtrkQual);
   fChain->SetBranchAddress("mtrkChi2", mtrkChi2, &b_mtrkChi2);
   fChain->SetBranchAddress("mtrkNdof", mtrkNdof, &b_mtrkNdof);
   fChain->SetBranchAddress("mtrkDz1", mtrkDz1, &b_mtrkDz1);
   fChain->SetBranchAddress("mtrkDzError1", mtrkDzError1, &b_mtrkDzError1);
   fChain->SetBranchAddress("mtrkDxy1", mtrkDxy1, &b_mtrkDxy1);
   fChain->SetBranchAddress("mtrkDxyError1", mtrkDxyError1, &b_mtrkDxyError1);
   fChain->SetBranchAddress("mtrkAlgo", mtrkAlgo, &b_mtrkAlgo);
   fChain->SetBranchAddress("mtrkPfType", mtrkPfType, &b_mtrkPfType);
   fChain->SetBranchAddress("mtrkPfCandPt", mtrkPfCandPt, &b_mtrkPfCandPt);
   fChain->SetBranchAddress("mtrkPfSumEcal", mtrkPfSumEcal, &b_mtrkPfSumEcal);
   fChain->SetBranchAddress("mtrkPfSumHcal", mtrkPfSumHcal, &b_mtrkPfSumHcal);
   Notify();
}

Bool_t ppTrack::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ppTrack::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ppTrack::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
void ppTrack::Close()
{
f->Close();
}
#endif // #ifdef ppTrack_cxx
