//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Mar 15 14:45:26 2014 by ROOT version 5.32/00
// from TTree HltTree/
// found on file: root://eoscms//eos/cms/store/caf/user/yjlee/pp2013/promptReco/PP2013_HiForest_PromptReco_JsonPP_Jet80_PPReco_forestv82.root
//////////////////////////////////////////////////////////

#ifndef HltTree_pp_h
#define HltTree_pp_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class HltTree_pp {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           NL1IsolEm;
   Float_t         L1IsolEmEt[4];   //[NL1IsolEm]
   Float_t         L1IsolEmE[4];   //[NL1IsolEm]
   Float_t         L1IsolEmEta[4];   //[NL1IsolEm]
   Float_t         L1IsolEmPhi[4];   //[NL1IsolEm]
   Int_t           NL1NIsolEm;
   Float_t         L1NIsolEmEt[4];   //[NL1NIsolEm]
   Float_t         L1NIsolEmE[4];   //[NL1NIsolEm]
   Float_t         L1NIsolEmEta[4];   //[NL1NIsolEm]
   Float_t         L1NIsolEmPhi[4];   //[NL1NIsolEm]
   Int_t           NL1Mu;
   Float_t         L1MuPt[4];   //[NL1Mu]
   Float_t         L1MuE[4];   //[NL1Mu]
   Float_t         L1MuEta[4];   //[NL1Mu]
   Float_t         L1MuPhi[4];   //[NL1Mu]
   Int_t           L1MuIsol[4];   //[NL1Mu]
   Int_t           L1MuMip[4];   //[NL1Mu]
   Int_t           L1MuFor[4];   //[NL1Mu]
   Int_t           L1MuRPC[4];   //[NL1Mu]
   Int_t           L1MuQal[4];   //[NL1Mu]
   Int_t           L1MuChg[4];   //[NL1Mu]
   Int_t           NL1CenJet;
   Float_t         L1CenJetEt[4];   //[NL1CenJet]
   Float_t         L1CenJetE[4];   //[NL1CenJet]
   Float_t         L1CenJetEta[4];   //[NL1CenJet]
   Float_t         L1CenJetPhi[4];   //[NL1CenJet]
   Int_t           NL1ForJet;
   Float_t         L1ForJetEt[4];   //[NL1ForJet]
   Float_t         L1ForJetE[4];   //[NL1ForJet]
   Float_t         L1ForJetEta[4];   //[NL1ForJet]
   Float_t         L1ForJetPhi[4];   //[NL1ForJet]
   Int_t           NL1Tau;
   Float_t         L1TauEt[4];   //[NL1Tau]
   Float_t         L1TauE[4];   //[NL1Tau]
   Float_t         L1TauEta[4];   //[NL1Tau]
   Float_t         L1TauPhi[4];   //[NL1Tau]
   Float_t         L1Met;
   Float_t         L1MetPhi;
   Float_t         L1EtTot;
   Float_t         L1Mht;
   Float_t         L1MhtPhi;
   Float_t         L1EtHad;
   Int_t           L1HfRing1EtSumPositiveEta;
   Int_t           L1HfRing2EtSumPositiveEta;
   Int_t           L1HfRing1EtSumNegativeEta;
   Int_t           L1HfRing2EtSumNegativeEta;
   Int_t           L1HfTowerCountPositiveEtaRing1;
   Int_t           L1HfTowerCountNegativeEtaRing1;
   Int_t           L1HfTowerCountPositiveEtaRing2;
   Int_t           L1HfTowerCountNegativeEtaRing2;
   Int_t           Run;
   Int_t           Event;
   Int_t           LumiBlock;
   Int_t           Bx;
   Int_t           Orbit;
   Double_t        AvgInstDelLumi;
   Int_t           HLTriggerFirstPath;
   Int_t           HLTriggerFirstPath_Prescl;
   Int_t           HLT_Activity_Ecal_SC7_v13;
   Int_t           HLT_Activity_Ecal_SC7_v13_Prescl;
   Int_t           HLT_BeamGas_HF_Beam1_v5;
   Int_t           HLT_BeamGas_HF_Beam1_v5_Prescl;
   Int_t           HLT_BeamGas_HF_Beam2_v5;
   Int_t           HLT_BeamGas_HF_Beam2_v5_Prescl;
   Int_t           HLT_BeamHalo_v13;
   Int_t           HLT_BeamHalo_v13_Prescl;
   Int_t           HLT_PAHcalUTCA_v1;
   Int_t           HLT_PAHcalUTCA_v1_Prescl;
   Int_t           HLT_PAHcalPhiSym_v1;
   Int_t           HLT_PAHcalPhiSym_v1_Prescl;
   Int_t           HLT_PAHcalNZS_v1;
   Int_t           HLT_PAHcalNZS_v1_Prescl;
   Int_t           HLT_GlobalRunHPDNoise_v8;
   Int_t           HLT_GlobalRunHPDNoise_v8_Prescl;
   Int_t           HLT_Physics_v5;
   Int_t           HLT_Physics_v5_Prescl;
   Int_t           DST_Physics_v5;
   Int_t           DST_Physics_v5_Prescl;
   Int_t           HLT_DTCalibration_v2;
   Int_t           HLT_DTCalibration_v2_Prescl;
   Int_t           HLT_EcalCalibration_v3;
   Int_t           HLT_EcalCalibration_v3_Prescl;
   Int_t           HLT_HcalCalibration_v3;
   Int_t           HLT_HcalCalibration_v3_Prescl;
   Int_t           HLT_TrackerCalibration_v3;
   Int_t           HLT_TrackerCalibration_v3_Prescl;
   Int_t           HLT_L1SingleMuOpen_AntiBPTX_v7;
   Int_t           HLT_L1SingleMuOpen_AntiBPTX_v7_Prescl;
   Int_t           HLT_L1TrackerCosmics_v7;
   Int_t           HLT_L1TrackerCosmics_v7_Prescl;
   Int_t           AlCa_PAEcalPi0EBonly_v1;
   Int_t           AlCa_PAEcalPi0EBonly_v1_Prescl;
   Int_t           AlCa_PAEcalPi0EEonly_v1;
   Int_t           AlCa_PAEcalPi0EEonly_v1_Prescl;
   Int_t           AlCa_PAEcalEtaEBonly_v1;
   Int_t           AlCa_PAEcalEtaEBonly_v1_Prescl;
   Int_t           AlCa_PAEcalEtaEEonly_v1;
   Int_t           AlCa_PAEcalEtaEEonly_v1_Prescl;
   Int_t           AlCa_EcalPhiSym_v13;
   Int_t           AlCa_EcalPhiSym_v13_Prescl;
   Int_t           AlCa_RPCMuonNoTriggers_v9;
   Int_t           AlCa_RPCMuonNoTriggers_v9_Prescl;
   Int_t           AlCa_RPCMuonNoHits_v9;
   Int_t           AlCa_RPCMuonNoHits_v9_Prescl;
   Int_t           AlCa_RPCMuonNormalisation_v9;
   Int_t           AlCa_RPCMuonNormalisation_v9_Prescl;
   Int_t           AlCa_LumiPixels_v8;
   Int_t           AlCa_LumiPixels_v8_Prescl;
   Int_t           AlCa_LumiPixels_ZeroBias_v4;
   Int_t           AlCa_LumiPixels_ZeroBias_v4_Prescl;
   Int_t           AlCa_LumiPixels_Random_v1;
   Int_t           AlCa_LumiPixels_Random_v1_Prescl;
   Int_t           HLT_PAL1SingleJet16_v1;
   Int_t           HLT_PAL1SingleJet16_v1_Prescl;
   Int_t           HLT_PAL1SingleJet36_v1;
   Int_t           HLT_PAL1SingleJet36_v1_Prescl;
   Int_t           HLT_PASingleForJet15_v1;
   Int_t           HLT_PASingleForJet15_v1_Prescl;
   Int_t           HLT_PASingleForJet25_v1;
   Int_t           HLT_PASingleForJet25_v1_Prescl;
   Int_t           HLT_PAJet20_NoJetID_v1;
   Int_t           HLT_PAJet20_NoJetID_v1_Prescl;
   Int_t           HLT_PAJet40_NoJetID_v1;
   Int_t           HLT_PAJet40_NoJetID_v1_Prescl;
   Int_t           HLT_PAJet60_NoJetID_v1;
   Int_t           HLT_PAJet60_NoJetID_v1_Prescl;
   Int_t           HLT_PAJet80_NoJetID_v1;
   Int_t           HLT_PAJet80_NoJetID_v1_Prescl;
   Int_t           HLT_PAJet100_NoJetID_v1;
   Int_t           HLT_PAJet100_NoJetID_v1_Prescl;
   Int_t           HLT_PAJet120_NoJetID_v1;
   Int_t           HLT_PAJet120_NoJetID_v1_Prescl;
   Int_t           HLT_PAForJet20Eta2_v1;
   Int_t           HLT_PAForJet20Eta2_v1_Prescl;
   Int_t           HLT_PAForJet40Eta2_v1;
   Int_t           HLT_PAForJet40Eta2_v1_Prescl;
   Int_t           HLT_PAForJet60Eta2_v1;
   Int_t           HLT_PAForJet60Eta2_v1_Prescl;
   Int_t           HLT_PAForJet80Eta2_v1;
   Int_t           HLT_PAForJet80Eta2_v1_Prescl;
   Int_t           HLT_PAForJet100Eta2_v1;
   Int_t           HLT_PAForJet100Eta2_v1_Prescl;
   Int_t           HLT_PAForJet20Eta3_v1;
   Int_t           HLT_PAForJet20Eta3_v1_Prescl;
   Int_t           HLT_PAForJet40Eta3_v1;
   Int_t           HLT_PAForJet40Eta3_v1_Prescl;
   Int_t           HLT_PAForJet60Eta3_v1;
   Int_t           HLT_PAForJet60Eta3_v1_Prescl;
   Int_t           HLT_PAForJet80Eta3_v1;
   Int_t           HLT_PAForJet80Eta3_v1_Prescl;
   Int_t           HLT_PAForJet100Eta3_v1;
   Int_t           HLT_PAForJet100Eta3_v1_Prescl;
   Int_t           HLT_PATripleJet20_20_20_v1;
   Int_t           HLT_PATripleJet20_20_20_v1_Prescl;
   Int_t           HLT_PATripleJet40_20_20_v1;
   Int_t           HLT_PATripleJet40_20_20_v1_Prescl;
   Int_t           HLT_PATripleJet60_20_20_v1;
   Int_t           HLT_PATripleJet60_20_20_v1_Prescl;
   Int_t           HLT_PATripleJet80_20_20_v1;
   Int_t           HLT_PATripleJet80_20_20_v1_Prescl;
   Int_t           HLT_PATripleJet100_20_20_v1;
   Int_t           HLT_PATripleJet100_20_20_v1_Prescl;
   Int_t           HLT_PAJet40ETM30_v1;
   Int_t           HLT_PAJet40ETM30_v1_Prescl;
   Int_t           HLT_PAJet60ETM30_v1;
   Int_t           HLT_PAJet60ETM30_v1_Prescl;
   Int_t           HLT_PAL1DoubleMu0_v1;
   Int_t           HLT_PAL1DoubleMu0_v1_Prescl;
   Int_t           HLT_PADimuon0_NoVertexing_v1;
   Int_t           HLT_PADimuon0_NoVertexing_v1_Prescl;
   Int_t           HLT_PAL1DoubleMu0_HighQ_v1;
   Int_t           HLT_PAL1DoubleMu0_HighQ_v1_Prescl;
   Int_t           HLT_PAL1DoubleMuOpen_v1;
   Int_t           HLT_PAL1DoubleMuOpen_v1_Prescl;
   Int_t           HLT_PAL2DoubleMu3_v1;
   Int_t           HLT_PAL2DoubleMu3_v1_Prescl;
   Int_t           HLT_PAMu3_v1;
   Int_t           HLT_PAMu3_v1_Prescl;
   Int_t           HLT_PAMu7_v1;
   Int_t           HLT_PAMu7_v1_Prescl;
   Int_t           HLT_PAMu12_v1;
   Int_t           HLT_PAMu12_v1_Prescl;
   Int_t           HLT_Mu15_eta2p1_v5;
   Int_t           HLT_Mu15_eta2p1_v5_Prescl;
   Int_t           HLT_PABTagMu_Jet20_Mu4_v1;
   Int_t           HLT_PABTagMu_Jet20_Mu4_v1_Prescl;
   Int_t           HLT_PAMu3PFJet20_v1;
   Int_t           HLT_PAMu3PFJet20_v1_Prescl;
   Int_t           HLT_PAMu3PFJet40_v1;
   Int_t           HLT_PAMu3PFJet40_v1_Prescl;
   Int_t           HLT_PAMu7PFJet20_v1;
   Int_t           HLT_PAMu7PFJet20_v1_Prescl;
   Int_t           HLT_PAPhoton10_NoCaloIdVL_v1;
   Int_t           HLT_PAPhoton10_NoCaloIdVL_v1_Prescl;
   Int_t           HLT_PAPhoton15_NoCaloIdVL_v1;
   Int_t           HLT_PAPhoton15_NoCaloIdVL_v1_Prescl;
   Int_t           HLT_PAPhoton20_NoCaloIdVL_v1;
   Int_t           HLT_PAPhoton20_NoCaloIdVL_v1_Prescl;
   Int_t           HLT_PAPhoton30_NoCaloIdVL_v1;
   Int_t           HLT_PAPhoton30_NoCaloIdVL_v1_Prescl;
   Int_t           HLT_PAPhoton40_NoCaloIdVL_v1;
   Int_t           HLT_PAPhoton40_NoCaloIdVL_v1_Prescl;
   Int_t           HLT_PAPhoton60_NoCaloIdVL_v1;
   Int_t           HLT_PAPhoton60_NoCaloIdVL_v1_Prescl;
   Int_t           HLT_Ele22_CaloIdL_CaloIsoVL_v6;
   Int_t           HLT_Ele22_CaloIdL_CaloIsoVL_v6_Prescl;
   Int_t           HLT_PAPhoton10_TightCaloIdVL_v1;
   Int_t           HLT_PAPhoton10_TightCaloIdVL_v1_Prescl;
   Int_t           HLT_PAPhoton15_TightCaloIdVL_v1;
   Int_t           HLT_PAPhoton15_TightCaloIdVL_v1_Prescl;
   Int_t           HLT_PAPhoton20_TightCaloIdVL_v1;
   Int_t           HLT_PAPhoton20_TightCaloIdVL_v1_Prescl;
   Int_t           HLT_PAPhoton30_TightCaloIdVL_v1;
   Int_t           HLT_PAPhoton30_TightCaloIdVL_v1_Prescl;
   Int_t           HLT_PAPhoton40_TightCaloIdVL_v1;
   Int_t           HLT_PAPhoton40_TightCaloIdVL_v1_Prescl;
   Int_t           HLT_PAPhoton10_Photon10_NoCaloIdVL_v1;
   Int_t           HLT_PAPhoton10_Photon10_NoCaloIdVL_v1_Prescl;
   Int_t           HLT_PAPhoton15_Photon10_NoCaloIdVL_v1;
   Int_t           HLT_PAPhoton15_Photon10_NoCaloIdVL_v1_Prescl;
   Int_t           HLT_PAPhoton20_Photon15_NoCaloIdVL_v1;
   Int_t           HLT_PAPhoton20_Photon15_NoCaloIdVL_v1_Prescl;
   Int_t           HLT_PAPhoton20_Photon20_NoCaloIdVL_v1;
   Int_t           HLT_PAPhoton20_Photon20_NoCaloIdVL_v1_Prescl;
   Int_t           HLT_PAPhoton30_Photon30_NoCaloIdVL_v1;
   Int_t           HLT_PAPhoton30_Photon30_NoCaloIdVL_v1_Prescl;
   Int_t           HLT_PAPhoton10_Photon10_TightCaloIdVL_v1;
   Int_t           HLT_PAPhoton10_Photon10_TightCaloIdVL_v1_Prescl;
   Int_t           HLT_PAPhoton15_Photon10_TightCaloIdVL_v1;
   Int_t           HLT_PAPhoton15_Photon10_TightCaloIdVL_v1_Prescl;
   Int_t           HLT_PAPhoton20_Photon15_TightCaloIdVL_v1;
   Int_t           HLT_PAPhoton20_Photon15_TightCaloIdVL_v1_Prescl;
   Int_t           HLT_PAPixelTracks_Multiplicity100_v2;
   Int_t           HLT_PAPixelTracks_Multiplicity100_v2_Prescl;
   Int_t           HLT_PAPixelTracks_Multiplicity130_v2;
   Int_t           HLT_PAPixelTracks_Multiplicity130_v2_Prescl;
   Int_t           HLT_PAPixelTracks_Multiplicity160_v2;
   Int_t           HLT_PAPixelTracks_Multiplicity160_v2_Prescl;
   Int_t           HLT_PAPixelTracks_Multiplicity190_v2;
   Int_t           HLT_PAPixelTracks_Multiplicity190_v2_Prescl;
   Int_t           HLT_PAPixelTracks_Multiplicity220_v2;
   Int_t           HLT_PAPixelTracks_Multiplicity220_v2_Prescl;
   Int_t           HLT_PAPixelTrackMultiplicity100_FullTrack12_v2;
   Int_t           HLT_PAPixelTrackMultiplicity100_FullTrack12_v2_Prescl;
   Int_t           HLT_PAPixelTrackMultiplicity130_FullTrack12_v2;
   Int_t           HLT_PAPixelTrackMultiplicity130_FullTrack12_v2_Prescl;
   Int_t           HLT_PAPixelTrackMultiplicity160_FullTrack12_v2;
   Int_t           HLT_PAPixelTrackMultiplicity160_FullTrack12_v2_Prescl;
   Int_t           HLT_PAFullTrack12_v2;
   Int_t           HLT_PAFullTrack12_v2_Prescl;
   Int_t           HLT_PAFullTrack20_v2;
   Int_t           HLT_PAFullTrack20_v2_Prescl;
   Int_t           HLT_PAFullTrack30_v2;
   Int_t           HLT_PAFullTrack30_v2_Prescl;
   Int_t           HLT_PAFullTrack50_v2;
   Int_t           HLT_PAFullTrack50_v2_Prescl;
   Int_t           HLT_PAPixelTrackMultiplicity140_Jet80_NoJetID_v2;
   Int_t           HLT_PAPixelTrackMultiplicity140_Jet80_NoJetID_v2_Prescl;
   Int_t           HLT_PAPixelTrackMultiplicity100_L2DoubleMu3_v1;
   Int_t           HLT_PAPixelTrackMultiplicity100_L2DoubleMu3_v1_Prescl;
   Int_t           HLT_PPPixelTracks_Multiplicity55_v1;
   Int_t           HLT_PPPixelTracks_Multiplicity55_v1_Prescl;
   Int_t           HLT_PPPixelTracks_Multiplicity70_v1;
   Int_t           HLT_PPPixelTracks_Multiplicity70_v1_Prescl;
   Int_t           HLT_PPPixelTracks_Multiplicity85_v1;
   Int_t           HLT_PPPixelTracks_Multiplicity85_v1_Prescl;
   Int_t           HLT_PPPixelTrackMultiplicity55_FullTrack12_v1;
   Int_t           HLT_PPPixelTrackMultiplicity55_FullTrack12_v1_Prescl;
   Int_t           HLT_PPPixelTrackMultiplicity70_FullTrack12_v1;
   Int_t           HLT_PPPixelTrackMultiplicity70_FullTrack12_v1_Prescl;
   Int_t           HLT_PPL1DoubleJetC36_v1;
   Int_t           HLT_PPL1DoubleJetC36_v1_Prescl;
   Int_t           HLT_PATech35_v1;
   Int_t           HLT_PATech35_v1_Prescl;
   Int_t           HLT_PATech35_HFSumET100_v3;
   Int_t           HLT_PATech35_HFSumET100_v3_Prescl;
   Int_t           HLT_PAHFSumET100_v3;
   Int_t           HLT_PAHFSumET100_v3_Prescl;
   Int_t           HLT_PAHFSumET140_v3;
   Int_t           HLT_PAHFSumET140_v3_Prescl;
   Int_t           HLT_PAHFSumET170_v3;
   Int_t           HLT_PAHFSumET170_v3_Prescl;
   Int_t           HLT_PAHFSumET210_v3;
   Int_t           HLT_PAHFSumET210_v3_Prescl;
   Int_t           HLT_PARomanPots_Tech52_v1;
   Int_t           HLT_PARomanPots_Tech52_v1_Prescl;
   Int_t           HLT_PAL1Tech53_MB_v1;
   Int_t           HLT_PAL1Tech53_MB_v1_Prescl;
   Int_t           HLT_PAL1Tech53_MB_SingleTrack_v1;
   Int_t           HLT_PAL1Tech53_MB_SingleTrack_v1_Prescl;
   Int_t           HLT_PAL1Tech54_ZeroBias_v1;
   Int_t           HLT_PAL1Tech54_ZeroBias_v1_Prescl;
   Int_t           HLT_PAT1minbias_Tech55_v1;
   Int_t           HLT_PAT1minbias_Tech55_v1_Prescl;
   Int_t           HLT_PAL1Tech_HBHEHO_totalOR_v1;
   Int_t           HLT_PAL1Tech_HBHEHO_totalOR_v1_Prescl;
   Int_t           HLT_PAL1Tech63_CASTORHaloMuon_v1;
   Int_t           HLT_PAL1Tech63_CASTORHaloMuon_v1_Prescl;
   Int_t           HLT_PACastorEmTotemLowMultiplicity_v1;
   Int_t           HLT_PACastorEmTotemLowMultiplicity_v1_Prescl;
   Int_t           HLT_PACastorEmNotHfCoincidencePm_v1;
   Int_t           HLT_PACastorEmNotHfCoincidencePm_v1_Prescl;
   Int_t           HLT_PACastorEmNotHfSingleChannel_v1;
   Int_t           HLT_PACastorEmNotHfSingleChannel_v1_Prescl;
   Int_t           HLT_PAL1CastorTotalTotemLowMultiplicity_v1;
   Int_t           HLT_PAL1CastorTotalTotemLowMultiplicity_v1_Prescl;
   Int_t           HLT_PAMinBiasHF_v1;
   Int_t           HLT_PAMinBiasHF_v1_Prescl;
   Int_t           HLT_PAMinBiasHF_OR_v1;
   Int_t           HLT_PAMinBiasHF_OR_v1_Prescl;
   Int_t           HLT_PAMinBiasBHC_v1;
   Int_t           HLT_PAMinBiasBHC_v1_Prescl;
   Int_t           HLT_PAMinBiasBHC_OR_v1;
   Int_t           HLT_PAMinBiasBHC_OR_v1_Prescl;
   Int_t           HLT_PAMinBiasHfOrBHC_v1;
   Int_t           HLT_PAMinBiasHfOrBHC_v1_Prescl;
   Int_t           HLT_PABptxPlusNotBptxMinus_v1;
   Int_t           HLT_PABptxPlusNotBptxMinus_v1_Prescl;
   Int_t           HLT_PABptxMinusNotBptxPlus_v1;
   Int_t           HLT_PABptxMinusNotBptxPlus_v1_Prescl;
   Int_t           HLT_PAZeroBias_v1;
   Int_t           HLT_PAZeroBias_v1_Prescl;
   Int_t           HLT_PAZeroBiasPixel_SingleTrack_v1;
   Int_t           HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl;
   Int_t           HLT_PAHFOR_SingleTrack_v1;
   Int_t           HLT_PAHFOR_SingleTrack_v1_Prescl;
   Int_t           HLT_PAZeroBiasPixel_DoubleTrack_v1;
   Int_t           HLT_PAZeroBiasPixel_DoubleTrack_v1_Prescl;
   Int_t           HLT_PADoubleMu4_Acoplanarity03_v1;
   Int_t           HLT_PADoubleMu4_Acoplanarity03_v1_Prescl;
   Int_t           HLT_PAExclDijet35_HFOR_v1;
   Int_t           HLT_PAExclDijet35_HFOR_v1_Prescl;
   Int_t           HLT_PAExclDijet35_HFAND_v1;
   Int_t           HLT_PAExclDijet35_HFAND_v1_Prescl;
   Int_t           HLT_PAL1DoubleEG3_FwdVeto_v1;
   Int_t           HLT_PAL1DoubleEG3_FwdVeto_v1_Prescl;
   Int_t           HLT_PAL1SingleJet52_TotemDiffractive_v1;
   Int_t           HLT_PAL1SingleJet52_TotemDiffractive_v1_Prescl;
   Int_t           HLT_PAL1SingleMu20_TotemDiffractive_v1;
   Int_t           HLT_PAL1SingleMu20_TotemDiffractive_v1_Prescl;
   Int_t           HLT_PAL1SingleEG20_TotemDiffractive_v1;
   Int_t           HLT_PAL1SingleEG20_TotemDiffractive_v1_Prescl;
   Int_t           HLT_PAL1DoubleJet20_TotemDiffractive_v1;
   Int_t           HLT_PAL1DoubleJet20_TotemDiffractive_v1_Prescl;
   Int_t           HLT_PAL1DoubleJetC36_TotemDiffractive_v1;
   Int_t           HLT_PAL1DoubleJetC36_TotemDiffractive_v1_Prescl;
   Int_t           HLT_PAL1DoubleMu5_TotemDiffractive_v1;
   Int_t           HLT_PAL1DoubleMu5_TotemDiffractive_v1_Prescl;
   Int_t           HLT_PAL1DoubleEG5_TotemDiffractive_v1;
   Int_t           HLT_PAL1DoubleEG5_TotemDiffractive_v1_Prescl;
   Int_t           HLT_PADoubleJet20_ForwardBackward_v1;
   Int_t           HLT_PADoubleJet20_ForwardBackward_v1_Prescl;
   Int_t           HLT_PAMu7_Ele7_CaloIdT_CaloIsoVL_v1;
   Int_t           HLT_PAMu7_Ele7_CaloIdT_CaloIsoVL_v1_Prescl;
   Int_t           HLT_PARandom_v1;
   Int_t           HLT_PARandom_v1_Prescl;
   Int_t           DQM_FEDIntegrity_v11;
   Int_t           DQM_FEDIntegrity_v11_Prescl;
   Int_t           HLT_LogMonitor_v4;
   Int_t           HLT_LogMonitor_v4_Prescl;
   Int_t           HLTriggerFinalPath;
   Int_t           HLTriggerFinalPath_Prescl;
   Int_t           L1_AlwaysTrue;
   Int_t           L1_AlwaysTrue_Prescl;
   Int_t           L1_BeamGas_Hf_BptxMinusPostQuiet;
   Int_t           L1_BeamGas_Hf_BptxMinusPostQuiet_Prescl;
   Int_t           L1_BeamGas_Hf_BptxPlusPostQuiet;
   Int_t           L1_BeamGas_Hf_BptxPlusPostQuiet_Prescl;
   Int_t           L1_BeamHalo;
   Int_t           L1_BeamHalo_Prescl;
   Int_t           L1_BptxMinus_NotBptxPlus;
   Int_t           L1_BptxMinus_NotBptxPlus_Prescl;
   Int_t           L1_BptxPlus_NotBptxMinus;
   Int_t           L1_BptxPlus_NotBptxMinus_Prescl;
   Int_t           L1_BscHighMultiplicity_BptxAND;
   Int_t           L1_BscHighMultiplicity_BptxAND_Prescl;
   Int_t           L1_BscMinBiasOR_BptxAND;
   Int_t           L1_BscMinBiasOR_BptxAND_Prescl;
   Int_t           L1_BscMinBiasThreshold1_BptxAND;
   Int_t           L1_BscMinBiasThreshold1_BptxAND_Prescl;
   Int_t           L1_BscMinBiasThreshold2_BptxAND;
   Int_t           L1_BscMinBiasThreshold2_BptxAND_Prescl;
   Int_t           L1_CastorEm_NotHcalHfCoincidencePm;
   Int_t           L1_CastorEm_NotHcalHfCoincidencePm_Prescl;
   Int_t           L1_CastorEm_NotHcalHfSingleChannel;
   Int_t           L1_CastorEm_NotHcalHfSingleChannel_Prescl;
   Int_t           L1_CastorEm_TotemLowMultiplicity;
   Int_t           L1_CastorEm_TotemLowMultiplicity_Prescl;
   Int_t           L1_CastorTotalEnergy_TotemLowMultiplicity;
   Int_t           L1_CastorTotalEnergy_TotemLowMultiplicity_Prescl;
   Int_t           L1_DoubleEG3_FwdVeto;
   Int_t           L1_DoubleEG3_FwdVeto_Prescl;
   Int_t           L1_DoubleEG5;
   Int_t           L1_DoubleEG5_Prescl;
   Int_t           L1_DoubleEG5_TotemDiffractive;
   Int_t           L1_DoubleEG5_TotemDiffractive_Prescl;
   Int_t           L1_DoubleEG6_HTT100;
   Int_t           L1_DoubleEG6_HTT100_Prescl;
   Int_t           L1_DoubleEG6_HTT125;
   Int_t           L1_DoubleEG6_HTT125_Prescl;
   Int_t           L1_DoubleEG_13_7;
   Int_t           L1_DoubleEG_13_7_Prescl;
   Int_t           L1_DoubleForJet16_EtaOpp;
   Int_t           L1_DoubleForJet16_EtaOpp_Prescl;
   Int_t           L1_DoubleJet20;
   Int_t           L1_DoubleJet20_Prescl;
   Int_t           L1_DoubleJet20_TotemDiffractive;
   Int_t           L1_DoubleJet20_TotemDiffractive_Prescl;
   Int_t           L1_DoubleJet24_v1;
   Int_t           L1_DoubleJet24_v1_Prescl;
   Int_t           L1_DoubleJet36_Central;
   Int_t           L1_DoubleJet36_Central_Prescl;
   Int_t           L1_DoubleJet52_Central;
   Int_t           L1_DoubleJet52_Central_Prescl;
   Int_t           L1_DoubleJetC36_TotemDiffractive;
   Int_t           L1_DoubleJetC36_TotemDiffractive_Prescl;
   Int_t           L1_DoubleJetC44_ETM30;
   Int_t           L1_DoubleJetC44_ETM30_Prescl;
   Int_t           L1_DoubleJetC56;
   Int_t           L1_DoubleJetC56_Prescl;
   Int_t           L1_DoubleJetC56_Eta1p74_WdEta4;
   Int_t           L1_DoubleJetC56_Eta1p74_WdEta4_Prescl;
   Int_t           L1_DoubleMu0;
   Int_t           L1_DoubleMu0_Prescl;
   Int_t           L1_DoubleMu0_HighQ_EtaCuts;
   Int_t           L1_DoubleMu0_HighQ_EtaCuts_Prescl;
   Int_t           L1_DoubleMu3p5_EG5;
   Int_t           L1_DoubleMu3p5_EG5_Prescl;
   Int_t           L1_DoubleMu5_EG5;
   Int_t           L1_DoubleMu5_EG5_Prescl;
   Int_t           L1_DoubleMu5_TotemDiffractive;
   Int_t           L1_DoubleMu5_TotemDiffractive_Prescl;
   Int_t           L1_DoubleMu5_v1;
   Int_t           L1_DoubleMu5_v1_Prescl;
   Int_t           L1_DoubleMuOpen_BptxAND;
   Int_t           L1_DoubleMuOpen_BptxAND_Prescl;
   Int_t           L1_DoubleMu_10_3p5;
   Int_t           L1_DoubleMu_10_3p5_Prescl;
   Int_t           L1_DoubleMu_10_Open;
   Int_t           L1_DoubleMu_10_Open_Prescl;
   Int_t           L1_DoubleMu_12_5;
   Int_t           L1_DoubleMu_12_5_Prescl;
   Int_t           L1_DoubleMu_3er_0er_HighQ_WdEta22;
   Int_t           L1_DoubleMu_3er_0er_HighQ_WdEta22_Prescl;
   Int_t           L1_DoubleMu_5er_0er_HighQ_WdEta22;
   Int_t           L1_DoubleMu_5er_0er_HighQ_WdEta22_Prescl;
   Int_t           L1_EG8_DoubleJetC20;
   Int_t           L1_EG8_DoubleJetC20_Prescl;
   Int_t           L1_ETM100;
   Int_t           L1_ETM100_Prescl;
   Int_t           L1_ETM30;
   Int_t           L1_ETM30_Prescl;
   Int_t           L1_ETM36;
   Int_t           L1_ETM36_Prescl;
   Int_t           L1_ETM40;
   Int_t           L1_ETM40_Prescl;
   Int_t           L1_ETM50;
   Int_t           L1_ETM50_Prescl;
   Int_t           L1_ETM70;
   Int_t           L1_ETM70_Prescl;
   Int_t           L1_ETT140;
   Int_t           L1_ETT140_Prescl;
   Int_t           L1_ETT20_BptxAND;
   Int_t           L1_ETT20_BptxAND_Prescl;
   Int_t           L1_ETT300;
   Int_t           L1_ETT300_Prescl;
   Int_t           L1_ETT40;
   Int_t           L1_ETT40_Prescl;
   Int_t           L1_ETT60;
   Int_t           L1_ETT60_Prescl;
   Int_t           L1_ETT80;
   Int_t           L1_ETT80_Prescl;
   Int_t           L1_HTT100;
   Int_t           L1_HTT100_Prescl;
   Int_t           L1_HTT125;
   Int_t           L1_HTT125_Prescl;
   Int_t           L1_HTT150;
   Int_t           L1_HTT150_Prescl;
   Int_t           L1_HTT175;
   Int_t           L1_HTT175_Prescl;
   Int_t           L1_HTT200;
   Int_t           L1_HTT200_Prescl;
   Int_t           L1_HTT50;
   Int_t           L1_HTT50_Prescl;
   Int_t           L1_HTT75;
   Int_t           L1_HTT75_Prescl;
   Int_t           L1_HcalHfCoincidencePm_BptxAND_v1;
   Int_t           L1_HcalHfCoincidencePm_BptxAND_v1_Prescl;
   Int_t           L1_HcalHfSingleChannel_BptxAND;
   Int_t           L1_HcalHfSingleChannel_BptxAND_Prescl;
   Int_t           L1_HcalHfSingleChannel_BptxAND_Instance1;
   Int_t           L1_HcalHfSingleChannel_BptxAND_Instance1_Prescl;
   Int_t           L1_Mu10er_JetC32;
   Int_t           L1_Mu10er_JetC32_Prescl;
   Int_t           L1_Mu12_EG7;
   Int_t           L1_Mu12_EG7_Prescl;
   Int_t           L1_Mu3_Jet16;
   Int_t           L1_Mu3_Jet16_Prescl;
   Int_t           L1_Mu3_Jet36;
   Int_t           L1_Mu3_Jet36_Prescl;
   Int_t           L1_Mu3_JetC16_WdEtaPhi2;
   Int_t           L1_Mu3_JetC16_WdEtaPhi2_Prescl;
   Int_t           L1_Mu3_JetC52_WdEtaPhi2;
   Int_t           L1_Mu3_JetC52_WdEtaPhi2_Prescl;
   Int_t           L1_Mu5_DoubleEG5;
   Int_t           L1_Mu5_DoubleEG5_Prescl;
   Int_t           L1_Mu5_DoubleEG6;
   Int_t           L1_Mu5_DoubleEG6_Prescl;
   Int_t           L1_Mu7_Jet16;
   Int_t           L1_Mu7_Jet16_Prescl;
   Int_t           L1_Mu8_DoubleJetC20;
   Int_t           L1_Mu8_DoubleJetC20_Prescl;
   Int_t           L1_MuOpen_EG12;
   Int_t           L1_MuOpen_EG12_Prescl;
   Int_t           L1_MuOpen_EG5;
   Int_t           L1_MuOpen_EG5_Prescl;
   Int_t           L1_QuadJetC32;
   Int_t           L1_QuadJetC32_Prescl;
   Int_t           L1_QuadJetC36;
   Int_t           L1_QuadJetC36_Prescl;
   Int_t           L1_QuadJetC40;
   Int_t           L1_QuadJetC40_Prescl;
   Int_t           L1_SingleEG12;
   Int_t           L1_SingleEG12_Prescl;
   Int_t           L1_SingleEG18er;
   Int_t           L1_SingleEG18er_Prescl;
   Int_t           L1_SingleEG20;
   Int_t           L1_SingleEG20_Prescl;
   Int_t           L1_SingleEG20_TotemDiffractive;
   Int_t           L1_SingleEG20_TotemDiffractive_Prescl;
   Int_t           L1_SingleEG22;
   Int_t           L1_SingleEG22_Prescl;
   Int_t           L1_SingleEG24;
   Int_t           L1_SingleEG24_Prescl;
   Int_t           L1_SingleEG30;
   Int_t           L1_SingleEG30_Prescl;
   Int_t           L1_SingleEG5_BptxAND;
   Int_t           L1_SingleEG5_BptxAND_Prescl;
   Int_t           L1_SingleEG7;
   Int_t           L1_SingleEG7_Prescl;
   Int_t           L1_SingleForJet16;
   Int_t           L1_SingleForJet16_Prescl;
   Int_t           L1_SingleIsoEG18er;
   Int_t           L1_SingleIsoEG18er_Prescl;
   Int_t           L1_SingleIsoEG20er;
   Int_t           L1_SingleIsoEG20er_Prescl;
   Int_t           L1_SingleJet128;
   Int_t           L1_SingleJet128_Prescl;
   Int_t           L1_SingleJet12_BptxAND;
   Int_t           L1_SingleJet12_BptxAND_Prescl;
   Int_t           L1_SingleJet16_BptxAND;
   Int_t           L1_SingleJet16_BptxAND_Prescl;
   Int_t           L1_SingleJet16_FwdVeto5;
   Int_t           L1_SingleJet16_FwdVeto5_Prescl;
   Int_t           L1_SingleJet20_Central_NotBptxOR;
   Int_t           L1_SingleJet20_Central_NotBptxOR_Prescl;
   Int_t           L1_SingleJet36;
   Int_t           L1_SingleJet36_Prescl;
   Int_t           L1_SingleJet36_FwdVeto5;
   Int_t           L1_SingleJet36_FwdVeto5_Prescl;
   Int_t           L1_SingleJet52;
   Int_t           L1_SingleJet52_Prescl;
   Int_t           L1_SingleJet52_TotemDiffractive;
   Int_t           L1_SingleJet52_TotemDiffractive_Prescl;
   Int_t           L1_SingleJet68;
   Int_t           L1_SingleJet68_Prescl;
   Int_t           L1_SingleJet92;
   Int_t           L1_SingleJet92_Prescl;
   Int_t           L1_SingleJetC32_NotBptxOR;
   Int_t           L1_SingleJetC32_NotBptxOR_Prescl;
   Int_t           L1_SingleMu12;
   Int_t           L1_SingleMu12_Prescl;
   Int_t           L1_SingleMu12er;
   Int_t           L1_SingleMu12er_Prescl;
   Int_t           L1_SingleMu14_Eta2p1;
   Int_t           L1_SingleMu14_Eta2p1_Prescl;
   Int_t           L1_SingleMu16;
   Int_t           L1_SingleMu16_Prescl;
   Int_t           L1_SingleMu16_Eta2p1;
   Int_t           L1_SingleMu16_Eta2p1_Prescl;
   Int_t           L1_SingleMu18er;
   Int_t           L1_SingleMu18er_Prescl;
   Int_t           L1_SingleMu20;
   Int_t           L1_SingleMu20_Prescl;
   Int_t           L1_SingleMu20_TotemDiffractive;
   Int_t           L1_SingleMu20_TotemDiffractive_Prescl;
   Int_t           L1_SingleMu20er;
   Int_t           L1_SingleMu20er_Prescl;
   Int_t           L1_SingleMu25er;
   Int_t           L1_SingleMu25er_Prescl;
   Int_t           L1_SingleMu3;
   Int_t           L1_SingleMu3_Prescl;
   Int_t           L1_SingleMu6_NotBptxOR;
   Int_t           L1_SingleMu6_NotBptxOR_Prescl;
   Int_t           L1_SingleMu7;
   Int_t           L1_SingleMu7_Prescl;
   Int_t           L1_SingleMuBeamHalo;
   Int_t           L1_SingleMuBeamHalo_Prescl;
   Int_t           L1_SingleMuOpen;
   Int_t           L1_SingleMuOpen_Prescl;
   Int_t           L1_TotemMinBias;
   Int_t           L1_TotemMinBias_Prescl;
   Int_t           L1_TripleEG7;
   Int_t           L1_TripleEG7_Prescl;
   Int_t           L1_TripleEG_12_7_5;
   Int_t           L1_TripleEG_12_7_5_Prescl;
   Int_t           L1_TripleJetC_52_28_28;
   Int_t           L1_TripleJetC_52_28_28_Prescl;
   Int_t           L1_TripleJet_64_44_24_VBF;
   Int_t           L1_TripleJet_64_44_24_VBF_Prescl;
   Int_t           L1_TripleJet_64_48_28_VBF;
   Int_t           L1_TripleJet_64_48_28_VBF_Prescl;
   Int_t           L1_ZdcCaloPlus_TotemDiffractive_QElastic;
   Int_t           L1_ZdcCaloPlus_TotemDiffractive_QElastic_Prescl;
   Int_t           L1_ZeroBias;
   Int_t           L1_ZeroBias_Prescl;
   Int_t           L1Tech_BPTX_PreBPTX_v0;
   Int_t           L1Tech_BPTX_PreBPTX_v0_Prescl;
   Int_t           L1Tech_BPTX_minus_v0;
   Int_t           L1Tech_BPTX_minus_v0_Prescl;
   Int_t           L1Tech_BPTX_minus_AND_not_plus_v0;
   Int_t           L1Tech_BPTX_minus_AND_not_plus_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_v0;
   Int_t           L1Tech_BPTX_plus_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_AND_NOT_minus_v0;
   Int_t           L1Tech_BPTX_plus_AND_NOT_minus_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_AND_minus_v0;
   Int_t           L1Tech_BPTX_plus_AND_minus_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_AND_minus_instance1_v0;
   Int_t           L1Tech_BPTX_plus_AND_minus_instance1_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_OR_minus_v0;
   Int_t           L1Tech_BPTX_plus_OR_minus_v0_Prescl;
   Int_t           L1Tech_BPTX_quiet_v0;
   Int_t           L1Tech_BPTX_quiet_v0_Prescl;
   Int_t           L1Tech_BSC_HighMultiplicity_v0;
   Int_t           L1Tech_BSC_HighMultiplicity_v0_Prescl;
   Int_t           L1Tech_BSC_halo_beam1_inner_v0;
   Int_t           L1Tech_BSC_halo_beam1_inner_v0_Prescl;
   Int_t           L1Tech_BSC_halo_beam1_outer_v0;
   Int_t           L1Tech_BSC_halo_beam1_outer_v0_Prescl;
   Int_t           L1Tech_BSC_halo_beam2_inner_v0;
   Int_t           L1Tech_BSC_halo_beam2_inner_v0_Prescl;
   Int_t           L1Tech_BSC_halo_beam2_outer_v0;
   Int_t           L1Tech_BSC_halo_beam2_outer_v0_Prescl;
   Int_t           L1Tech_BSC_minBias_OR_v0;
   Int_t           L1Tech_BSC_minBias_OR_v0_Prescl;
   Int_t           L1Tech_BSC_minBias_inner_threshold1_v0;
   Int_t           L1Tech_BSC_minBias_inner_threshold1_v0_Prescl;
   Int_t           L1Tech_BSC_minBias_inner_threshold2_v0;
   Int_t           L1Tech_BSC_minBias_inner_threshold2_v0_Prescl;
   Int_t           L1Tech_BSC_minBias_threshold1_v0;
   Int_t           L1Tech_BSC_minBias_threshold1_v0_Prescl;
   Int_t           L1Tech_BSC_minBias_threshold2_v0;
   Int_t           L1Tech_BSC_minBias_threshold2_v0_Prescl;
   Int_t           L1Tech_BSC_splash_beam1_v0;
   Int_t           L1Tech_BSC_splash_beam1_v0_Prescl;
   Int_t           L1Tech_BSC_splash_beam2_v0;
   Int_t           L1Tech_BSC_splash_beam2_v0_Prescl;
   Int_t           L1Tech_CASTOR_0_v0;
   Int_t           L1Tech_CASTOR_0_v0_Prescl;
   Int_t           L1Tech_CASTOR_EM_v0;
   Int_t           L1Tech_CASTOR_EM_v0_Prescl;
   Int_t           L1Tech_CASTOR_HaloMuon_v0;
   Int_t           L1Tech_CASTOR_HaloMuon_v0_Prescl;
   Int_t           L1Tech_CASTOR_TotalEnergy_v0;
   Int_t           L1Tech_CASTOR_TotalEnergy_v0_Prescl;
   Int_t           L1Tech_DT_GlobalOR_v0;
   Int_t           L1Tech_DT_GlobalOR_v0_Prescl;
   Int_t           L1Tech_FSC_St3Sect45_downLeft_v0;
   Int_t           L1Tech_FSC_St3Sect45_downLeft_v0_Prescl;
   Int_t           L1Tech_FSC_St3Sect45_downRight_v0;
   Int_t           L1Tech_FSC_St3Sect45_downRight_v0_Prescl;
   Int_t           L1Tech_FSC_St3Sect45_uppLeft_v0;
   Int_t           L1Tech_FSC_St3Sect45_uppLeft_v0_Prescl;
   Int_t           L1Tech_FSC_St3Sect45_uppRight_v0;
   Int_t           L1Tech_FSC_St3Sect45_uppRight_v0_Prescl;
   Int_t           L1Tech_HCAL_HBHE_totalOR_v0;
   Int_t           L1Tech_HCAL_HBHE_totalOR_v0_Prescl;
   Int_t           L1Tech_HCAL_HF_MMP_or_MPP_v1;
   Int_t           L1Tech_HCAL_HF_MMP_or_MPP_v1_Prescl;
   Int_t           L1Tech_HCAL_HF_coincidence_PM_v2;
   Int_t           L1Tech_HCAL_HF_coincidence_PM_v2_Prescl;
   Int_t           L1Tech_HCAL_HF_single_channel_v0;
   Int_t           L1Tech_HCAL_HF_single_channel_v0_Prescl;
   Int_t           L1Tech_HCAL_HO_totalOR_v0;
   Int_t           L1Tech_HCAL_HO_totalOR_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_RB0_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_RB0_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_RBminus1_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_RBminus1_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_RBminus2_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_RBminus2_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_RBplus1_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_RBplus1_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_RBplus2_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_RBplus2_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_RBst1_collisions_v0;
   Int_t           L1Tech_RPC_TTU_RBst1_collisions_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_barrel_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_barrel_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_pointing_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_pointing_Cosmics_v0_Prescl;
   Int_t           L1Tech_TOTEM_Diffractive_v0;
   Int_t           L1Tech_TOTEM_Diffractive_v0_Prescl;
   Int_t           L1Tech_TOTEM_LowMultiplicity_v0;
   Int_t           L1Tech_TOTEM_LowMultiplicity_v0_Prescl;
   Int_t           L1Tech_TOTEM_MinBias_v0;
   Int_t           L1Tech_TOTEM_MinBias_v0_Prescl;
   Int_t           L1Tech_TOTEM_ZeroBias_v0;
   Int_t           L1Tech_TOTEM_ZeroBias_v0_Prescl;
   Int_t           L1Tech_ZDC_Scint_loose_vertex_v0;
   Int_t           L1Tech_ZDC_Scint_loose_vertex_v0_Prescl;
   Int_t           L1Tech_ZDC_Scint_minus_v0;
   Int_t           L1Tech_ZDC_Scint_minus_v0_Prescl;
   Int_t           L1Tech_ZDC_Scint_plus_v0;
   Int_t           L1Tech_ZDC_Scint_plus_v0_Prescl;
   Int_t           L1Tech_ZDC_Scint_tight_vertex_v0;
   Int_t           L1Tech_ZDC_Scint_tight_vertex_v0_Prescl;

   // List of branches
   TBranch        *b_NL1IsolEm;   //!
   TBranch        *b_L1IsolEmEt;   //!
   TBranch        *b_L1IsolEmE;   //!
   TBranch        *b_L1IsolEmEta;   //!
   TBranch        *b_L1IsolEmPhi;   //!
   TBranch        *b_NL1NIsolEm;   //!
   TBranch        *b_L1NIsolEmEt;   //!
   TBranch        *b_L1NIsolEmE;   //!
   TBranch        *b_L1NIsolEmEta;   //!
   TBranch        *b_L1NIsolEmPhi;   //!
   TBranch        *b_NL1Mu;   //!
   TBranch        *b_L1MuPt;   //!
   TBranch        *b_L1MuE;   //!
   TBranch        *b_L1MuEta;   //!
   TBranch        *b_L1MuPhi;   //!
   TBranch        *b_L1MuIsol;   //!
   TBranch        *b_L1MuMip;   //!
   TBranch        *b_L1MuFor;   //!
   TBranch        *b_L1MuRPC;   //!
   TBranch        *b_L1MuQal;   //!
   TBranch        *b_L1MuChg;   //!
   TBranch        *b_NL1CenJet;   //!
   TBranch        *b_L1CenJetEt;   //!
   TBranch        *b_L1CenJetE;   //!
   TBranch        *b_L1CenJetEta;   //!
   TBranch        *b_L1CenJetPhi;   //!
   TBranch        *b_NL1ForJet;   //!
   TBranch        *b_L1ForJetEt;   //!
   TBranch        *b_L1ForJetE;   //!
   TBranch        *b_L1ForJetEta;   //!
   TBranch        *b_L1ForJetPhi;   //!
   TBranch        *b_NL1Tau;   //!
   TBranch        *b_L1TauEt;   //!
   TBranch        *b_L1TauE;   //!
   TBranch        *b_L1TauEta;   //!
   TBranch        *b_L1TauPhi;   //!
   TBranch        *b_L1Met;   //!
   TBranch        *b_L1MetPhi;   //!
   TBranch        *b_L1EtTot;   //!
   TBranch        *b_L1Mht;   //!
   TBranch        *b_L1MhtPhi;   //!
   TBranch        *b_L1EtHad;   //!
   TBranch        *b_L1HfRing1EtSumPositiveEta;   //!
   TBranch        *b_L1HfRing2EtSumPositiveEta;   //!
   TBranch        *b_L1HfRing1EtSumNegativeEta;   //!
   TBranch        *b_L1HfRing2EtSumNegativeEta;   //!
   TBranch        *b_L1HfTowerCountPositiveEtaRing1;   //!
   TBranch        *b_L1HfTowerCountNegativeEtaRing1;   //!
   TBranch        *b_L1HfTowerCountPositiveEtaRing2;   //!
   TBranch        *b_L1HfTowerCountNegativeEtaRing2;   //!
   TBranch        *b_Run;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_LumiBlock;   //!
   TBranch        *b_Bx;   //!
   TBranch        *b_Orbit;   //!
   TBranch        *b_AvgInstDelLumi;   //!
   TBranch        *b_HLTriggerFirstPath;   //!
   TBranch        *b_HLTriggerFirstPath_Prescl;   //!
   TBranch        *b_HLT_Activity_Ecal_SC7_v13;   //!
   TBranch        *b_HLT_Activity_Ecal_SC7_v13_Prescl;   //!
   TBranch        *b_HLT_BeamGas_HF_Beam1_v5;   //!
   TBranch        *b_HLT_BeamGas_HF_Beam1_v5_Prescl;   //!
   TBranch        *b_HLT_BeamGas_HF_Beam2_v5;   //!
   TBranch        *b_HLT_BeamGas_HF_Beam2_v5_Prescl;   //!
   TBranch        *b_HLT_BeamHalo_v13;   //!
   TBranch        *b_HLT_BeamHalo_v13_Prescl;   //!
   TBranch        *b_HLT_PAHcalUTCA_v1;   //!
   TBranch        *b_HLT_PAHcalUTCA_v1_Prescl;   //!
   TBranch        *b_HLT_PAHcalPhiSym_v1;   //!
   TBranch        *b_HLT_PAHcalPhiSym_v1_Prescl;   //!
   TBranch        *b_HLT_PAHcalNZS_v1;   //!
   TBranch        *b_HLT_PAHcalNZS_v1_Prescl;   //!
   TBranch        *b_HLT_GlobalRunHPDNoise_v8;   //!
   TBranch        *b_HLT_GlobalRunHPDNoise_v8_Prescl;   //!
   TBranch        *b_HLT_Physics_v5;   //!
   TBranch        *b_HLT_Physics_v5_Prescl;   //!
   TBranch        *b_DST_Physics_v5;   //!
   TBranch        *b_DST_Physics_v5_Prescl;   //!
   TBranch        *b_HLT_DTCalibration_v2;   //!
   TBranch        *b_HLT_DTCalibration_v2_Prescl;   //!
   TBranch        *b_HLT_EcalCalibration_v3;   //!
   TBranch        *b_HLT_EcalCalibration_v3_Prescl;   //!
   TBranch        *b_HLT_HcalCalibration_v3;   //!
   TBranch        *b_HLT_HcalCalibration_v3_Prescl;   //!
   TBranch        *b_HLT_TrackerCalibration_v3;   //!
   TBranch        *b_HLT_TrackerCalibration_v3_Prescl;   //!
   TBranch        *b_HLT_L1SingleMuOpen_AntiBPTX_v7;   //!
   TBranch        *b_HLT_L1SingleMuOpen_AntiBPTX_v7_Prescl;   //!
   TBranch        *b_HLT_L1TrackerCosmics_v7;   //!
   TBranch        *b_HLT_L1TrackerCosmics_v7_Prescl;   //!
   TBranch        *b_AlCa_PAEcalPi0EBonly_v1;   //!
   TBranch        *b_AlCa_PAEcalPi0EBonly_v1_Prescl;   //!
   TBranch        *b_AlCa_PAEcalPi0EEonly_v1;   //!
   TBranch        *b_AlCa_PAEcalPi0EEonly_v1_Prescl;   //!
   TBranch        *b_AlCa_PAEcalEtaEBonly_v1;   //!
   TBranch        *b_AlCa_PAEcalEtaEBonly_v1_Prescl;   //!
   TBranch        *b_AlCa_PAEcalEtaEEonly_v1;   //!
   TBranch        *b_AlCa_PAEcalEtaEEonly_v1_Prescl;   //!
   TBranch        *b_AlCa_EcalPhiSym_v13;   //!
   TBranch        *b_AlCa_EcalPhiSym_v13_Prescl;   //!
   TBranch        *b_AlCa_RPCMuonNoTriggers_v9;   //!
   TBranch        *b_AlCa_RPCMuonNoTriggers_v9_Prescl;   //!
   TBranch        *b_AlCa_RPCMuonNoHits_v9;   //!
   TBranch        *b_AlCa_RPCMuonNoHits_v9_Prescl;   //!
   TBranch        *b_AlCa_RPCMuonNormalisation_v9;   //!
   TBranch        *b_AlCa_RPCMuonNormalisation_v9_Prescl;   //!
   TBranch        *b_AlCa_LumiPixels_v8;   //!
   TBranch        *b_AlCa_LumiPixels_v8_Prescl;   //!
   TBranch        *b_AlCa_LumiPixels_ZeroBias_v4;   //!
   TBranch        *b_AlCa_LumiPixels_ZeroBias_v4_Prescl;   //!
   TBranch        *b_AlCa_LumiPixels_Random_v1;   //!
   TBranch        *b_AlCa_LumiPixels_Random_v1_Prescl;   //!
   TBranch        *b_HLT_PAL1SingleJet16_v1;   //!
   TBranch        *b_HLT_PAL1SingleJet16_v1_Prescl;   //!
   TBranch        *b_HLT_PAL1SingleJet36_v1;   //!
   TBranch        *b_HLT_PAL1SingleJet36_v1_Prescl;   //!
   TBranch        *b_HLT_PASingleForJet15_v1;   //!
   TBranch        *b_HLT_PASingleForJet15_v1_Prescl;   //!
   TBranch        *b_HLT_PASingleForJet25_v1;   //!
   TBranch        *b_HLT_PASingleForJet25_v1_Prescl;   //!
   TBranch        *b_HLT_PAJet20_NoJetID_v1;   //!
   TBranch        *b_HLT_PAJet20_NoJetID_v1_Prescl;   //!
   TBranch        *b_HLT_PAJet40_NoJetID_v1;   //!
   TBranch        *b_HLT_PAJet40_NoJetID_v1_Prescl;   //!
   TBranch        *b_HLT_PAJet60_NoJetID_v1;   //!
   TBranch        *b_HLT_PAJet60_NoJetID_v1_Prescl;   //!
   TBranch        *b_HLT_PAJet80_NoJetID_v1;   //!
   TBranch        *b_HLT_PAJet80_NoJetID_v1_Prescl;   //!
   TBranch        *b_HLT_PAJet100_NoJetID_v1;   //!
   TBranch        *b_HLT_PAJet100_NoJetID_v1_Prescl;   //!
   TBranch        *b_HLT_PAJet120_NoJetID_v1;   //!
   TBranch        *b_HLT_PAJet120_NoJetID_v1_Prescl;   //!
   TBranch        *b_HLT_PAForJet20Eta2_v1;   //!
   TBranch        *b_HLT_PAForJet20Eta2_v1_Prescl;   //!
   TBranch        *b_HLT_PAForJet40Eta2_v1;   //!
   TBranch        *b_HLT_PAForJet40Eta2_v1_Prescl;   //!
   TBranch        *b_HLT_PAForJet60Eta2_v1;   //!
   TBranch        *b_HLT_PAForJet60Eta2_v1_Prescl;   //!
   TBranch        *b_HLT_PAForJet80Eta2_v1;   //!
   TBranch        *b_HLT_PAForJet80Eta2_v1_Prescl;   //!
   TBranch        *b_HLT_PAForJet100Eta2_v1;   //!
   TBranch        *b_HLT_PAForJet100Eta2_v1_Prescl;   //!
   TBranch        *b_HLT_PAForJet20Eta3_v1;   //!
   TBranch        *b_HLT_PAForJet20Eta3_v1_Prescl;   //!
   TBranch        *b_HLT_PAForJet40Eta3_v1;   //!
   TBranch        *b_HLT_PAForJet40Eta3_v1_Prescl;   //!
   TBranch        *b_HLT_PAForJet60Eta3_v1;   //!
   TBranch        *b_HLT_PAForJet60Eta3_v1_Prescl;   //!
   TBranch        *b_HLT_PAForJet80Eta3_v1;   //!
   TBranch        *b_HLT_PAForJet80Eta3_v1_Prescl;   //!
   TBranch        *b_HLT_PAForJet100Eta3_v1;   //!
   TBranch        *b_HLT_PAForJet100Eta3_v1_Prescl;   //!
   TBranch        *b_HLT_PATripleJet20_20_20_v1;   //!
   TBranch        *b_HLT_PATripleJet20_20_20_v1_Prescl;   //!
   TBranch        *b_HLT_PATripleJet40_20_20_v1;   //!
   TBranch        *b_HLT_PATripleJet40_20_20_v1_Prescl;   //!
   TBranch        *b_HLT_PATripleJet60_20_20_v1;   //!
   TBranch        *b_HLT_PATripleJet60_20_20_v1_Prescl;   //!
   TBranch        *b_HLT_PATripleJet80_20_20_v1;   //!
   TBranch        *b_HLT_PATripleJet80_20_20_v1_Prescl;   //!
   TBranch        *b_HLT_PATripleJet100_20_20_v1;   //!
   TBranch        *b_HLT_PATripleJet100_20_20_v1_Prescl;   //!
   TBranch        *b_HLT_PAJet40ETM30_v1;   //!
   TBranch        *b_HLT_PAJet40ETM30_v1_Prescl;   //!
   TBranch        *b_HLT_PAJet60ETM30_v1;   //!
   TBranch        *b_HLT_PAJet60ETM30_v1_Prescl;   //!
   TBranch        *b_HLT_PAL1DoubleMu0_v1;   //!
   TBranch        *b_HLT_PAL1DoubleMu0_v1_Prescl;   //!
   TBranch        *b_HLT_PADimuon0_NoVertexing_v1;   //!
   TBranch        *b_HLT_PADimuon0_NoVertexing_v1_Prescl;   //!
   TBranch        *b_HLT_PAL1DoubleMu0_HighQ_v1;   //!
   TBranch        *b_HLT_PAL1DoubleMu0_HighQ_v1_Prescl;   //!
   TBranch        *b_HLT_PAL1DoubleMuOpen_v1;   //!
   TBranch        *b_HLT_PAL1DoubleMuOpen_v1_Prescl;   //!
   TBranch        *b_HLT_PAL2DoubleMu3_v1;   //!
   TBranch        *b_HLT_PAL2DoubleMu3_v1_Prescl;   //!
   TBranch        *b_HLT_PAMu3_v1;   //!
   TBranch        *b_HLT_PAMu3_v1_Prescl;   //!
   TBranch        *b_HLT_PAMu7_v1;   //!
   TBranch        *b_HLT_PAMu7_v1_Prescl;   //!
   TBranch        *b_HLT_PAMu12_v1;   //!
   TBranch        *b_HLT_PAMu12_v1_Prescl;   //!
   TBranch        *b_HLT_Mu15_eta2p1_v5;   //!
   TBranch        *b_HLT_Mu15_eta2p1_v5_Prescl;   //!
   TBranch        *b_HLT_PABTagMu_Jet20_Mu4_v1;   //!
   TBranch        *b_HLT_PABTagMu_Jet20_Mu4_v1_Prescl;   //!
   TBranch        *b_HLT_PAMu3PFJet20_v1;   //!
   TBranch        *b_HLT_PAMu3PFJet20_v1_Prescl;   //!
   TBranch        *b_HLT_PAMu3PFJet40_v1;   //!
   TBranch        *b_HLT_PAMu3PFJet40_v1_Prescl;   //!
   TBranch        *b_HLT_PAMu7PFJet20_v1;   //!
   TBranch        *b_HLT_PAMu7PFJet20_v1_Prescl;   //!
   TBranch        *b_HLT_PAPhoton10_NoCaloIdVL_v1;   //!
   TBranch        *b_HLT_PAPhoton10_NoCaloIdVL_v1_Prescl;   //!
   TBranch        *b_HLT_PAPhoton15_NoCaloIdVL_v1;   //!
   TBranch        *b_HLT_PAPhoton15_NoCaloIdVL_v1_Prescl;   //!
   TBranch        *b_HLT_PAPhoton20_NoCaloIdVL_v1;   //!
   TBranch        *b_HLT_PAPhoton20_NoCaloIdVL_v1_Prescl;   //!
   TBranch        *b_HLT_PAPhoton30_NoCaloIdVL_v1;   //!
   TBranch        *b_HLT_PAPhoton30_NoCaloIdVL_v1_Prescl;   //!
   TBranch        *b_HLT_PAPhoton40_NoCaloIdVL_v1;   //!
   TBranch        *b_HLT_PAPhoton40_NoCaloIdVL_v1_Prescl;   //!
   TBranch        *b_HLT_PAPhoton60_NoCaloIdVL_v1;   //!
   TBranch        *b_HLT_PAPhoton60_NoCaloIdVL_v1_Prescl;   //!
   TBranch        *b_HLT_Ele22_CaloIdL_CaloIsoVL_v6;   //!
   TBranch        *b_HLT_Ele22_CaloIdL_CaloIsoVL_v6_Prescl;   //!
   TBranch        *b_HLT_PAPhoton10_TightCaloIdVL_v1;   //!
   TBranch        *b_HLT_PAPhoton10_TightCaloIdVL_v1_Prescl;   //!
   TBranch        *b_HLT_PAPhoton15_TightCaloIdVL_v1;   //!
   TBranch        *b_HLT_PAPhoton15_TightCaloIdVL_v1_Prescl;   //!
   TBranch        *b_HLT_PAPhoton20_TightCaloIdVL_v1;   //!
   TBranch        *b_HLT_PAPhoton20_TightCaloIdVL_v1_Prescl;   //!
   TBranch        *b_HLT_PAPhoton30_TightCaloIdVL_v1;   //!
   TBranch        *b_HLT_PAPhoton30_TightCaloIdVL_v1_Prescl;   //!
   TBranch        *b_HLT_PAPhoton40_TightCaloIdVL_v1;   //!
   TBranch        *b_HLT_PAPhoton40_TightCaloIdVL_v1_Prescl;   //!
   TBranch        *b_HLT_PAPhoton10_Photon10_NoCaloIdVL_v1;   //!
   TBranch        *b_HLT_PAPhoton10_Photon10_NoCaloIdVL_v1_Prescl;   //!
   TBranch        *b_HLT_PAPhoton15_Photon10_NoCaloIdVL_v1;   //!
   TBranch        *b_HLT_PAPhoton15_Photon10_NoCaloIdVL_v1_Prescl;   //!
   TBranch        *b_HLT_PAPhoton20_Photon15_NoCaloIdVL_v1;   //!
   TBranch        *b_HLT_PAPhoton20_Photon15_NoCaloIdVL_v1_Prescl;   //!
   TBranch        *b_HLT_PAPhoton20_Photon20_NoCaloIdVL_v1;   //!
   TBranch        *b_HLT_PAPhoton20_Photon20_NoCaloIdVL_v1_Prescl;   //!
   TBranch        *b_HLT_PAPhoton30_Photon30_NoCaloIdVL_v1;   //!
   TBranch        *b_HLT_PAPhoton30_Photon30_NoCaloIdVL_v1_Prescl;   //!
   TBranch        *b_HLT_PAPhoton10_Photon10_TightCaloIdVL_v1;   //!
   TBranch        *b_HLT_PAPhoton10_Photon10_TightCaloIdVL_v1_Prescl;   //!
   TBranch        *b_HLT_PAPhoton15_Photon10_TightCaloIdVL_v1;   //!
   TBranch        *b_HLT_PAPhoton15_Photon10_TightCaloIdVL_v1_Prescl;   //!
   TBranch        *b_HLT_PAPhoton20_Photon15_TightCaloIdVL_v1;   //!
   TBranch        *b_HLT_PAPhoton20_Photon15_TightCaloIdVL_v1_Prescl;   //!
   TBranch        *b_HLT_PAPixelTracks_Multiplicity100_v2;   //!
   TBranch        *b_HLT_PAPixelTracks_Multiplicity100_v2_Prescl;   //!
   TBranch        *b_HLT_PAPixelTracks_Multiplicity130_v2;   //!
   TBranch        *b_HLT_PAPixelTracks_Multiplicity130_v2_Prescl;   //!
   TBranch        *b_HLT_PAPixelTracks_Multiplicity160_v2;   //!
   TBranch        *b_HLT_PAPixelTracks_Multiplicity160_v2_Prescl;   //!
   TBranch        *b_HLT_PAPixelTracks_Multiplicity190_v2;   //!
   TBranch        *b_HLT_PAPixelTracks_Multiplicity190_v2_Prescl;   //!
   TBranch        *b_HLT_PAPixelTracks_Multiplicity220_v2;   //!
   TBranch        *b_HLT_PAPixelTracks_Multiplicity220_v2_Prescl;   //!
   TBranch        *b_HLT_PAPixelTrackMultiplicity100_FullTrack12_v2;   //!
   TBranch        *b_HLT_PAPixelTrackMultiplicity100_FullTrack12_v2_Prescl;   //!
   TBranch        *b_HLT_PAPixelTrackMultiplicity130_FullTrack12_v2;   //!
   TBranch        *b_HLT_PAPixelTrackMultiplicity130_FullTrack12_v2_Prescl;   //!
   TBranch        *b_HLT_PAPixelTrackMultiplicity160_FullTrack12_v2;   //!
   TBranch        *b_HLT_PAPixelTrackMultiplicity160_FullTrack12_v2_Prescl;   //!
   TBranch        *b_HLT_PAFullTrack12_v2;   //!
   TBranch        *b_HLT_PAFullTrack12_v2_Prescl;   //!
   TBranch        *b_HLT_PAFullTrack20_v2;   //!
   TBranch        *b_HLT_PAFullTrack20_v2_Prescl;   //!
   TBranch        *b_HLT_PAFullTrack30_v2;   //!
   TBranch        *b_HLT_PAFullTrack30_v2_Prescl;   //!
   TBranch        *b_HLT_PAFullTrack50_v2;   //!
   TBranch        *b_HLT_PAFullTrack50_v2_Prescl;   //!
   TBranch        *b_HLT_PAPixelTrackMultiplicity140_Jet80_NoJetID_v2;   //!
   TBranch        *b_HLT_PAPixelTrackMultiplicity140_Jet80_NoJetID_v2_Prescl;   //!
   TBranch        *b_HLT_PAPixelTrackMultiplicity100_L2DoubleMu3_v1;   //!
   TBranch        *b_HLT_PAPixelTrackMultiplicity100_L2DoubleMu3_v1_Prescl;   //!
   TBranch        *b_HLT_PPPixelTracks_Multiplicity55_v1;   //!
   TBranch        *b_HLT_PPPixelTracks_Multiplicity55_v1_Prescl;   //!
   TBranch        *b_HLT_PPPixelTracks_Multiplicity70_v1;   //!
   TBranch        *b_HLT_PPPixelTracks_Multiplicity70_v1_Prescl;   //!
   TBranch        *b_HLT_PPPixelTracks_Multiplicity85_v1;   //!
   TBranch        *b_HLT_PPPixelTracks_Multiplicity85_v1_Prescl;   //!
   TBranch        *b_HLT_PPPixelTrackMultiplicity55_FullTrack12_v1;   //!
   TBranch        *b_HLT_PPPixelTrackMultiplicity55_FullTrack12_v1_Prescl;   //!
   TBranch        *b_HLT_PPPixelTrackMultiplicity70_FullTrack12_v1;   //!
   TBranch        *b_HLT_PPPixelTrackMultiplicity70_FullTrack12_v1_Prescl;   //!
   TBranch        *b_HLT_PPL1DoubleJetC36_v1;   //!
   TBranch        *b_HLT_PPL1DoubleJetC36_v1_Prescl;   //!
   TBranch        *b_HLT_PATech35_v1;   //!
   TBranch        *b_HLT_PATech35_v1_Prescl;   //!
   TBranch        *b_HLT_PATech35_HFSumET100_v3;   //!
   TBranch        *b_HLT_PATech35_HFSumET100_v3_Prescl;   //!
   TBranch        *b_HLT_PAHFSumET100_v3;   //!
   TBranch        *b_HLT_PAHFSumET100_v3_Prescl;   //!
   TBranch        *b_HLT_PAHFSumET140_v3;   //!
   TBranch        *b_HLT_PAHFSumET140_v3_Prescl;   //!
   TBranch        *b_HLT_PAHFSumET170_v3;   //!
   TBranch        *b_HLT_PAHFSumET170_v3_Prescl;   //!
   TBranch        *b_HLT_PAHFSumET210_v3;   //!
   TBranch        *b_HLT_PAHFSumET210_v3_Prescl;   //!
   TBranch        *b_HLT_PARomanPots_Tech52_v1;   //!
   TBranch        *b_HLT_PARomanPots_Tech52_v1_Prescl;   //!
   TBranch        *b_HLT_PAL1Tech53_MB_v1;   //!
   TBranch        *b_HLT_PAL1Tech53_MB_v1_Prescl;   //!
   TBranch        *b_HLT_PAL1Tech53_MB_SingleTrack_v1;   //!
   TBranch        *b_HLT_PAL1Tech53_MB_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_PAL1Tech54_ZeroBias_v1;   //!
   TBranch        *b_HLT_PAL1Tech54_ZeroBias_v1_Prescl;   //!
   TBranch        *b_HLT_PAT1minbias_Tech55_v1;   //!
   TBranch        *b_HLT_PAT1minbias_Tech55_v1_Prescl;   //!
   TBranch        *b_HLT_PAL1Tech_HBHEHO_totalOR_v1;   //!
   TBranch        *b_HLT_PAL1Tech_HBHEHO_totalOR_v1_Prescl;   //!
   TBranch        *b_HLT_PAL1Tech63_CASTORHaloMuon_v1;   //!
   TBranch        *b_HLT_PAL1Tech63_CASTORHaloMuon_v1_Prescl;   //!
   TBranch        *b_HLT_PACastorEmTotemLowMultiplicity_v1;   //!
   TBranch        *b_HLT_PACastorEmTotemLowMultiplicity_v1_Prescl;   //!
   TBranch        *b_HLT_PACastorEmNotHfCoincidencePm_v1;   //!
   TBranch        *b_HLT_PACastorEmNotHfCoincidencePm_v1_Prescl;   //!
   TBranch        *b_HLT_PACastorEmNotHfSingleChannel_v1;   //!
   TBranch        *b_HLT_PACastorEmNotHfSingleChannel_v1_Prescl;   //!
   TBranch        *b_HLT_PAL1CastorTotalTotemLowMultiplicity_v1;   //!
   TBranch        *b_HLT_PAL1CastorTotalTotemLowMultiplicity_v1_Prescl;   //!
   TBranch        *b_HLT_PAMinBiasHF_v1;   //!
   TBranch        *b_HLT_PAMinBiasHF_v1_Prescl;   //!
   TBranch        *b_HLT_PAMinBiasHF_OR_v1;   //!
   TBranch        *b_HLT_PAMinBiasHF_OR_v1_Prescl;   //!
   TBranch        *b_HLT_PAMinBiasBHC_v1;   //!
   TBranch        *b_HLT_PAMinBiasBHC_v1_Prescl;   //!
   TBranch        *b_HLT_PAMinBiasBHC_OR_v1;   //!
   TBranch        *b_HLT_PAMinBiasBHC_OR_v1_Prescl;   //!
   TBranch        *b_HLT_PAMinBiasHfOrBHC_v1;   //!
   TBranch        *b_HLT_PAMinBiasHfOrBHC_v1_Prescl;   //!
   TBranch        *b_HLT_PABptxPlusNotBptxMinus_v1;   //!
   TBranch        *b_HLT_PABptxPlusNotBptxMinus_v1_Prescl;   //!
   TBranch        *b_HLT_PABptxMinusNotBptxPlus_v1;   //!
   TBranch        *b_HLT_PABptxMinusNotBptxPlus_v1_Prescl;   //!
   TBranch        *b_HLT_PAZeroBias_v1;   //!
   TBranch        *b_HLT_PAZeroBias_v1_Prescl;   //!
   TBranch        *b_HLT_PAZeroBiasPixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_PAHFOR_SingleTrack_v1;   //!
   TBranch        *b_HLT_PAHFOR_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_PAZeroBiasPixel_DoubleTrack_v1;   //!
   TBranch        *b_HLT_PAZeroBiasPixel_DoubleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_PADoubleMu4_Acoplanarity03_v1;   //!
   TBranch        *b_HLT_PADoubleMu4_Acoplanarity03_v1_Prescl;   //!
   TBranch        *b_HLT_PAExclDijet35_HFOR_v1;   //!
   TBranch        *b_HLT_PAExclDijet35_HFOR_v1_Prescl;   //!
   TBranch        *b_HLT_PAExclDijet35_HFAND_v1;   //!
   TBranch        *b_HLT_PAExclDijet35_HFAND_v1_Prescl;   //!
   TBranch        *b_HLT_PAL1DoubleEG3_FwdVeto_v1;   //!
   TBranch        *b_HLT_PAL1DoubleEG3_FwdVeto_v1_Prescl;   //!
   TBranch        *b_HLT_PAL1SingleJet52_TotemDiffractive_v1;   //!
   TBranch        *b_HLT_PAL1SingleJet52_TotemDiffractive_v1_Prescl;   //!
   TBranch        *b_HLT_PAL1SingleMu20_TotemDiffractive_v1;   //!
   TBranch        *b_HLT_PAL1SingleMu20_TotemDiffractive_v1_Prescl;   //!
   TBranch        *b_HLT_PAL1SingleEG20_TotemDiffractive_v1;   //!
   TBranch        *b_HLT_PAL1SingleEG20_TotemDiffractive_v1_Prescl;   //!
   TBranch        *b_HLT_PAL1DoubleJet20_TotemDiffractive_v1;   //!
   TBranch        *b_HLT_PAL1DoubleJet20_TotemDiffractive_v1_Prescl;   //!
   TBranch        *b_HLT_PAL1DoubleJetC36_TotemDiffractive_v1;   //!
   TBranch        *b_HLT_PAL1DoubleJetC36_TotemDiffractive_v1_Prescl;   //!
   TBranch        *b_HLT_PAL1DoubleMu5_TotemDiffractive_v1;   //!
   TBranch        *b_HLT_PAL1DoubleMu5_TotemDiffractive_v1_Prescl;   //!
   TBranch        *b_HLT_PAL1DoubleEG5_TotemDiffractive_v1;   //!
   TBranch        *b_HLT_PAL1DoubleEG5_TotemDiffractive_v1_Prescl;   //!
   TBranch        *b_HLT_PADoubleJet20_ForwardBackward_v1;   //!
   TBranch        *b_HLT_PADoubleJet20_ForwardBackward_v1_Prescl;   //!
   TBranch        *b_HLT_PAMu7_Ele7_CaloIdT_CaloIsoVL_v1;   //!
   TBranch        *b_HLT_PAMu7_Ele7_CaloIdT_CaloIsoVL_v1_Prescl;   //!
   TBranch        *b_HLT_PARandom_v1;   //!
   TBranch        *b_HLT_PARandom_v1_Prescl;   //!
   TBranch        *b_DQM_FEDIntegrity_v11;   //!
   TBranch        *b_DQM_FEDIntegrity_v11_Prescl;   //!
   TBranch        *b_HLT_LogMonitor_v4;   //!
   TBranch        *b_HLT_LogMonitor_v4_Prescl;   //!
   TBranch        *b_HLTriggerFinalPath;   //!
   TBranch        *b_HLTriggerFinalPath_Prescl;   //!
   TBranch        *b_L1_AlwaysTrue;   //!
   TBranch        *b_L1_AlwaysTrue_Prescl;   //!
   TBranch        *b_L1_BeamGas_Hf_BptxMinusPostQuiet;   //!
   TBranch        *b_L1_BeamGas_Hf_BptxMinusPostQuiet_Prescl;   //!
   TBranch        *b_L1_BeamGas_Hf_BptxPlusPostQuiet;   //!
   TBranch        *b_L1_BeamGas_Hf_BptxPlusPostQuiet_Prescl;   //!
   TBranch        *b_L1_BeamHalo;   //!
   TBranch        *b_L1_BeamHalo_Prescl;   //!
   TBranch        *b_L1_BptxMinus_NotBptxPlus;   //!
   TBranch        *b_L1_BptxMinus_NotBptxPlus_Prescl;   //!
   TBranch        *b_L1_BptxPlus_NotBptxMinus;   //!
   TBranch        *b_L1_BptxPlus_NotBptxMinus_Prescl;   //!
   TBranch        *b_L1_BscHighMultiplicity_BptxAND;   //!
   TBranch        *b_L1_BscHighMultiplicity_BptxAND_Prescl;   //!
   TBranch        *b_L1_BscMinBiasOR_BptxAND;   //!
   TBranch        *b_L1_BscMinBiasOR_BptxAND_Prescl;   //!
   TBranch        *b_L1_BscMinBiasThreshold1_BptxAND;   //!
   TBranch        *b_L1_BscMinBiasThreshold1_BptxAND_Prescl;   //!
   TBranch        *b_L1_BscMinBiasThreshold2_BptxAND;   //!
   TBranch        *b_L1_BscMinBiasThreshold2_BptxAND_Prescl;   //!
   TBranch        *b_L1_CastorEm_NotHcalHfCoincidencePm;   //!
   TBranch        *b_L1_CastorEm_NotHcalHfCoincidencePm_Prescl;   //!
   TBranch        *b_L1_CastorEm_NotHcalHfSingleChannel;   //!
   TBranch        *b_L1_CastorEm_NotHcalHfSingleChannel_Prescl;   //!
   TBranch        *b_L1_CastorEm_TotemLowMultiplicity;   //!
   TBranch        *b_L1_CastorEm_TotemLowMultiplicity_Prescl;   //!
   TBranch        *b_L1_CastorTotalEnergy_TotemLowMultiplicity;   //!
   TBranch        *b_L1_CastorTotalEnergy_TotemLowMultiplicity_Prescl;   //!
   TBranch        *b_L1_DoubleEG3_FwdVeto;   //!
   TBranch        *b_L1_DoubleEG3_FwdVeto_Prescl;   //!
   TBranch        *b_L1_DoubleEG5;   //!
   TBranch        *b_L1_DoubleEG5_Prescl;   //!
   TBranch        *b_L1_DoubleEG5_TotemDiffractive;   //!
   TBranch        *b_L1_DoubleEG5_TotemDiffractive_Prescl;   //!
   TBranch        *b_L1_DoubleEG6_HTT100;   //!
   TBranch        *b_L1_DoubleEG6_HTT100_Prescl;   //!
   TBranch        *b_L1_DoubleEG6_HTT125;   //!
   TBranch        *b_L1_DoubleEG6_HTT125_Prescl;   //!
   TBranch        *b_L1_DoubleEG_13_7;   //!
   TBranch        *b_L1_DoubleEG_13_7_Prescl;   //!
   TBranch        *b_L1_DoubleForJet16_EtaOpp;   //!
   TBranch        *b_L1_DoubleForJet16_EtaOpp_Prescl;   //!
   TBranch        *b_L1_DoubleJet20;   //!
   TBranch        *b_L1_DoubleJet20_Prescl;   //!
   TBranch        *b_L1_DoubleJet20_TotemDiffractive;   //!
   TBranch        *b_L1_DoubleJet20_TotemDiffractive_Prescl;   //!
   TBranch        *b_L1_DoubleJet24_v1;   //!
   TBranch        *b_L1_DoubleJet24_v1_Prescl;   //!
   TBranch        *b_L1_DoubleJet36_Central;   //!
   TBranch        *b_L1_DoubleJet36_Central_Prescl;   //!
   TBranch        *b_L1_DoubleJet52_Central;   //!
   TBranch        *b_L1_DoubleJet52_Central_Prescl;   //!
   TBranch        *b_L1_DoubleJetC36_TotemDiffractive;   //!
   TBranch        *b_L1_DoubleJetC36_TotemDiffractive_Prescl;   //!
   TBranch        *b_L1_DoubleJetC44_ETM30;   //!
   TBranch        *b_L1_DoubleJetC44_ETM30_Prescl;   //!
   TBranch        *b_L1_DoubleJetC56;   //!
   TBranch        *b_L1_DoubleJetC56_Prescl;   //!
   TBranch        *b_L1_DoubleJetC56_Eta1p74_WdEta4;   //!
   TBranch        *b_L1_DoubleJetC56_Eta1p74_WdEta4_Prescl;   //!
   TBranch        *b_L1_DoubleMu0;   //!
   TBranch        *b_L1_DoubleMu0_Prescl;   //!
   TBranch        *b_L1_DoubleMu0_HighQ_EtaCuts;   //!
   TBranch        *b_L1_DoubleMu0_HighQ_EtaCuts_Prescl;   //!
   TBranch        *b_L1_DoubleMu3p5_EG5;   //!
   TBranch        *b_L1_DoubleMu3p5_EG5_Prescl;   //!
   TBranch        *b_L1_DoubleMu5_EG5;   //!
   TBranch        *b_L1_DoubleMu5_EG5_Prescl;   //!
   TBranch        *b_L1_DoubleMu5_TotemDiffractive;   //!
   TBranch        *b_L1_DoubleMu5_TotemDiffractive_Prescl;   //!
   TBranch        *b_L1_DoubleMu5_v1;   //!
   TBranch        *b_L1_DoubleMu5_v1_Prescl;   //!
   TBranch        *b_L1_DoubleMuOpen_BptxAND;   //!
   TBranch        *b_L1_DoubleMuOpen_BptxAND_Prescl;   //!
   TBranch        *b_L1_DoubleMu_10_3p5;   //!
   TBranch        *b_L1_DoubleMu_10_3p5_Prescl;   //!
   TBranch        *b_L1_DoubleMu_10_Open;   //!
   TBranch        *b_L1_DoubleMu_10_Open_Prescl;   //!
   TBranch        *b_L1_DoubleMu_12_5;   //!
   TBranch        *b_L1_DoubleMu_12_5_Prescl;   //!
   TBranch        *b_L1_DoubleMu_3er_0er_HighQ_WdEta22;   //!
   TBranch        *b_L1_DoubleMu_3er_0er_HighQ_WdEta22_Prescl;   //!
   TBranch        *b_L1_DoubleMu_5er_0er_HighQ_WdEta22;   //!
   TBranch        *b_L1_DoubleMu_5er_0er_HighQ_WdEta22_Prescl;   //!
   TBranch        *b_L1_EG8_DoubleJetC20;   //!
   TBranch        *b_L1_EG8_DoubleJetC20_Prescl;   //!
   TBranch        *b_L1_ETM100;   //!
   TBranch        *b_L1_ETM100_Prescl;   //!
   TBranch        *b_L1_ETM30;   //!
   TBranch        *b_L1_ETM30_Prescl;   //!
   TBranch        *b_L1_ETM36;   //!
   TBranch        *b_L1_ETM36_Prescl;   //!
   TBranch        *b_L1_ETM40;   //!
   TBranch        *b_L1_ETM40_Prescl;   //!
   TBranch        *b_L1_ETM50;   //!
   TBranch        *b_L1_ETM50_Prescl;   //!
   TBranch        *b_L1_ETM70;   //!
   TBranch        *b_L1_ETM70_Prescl;   //!
   TBranch        *b_L1_ETT140;   //!
   TBranch        *b_L1_ETT140_Prescl;   //!
   TBranch        *b_L1_ETT20_BptxAND;   //!
   TBranch        *b_L1_ETT20_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT300;   //!
   TBranch        *b_L1_ETT300_Prescl;   //!
   TBranch        *b_L1_ETT40;   //!
   TBranch        *b_L1_ETT40_Prescl;   //!
   TBranch        *b_L1_ETT60;   //!
   TBranch        *b_L1_ETT60_Prescl;   //!
   TBranch        *b_L1_ETT80;   //!
   TBranch        *b_L1_ETT80_Prescl;   //!
   TBranch        *b_L1_HTT100;   //!
   TBranch        *b_L1_HTT100_Prescl;   //!
   TBranch        *b_L1_HTT125;   //!
   TBranch        *b_L1_HTT125_Prescl;   //!
   TBranch        *b_L1_HTT150;   //!
   TBranch        *b_L1_HTT150_Prescl;   //!
   TBranch        *b_L1_HTT175;   //!
   TBranch        *b_L1_HTT175_Prescl;   //!
   TBranch        *b_L1_HTT200;   //!
   TBranch        *b_L1_HTT200_Prescl;   //!
   TBranch        *b_L1_HTT50;   //!
   TBranch        *b_L1_HTT50_Prescl;   //!
   TBranch        *b_L1_HTT75;   //!
   TBranch        *b_L1_HTT75_Prescl;   //!
   TBranch        *b_L1_HcalHfCoincidencePm_BptxAND_v1;   //!
   TBranch        *b_L1_HcalHfCoincidencePm_BptxAND_v1_Prescl;   //!
   TBranch        *b_L1_HcalHfSingleChannel_BptxAND;   //!
   TBranch        *b_L1_HcalHfSingleChannel_BptxAND_Prescl;   //!
   TBranch        *b_L1_HcalHfSingleChannel_BptxAND_Instance1;   //!
   TBranch        *b_L1_HcalHfSingleChannel_BptxAND_Instance1_Prescl;   //!
   TBranch        *b_L1_Mu10er_JetC32;   //!
   TBranch        *b_L1_Mu10er_JetC32_Prescl;   //!
   TBranch        *b_L1_Mu12_EG7;   //!
   TBranch        *b_L1_Mu12_EG7_Prescl;   //!
   TBranch        *b_L1_Mu3_Jet16;   //!
   TBranch        *b_L1_Mu3_Jet16_Prescl;   //!
   TBranch        *b_L1_Mu3_Jet36;   //!
   TBranch        *b_L1_Mu3_Jet36_Prescl;   //!
   TBranch        *b_L1_Mu3_JetC16_WdEtaPhi2;   //!
   TBranch        *b_L1_Mu3_JetC16_WdEtaPhi2_Prescl;   //!
   TBranch        *b_L1_Mu3_JetC52_WdEtaPhi2;   //!
   TBranch        *b_L1_Mu3_JetC52_WdEtaPhi2_Prescl;   //!
   TBranch        *b_L1_Mu5_DoubleEG5;   //!
   TBranch        *b_L1_Mu5_DoubleEG5_Prescl;   //!
   TBranch        *b_L1_Mu5_DoubleEG6;   //!
   TBranch        *b_L1_Mu5_DoubleEG6_Prescl;   //!
   TBranch        *b_L1_Mu7_Jet16;   //!
   TBranch        *b_L1_Mu7_Jet16_Prescl;   //!
   TBranch        *b_L1_Mu8_DoubleJetC20;   //!
   TBranch        *b_L1_Mu8_DoubleJetC20_Prescl;   //!
   TBranch        *b_L1_MuOpen_EG12;   //!
   TBranch        *b_L1_MuOpen_EG12_Prescl;   //!
   TBranch        *b_L1_MuOpen_EG5;   //!
   TBranch        *b_L1_MuOpen_EG5_Prescl;   //!
   TBranch        *b_L1_QuadJetC32;   //!
   TBranch        *b_L1_QuadJetC32_Prescl;   //!
   TBranch        *b_L1_QuadJetC36;   //!
   TBranch        *b_L1_QuadJetC36_Prescl;   //!
   TBranch        *b_L1_QuadJetC40;   //!
   TBranch        *b_L1_QuadJetC40_Prescl;   //!
   TBranch        *b_L1_SingleEG12;   //!
   TBranch        *b_L1_SingleEG12_Prescl;   //!
   TBranch        *b_L1_SingleEG18er;   //!
   TBranch        *b_L1_SingleEG18er_Prescl;   //!
   TBranch        *b_L1_SingleEG20;   //!
   TBranch        *b_L1_SingleEG20_Prescl;   //!
   TBranch        *b_L1_SingleEG20_TotemDiffractive;   //!
   TBranch        *b_L1_SingleEG20_TotemDiffractive_Prescl;   //!
   TBranch        *b_L1_SingleEG22;   //!
   TBranch        *b_L1_SingleEG22_Prescl;   //!
   TBranch        *b_L1_SingleEG24;   //!
   TBranch        *b_L1_SingleEG24_Prescl;   //!
   TBranch        *b_L1_SingleEG30;   //!
   TBranch        *b_L1_SingleEG30_Prescl;   //!
   TBranch        *b_L1_SingleEG5_BptxAND;   //!
   TBranch        *b_L1_SingleEG5_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleEG7;   //!
   TBranch        *b_L1_SingleEG7_Prescl;   //!
   TBranch        *b_L1_SingleForJet16;   //!
   TBranch        *b_L1_SingleForJet16_Prescl;   //!
   TBranch        *b_L1_SingleIsoEG18er;   //!
   TBranch        *b_L1_SingleIsoEG18er_Prescl;   //!
   TBranch        *b_L1_SingleIsoEG20er;   //!
   TBranch        *b_L1_SingleIsoEG20er_Prescl;   //!
   TBranch        *b_L1_SingleJet128;   //!
   TBranch        *b_L1_SingleJet128_Prescl;   //!
   TBranch        *b_L1_SingleJet12_BptxAND;   //!
   TBranch        *b_L1_SingleJet12_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet16_BptxAND;   //!
   TBranch        *b_L1_SingleJet16_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet16_FwdVeto5;   //!
   TBranch        *b_L1_SingleJet16_FwdVeto5_Prescl;   //!
   TBranch        *b_L1_SingleJet20_Central_NotBptxOR;   //!
   TBranch        *b_L1_SingleJet20_Central_NotBptxOR_Prescl;   //!
   TBranch        *b_L1_SingleJet36;   //!
   TBranch        *b_L1_SingleJet36_Prescl;   //!
   TBranch        *b_L1_SingleJet36_FwdVeto5;   //!
   TBranch        *b_L1_SingleJet36_FwdVeto5_Prescl;   //!
   TBranch        *b_L1_SingleJet52;   //!
   TBranch        *b_L1_SingleJet52_Prescl;   //!
   TBranch        *b_L1_SingleJet52_TotemDiffractive;   //!
   TBranch        *b_L1_SingleJet52_TotemDiffractive_Prescl;   //!
   TBranch        *b_L1_SingleJet68;   //!
   TBranch        *b_L1_SingleJet68_Prescl;   //!
   TBranch        *b_L1_SingleJet92;   //!
   TBranch        *b_L1_SingleJet92_Prescl;   //!
   TBranch        *b_L1_SingleJetC32_NotBptxOR;   //!
   TBranch        *b_L1_SingleJetC32_NotBptxOR_Prescl;   //!
   TBranch        *b_L1_SingleMu12;   //!
   TBranch        *b_L1_SingleMu12_Prescl;   //!
   TBranch        *b_L1_SingleMu12er;   //!
   TBranch        *b_L1_SingleMu12er_Prescl;   //!
   TBranch        *b_L1_SingleMu14_Eta2p1;   //!
   TBranch        *b_L1_SingleMu14_Eta2p1_Prescl;   //!
   TBranch        *b_L1_SingleMu16;   //!
   TBranch        *b_L1_SingleMu16_Prescl;   //!
   TBranch        *b_L1_SingleMu16_Eta2p1;   //!
   TBranch        *b_L1_SingleMu16_Eta2p1_Prescl;   //!
   TBranch        *b_L1_SingleMu18er;   //!
   TBranch        *b_L1_SingleMu18er_Prescl;   //!
   TBranch        *b_L1_SingleMu20;   //!
   TBranch        *b_L1_SingleMu20_Prescl;   //!
   TBranch        *b_L1_SingleMu20_TotemDiffractive;   //!
   TBranch        *b_L1_SingleMu20_TotemDiffractive_Prescl;   //!
   TBranch        *b_L1_SingleMu20er;   //!
   TBranch        *b_L1_SingleMu20er_Prescl;   //!
   TBranch        *b_L1_SingleMu25er;   //!
   TBranch        *b_L1_SingleMu25er_Prescl;   //!
   TBranch        *b_L1_SingleMu3;   //!
   TBranch        *b_L1_SingleMu3_Prescl;   //!
   TBranch        *b_L1_SingleMu6_NotBptxOR;   //!
   TBranch        *b_L1_SingleMu6_NotBptxOR_Prescl;   //!
   TBranch        *b_L1_SingleMu7;   //!
   TBranch        *b_L1_SingleMu7_Prescl;   //!
   TBranch        *b_L1_SingleMuBeamHalo;   //!
   TBranch        *b_L1_SingleMuBeamHalo_Prescl;   //!
   TBranch        *b_L1_SingleMuOpen;   //!
   TBranch        *b_L1_SingleMuOpen_Prescl;   //!
   TBranch        *b_L1_TotemMinBias;   //!
   TBranch        *b_L1_TotemMinBias_Prescl;   //!
   TBranch        *b_L1_TripleEG7;   //!
   TBranch        *b_L1_TripleEG7_Prescl;   //!
   TBranch        *b_L1_TripleEG_12_7_5;   //!
   TBranch        *b_L1_TripleEG_12_7_5_Prescl;   //!
   TBranch        *b_L1_TripleJetC_52_28_28;   //!
   TBranch        *b_L1_TripleJetC_52_28_28_Prescl;   //!
   TBranch        *b_L1_TripleJet_64_44_24_VBF;   //!
   TBranch        *b_L1_TripleJet_64_44_24_VBF_Prescl;   //!
   TBranch        *b_L1_TripleJet_64_48_28_VBF;   //!
   TBranch        *b_L1_TripleJet_64_48_28_VBF_Prescl;   //!
   TBranch        *b_L1_ZdcCaloPlus_TotemDiffractive_QElastic;   //!
   TBranch        *b_L1_ZdcCaloPlus_TotemDiffractive_QElastic_Prescl;   //!
   TBranch        *b_L1_ZeroBias;   //!
   TBranch        *b_L1_ZeroBias_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_PreBPTX_v0;   //!
   TBranch        *b_L1Tech_BPTX_PreBPTX_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_minus_v0;   //!
   TBranch        *b_L1Tech_BPTX_minus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_minus_AND_not_plus_v0;   //!
   TBranch        *b_L1Tech_BPTX_minus_AND_not_plus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_NOT_minus_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_NOT_minus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_minus_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_minus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_minus_instance1_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_minus_instance1_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_OR_minus_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_OR_minus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_quiet_v0;   //!
   TBranch        *b_L1Tech_BPTX_quiet_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_HighMultiplicity_v0;   //!
   TBranch        *b_L1Tech_BSC_HighMultiplicity_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_halo_beam1_inner_v0;   //!
   TBranch        *b_L1Tech_BSC_halo_beam1_inner_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_halo_beam1_outer_v0;   //!
   TBranch        *b_L1Tech_BSC_halo_beam1_outer_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_halo_beam2_inner_v0;   //!
   TBranch        *b_L1Tech_BSC_halo_beam2_inner_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_halo_beam2_outer_v0;   //!
   TBranch        *b_L1Tech_BSC_halo_beam2_outer_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_minBias_OR_v0;   //!
   TBranch        *b_L1Tech_BSC_minBias_OR_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_minBias_inner_threshold1_v0;   //!
   TBranch        *b_L1Tech_BSC_minBias_inner_threshold1_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_minBias_inner_threshold2_v0;   //!
   TBranch        *b_L1Tech_BSC_minBias_inner_threshold2_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_minBias_threshold1_v0;   //!
   TBranch        *b_L1Tech_BSC_minBias_threshold1_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_minBias_threshold2_v0;   //!
   TBranch        *b_L1Tech_BSC_minBias_threshold2_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_splash_beam1_v0;   //!
   TBranch        *b_L1Tech_BSC_splash_beam1_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_splash_beam2_v0;   //!
   TBranch        *b_L1Tech_BSC_splash_beam2_v0_Prescl;   //!
   TBranch        *b_L1Tech_CASTOR_0_v0;   //!
   TBranch        *b_L1Tech_CASTOR_0_v0_Prescl;   //!
   TBranch        *b_L1Tech_CASTOR_EM_v0;   //!
   TBranch        *b_L1Tech_CASTOR_EM_v0_Prescl;   //!
   TBranch        *b_L1Tech_CASTOR_HaloMuon_v0;   //!
   TBranch        *b_L1Tech_CASTOR_HaloMuon_v0_Prescl;   //!
   TBranch        *b_L1Tech_CASTOR_TotalEnergy_v0;   //!
   TBranch        *b_L1Tech_CASTOR_TotalEnergy_v0_Prescl;   //!
   TBranch        *b_L1Tech_DT_GlobalOR_v0;   //!
   TBranch        *b_L1Tech_DT_GlobalOR_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St3Sect45_downLeft_v0;   //!
   TBranch        *b_L1Tech_FSC_St3Sect45_downLeft_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St3Sect45_downRight_v0;   //!
   TBranch        *b_L1Tech_FSC_St3Sect45_downRight_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St3Sect45_uppLeft_v0;   //!
   TBranch        *b_L1Tech_FSC_St3Sect45_uppLeft_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St3Sect45_uppRight_v0;   //!
   TBranch        *b_L1Tech_FSC_St3Sect45_uppRight_v0_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HBHE_totalOR_v0;   //!
   TBranch        *b_L1Tech_HCAL_HBHE_totalOR_v0_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HF_MMP_or_MPP_v1;   //!
   TBranch        *b_L1Tech_HCAL_HF_MMP_or_MPP_v1_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HF_coincidence_PM_v2;   //!
   TBranch        *b_L1Tech_HCAL_HF_coincidence_PM_v2_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HF_single_channel_v0;   //!
   TBranch        *b_L1Tech_HCAL_HF_single_channel_v0_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HO_totalOR_v0;   //!
   TBranch        *b_L1Tech_HCAL_HO_totalOR_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_RB0_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_RB0_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBminus1_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBminus1_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBminus2_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBminus2_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBplus1_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBplus1_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBplus2_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBplus2_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBst1_collisions_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBst1_collisions_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_barrel_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_barrel_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_pointing_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_pointing_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_TOTEM_Diffractive_v0;   //!
   TBranch        *b_L1Tech_TOTEM_Diffractive_v0_Prescl;   //!
   TBranch        *b_L1Tech_TOTEM_LowMultiplicity_v0;   //!
   TBranch        *b_L1Tech_TOTEM_LowMultiplicity_v0_Prescl;   //!
   TBranch        *b_L1Tech_TOTEM_MinBias_v0;   //!
   TBranch        *b_L1Tech_TOTEM_MinBias_v0_Prescl;   //!
   TBranch        *b_L1Tech_TOTEM_ZeroBias_v0;   //!
   TBranch        *b_L1Tech_TOTEM_ZeroBias_v0_Prescl;   //!
   TBranch        *b_L1Tech_ZDC_Scint_loose_vertex_v0;   //!
   TBranch        *b_L1Tech_ZDC_Scint_loose_vertex_v0_Prescl;   //!
   TBranch        *b_L1Tech_ZDC_Scint_minus_v0;   //!
   TBranch        *b_L1Tech_ZDC_Scint_minus_v0_Prescl;   //!
   TBranch        *b_L1Tech_ZDC_Scint_plus_v0;   //!
   TBranch        *b_L1Tech_ZDC_Scint_plus_v0_Prescl;   //!
   TBranch        *b_L1Tech_ZDC_Scint_tight_vertex_v0;   //!
   TBranch        *b_L1Tech_ZDC_Scint_tight_vertex_v0_Prescl;   //!

   HltTree_pp(TString infile,TTree *tree=0);
   virtual ~HltTree_pp();
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

#ifdef HltTree_pp_cxx

HltTree_pp::HltTree_pp(TString infile,TTree *tree) 
{
      f = TFile::Open(infile.Data());
      tree=(TTree*)f->Get("hltanalysis/HltTree");
	  Init(tree);
	  // f->Close();
}


HltTree_pp::~HltTree_pp()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}


Int_t HltTree_pp::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t HltTree_pp::LoadTree(Long64_t entry)
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

void HltTree_pp::Init(TTree *tree)
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

   fChain->SetBranchAddress("NL1IsolEm", &NL1IsolEm, &b_NL1IsolEm);
   fChain->SetBranchAddress("L1IsolEmEt", L1IsolEmEt, &b_L1IsolEmEt);
   fChain->SetBranchAddress("L1IsolEmE", L1IsolEmE, &b_L1IsolEmE);
   fChain->SetBranchAddress("L1IsolEmEta", L1IsolEmEta, &b_L1IsolEmEta);
   fChain->SetBranchAddress("L1IsolEmPhi", L1IsolEmPhi, &b_L1IsolEmPhi);
   fChain->SetBranchAddress("NL1NIsolEm", &NL1NIsolEm, &b_NL1NIsolEm);
   fChain->SetBranchAddress("L1NIsolEmEt", L1NIsolEmEt, &b_L1NIsolEmEt);
   fChain->SetBranchAddress("L1NIsolEmE", L1NIsolEmE, &b_L1NIsolEmE);
   fChain->SetBranchAddress("L1NIsolEmEta", L1NIsolEmEta, &b_L1NIsolEmEta);
   fChain->SetBranchAddress("L1NIsolEmPhi", L1NIsolEmPhi, &b_L1NIsolEmPhi);
   fChain->SetBranchAddress("NL1Mu", &NL1Mu, &b_NL1Mu);
   fChain->SetBranchAddress("L1MuPt", L1MuPt, &b_L1MuPt);
   fChain->SetBranchAddress("L1MuE", L1MuE, &b_L1MuE);
   fChain->SetBranchAddress("L1MuEta", L1MuEta, &b_L1MuEta);
   fChain->SetBranchAddress("L1MuPhi", L1MuPhi, &b_L1MuPhi);
   fChain->SetBranchAddress("L1MuIsol", L1MuIsol, &b_L1MuIsol);
   fChain->SetBranchAddress("L1MuMip", L1MuMip, &b_L1MuMip);
   fChain->SetBranchAddress("L1MuFor", L1MuFor, &b_L1MuFor);
   fChain->SetBranchAddress("L1MuRPC", L1MuRPC, &b_L1MuRPC);
   fChain->SetBranchAddress("L1MuQal", L1MuQal, &b_L1MuQal);
   fChain->SetBranchAddress("L1MuChg", L1MuChg, &b_L1MuChg);
   fChain->SetBranchAddress("NL1CenJet", &NL1CenJet, &b_NL1CenJet);
   fChain->SetBranchAddress("L1CenJetEt", L1CenJetEt, &b_L1CenJetEt);
   fChain->SetBranchAddress("L1CenJetE", L1CenJetE, &b_L1CenJetE);
   fChain->SetBranchAddress("L1CenJetEta", L1CenJetEta, &b_L1CenJetEta);
   fChain->SetBranchAddress("L1CenJetPhi", L1CenJetPhi, &b_L1CenJetPhi);
   fChain->SetBranchAddress("NL1ForJet", &NL1ForJet, &b_NL1ForJet);
   fChain->SetBranchAddress("L1ForJetEt", L1ForJetEt, &b_L1ForJetEt);
   fChain->SetBranchAddress("L1ForJetE", L1ForJetE, &b_L1ForJetE);
   fChain->SetBranchAddress("L1ForJetEta", L1ForJetEta, &b_L1ForJetEta);
   fChain->SetBranchAddress("L1ForJetPhi", L1ForJetPhi, &b_L1ForJetPhi);
   fChain->SetBranchAddress("NL1Tau", &NL1Tau, &b_NL1Tau);
   fChain->SetBranchAddress("L1TauEt", L1TauEt, &b_L1TauEt);
   fChain->SetBranchAddress("L1TauE", L1TauE, &b_L1TauE);
   fChain->SetBranchAddress("L1TauEta", L1TauEta, &b_L1TauEta);
   fChain->SetBranchAddress("L1TauPhi", L1TauPhi, &b_L1TauPhi);
   fChain->SetBranchAddress("L1Met", &L1Met, &b_L1Met);
   fChain->SetBranchAddress("L1MetPhi", &L1MetPhi, &b_L1MetPhi);
   fChain->SetBranchAddress("L1EtTot", &L1EtTot, &b_L1EtTot);
   fChain->SetBranchAddress("L1Mht", &L1Mht, &b_L1Mht);
   fChain->SetBranchAddress("L1MhtPhi", &L1MhtPhi, &b_L1MhtPhi);
   fChain->SetBranchAddress("L1EtHad", &L1EtHad, &b_L1EtHad);
   fChain->SetBranchAddress("L1HfRing1EtSumPositiveEta", &L1HfRing1EtSumPositiveEta, &b_L1HfRing1EtSumPositiveEta);
   fChain->SetBranchAddress("L1HfRing2EtSumPositiveEta", &L1HfRing2EtSumPositiveEta, &b_L1HfRing2EtSumPositiveEta);
   fChain->SetBranchAddress("L1HfRing1EtSumNegativeEta", &L1HfRing1EtSumNegativeEta, &b_L1HfRing1EtSumNegativeEta);
   fChain->SetBranchAddress("L1HfRing2EtSumNegativeEta", &L1HfRing2EtSumNegativeEta, &b_L1HfRing2EtSumNegativeEta);
   fChain->SetBranchAddress("L1HfTowerCountPositiveEtaRing1", &L1HfTowerCountPositiveEtaRing1, &b_L1HfTowerCountPositiveEtaRing1);
   fChain->SetBranchAddress("L1HfTowerCountNegativeEtaRing1", &L1HfTowerCountNegativeEtaRing1, &b_L1HfTowerCountNegativeEtaRing1);
   fChain->SetBranchAddress("L1HfTowerCountPositiveEtaRing2", &L1HfTowerCountPositiveEtaRing2, &b_L1HfTowerCountPositiveEtaRing2);
   fChain->SetBranchAddress("L1HfTowerCountNegativeEtaRing2", &L1HfTowerCountNegativeEtaRing2, &b_L1HfTowerCountNegativeEtaRing2);
   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("LumiBlock", &LumiBlock, &b_LumiBlock);
   fChain->SetBranchAddress("Bx", &Bx, &b_Bx);
   fChain->SetBranchAddress("Orbit", &Orbit, &b_Orbit);
   fChain->SetBranchAddress("AvgInstDelLumi", &AvgInstDelLumi, &b_AvgInstDelLumi);
   fChain->SetBranchAddress("HLTriggerFirstPath", &HLTriggerFirstPath, &b_HLTriggerFirstPath);
   fChain->SetBranchAddress("HLTriggerFirstPath_Prescl", &HLTriggerFirstPath_Prescl, &b_HLTriggerFirstPath_Prescl);
   fChain->SetBranchAddress("HLT_Activity_Ecal_SC7_v13", &HLT_Activity_Ecal_SC7_v13, &b_HLT_Activity_Ecal_SC7_v13);
   fChain->SetBranchAddress("HLT_Activity_Ecal_SC7_v13_Prescl", &HLT_Activity_Ecal_SC7_v13_Prescl, &b_HLT_Activity_Ecal_SC7_v13_Prescl);
   fChain->SetBranchAddress("HLT_BeamGas_HF_Beam1_v5", &HLT_BeamGas_HF_Beam1_v5, &b_HLT_BeamGas_HF_Beam1_v5);
   fChain->SetBranchAddress("HLT_BeamGas_HF_Beam1_v5_Prescl", &HLT_BeamGas_HF_Beam1_v5_Prescl, &b_HLT_BeamGas_HF_Beam1_v5_Prescl);
   fChain->SetBranchAddress("HLT_BeamGas_HF_Beam2_v5", &HLT_BeamGas_HF_Beam2_v5, &b_HLT_BeamGas_HF_Beam2_v5);
   fChain->SetBranchAddress("HLT_BeamGas_HF_Beam2_v5_Prescl", &HLT_BeamGas_HF_Beam2_v5_Prescl, &b_HLT_BeamGas_HF_Beam2_v5_Prescl);
   fChain->SetBranchAddress("HLT_BeamHalo_v13", &HLT_BeamHalo_v13, &b_HLT_BeamHalo_v13);
   fChain->SetBranchAddress("HLT_BeamHalo_v13_Prescl", &HLT_BeamHalo_v13_Prescl, &b_HLT_BeamHalo_v13_Prescl);
   fChain->SetBranchAddress("HLT_PAHcalUTCA_v1", &HLT_PAHcalUTCA_v1, &b_HLT_PAHcalUTCA_v1);
   fChain->SetBranchAddress("HLT_PAHcalUTCA_v1_Prescl", &HLT_PAHcalUTCA_v1_Prescl, &b_HLT_PAHcalUTCA_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAHcalPhiSym_v1", &HLT_PAHcalPhiSym_v1, &b_HLT_PAHcalPhiSym_v1);
   fChain->SetBranchAddress("HLT_PAHcalPhiSym_v1_Prescl", &HLT_PAHcalPhiSym_v1_Prescl, &b_HLT_PAHcalPhiSym_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAHcalNZS_v1", &HLT_PAHcalNZS_v1, &b_HLT_PAHcalNZS_v1);
   fChain->SetBranchAddress("HLT_PAHcalNZS_v1_Prescl", &HLT_PAHcalNZS_v1_Prescl, &b_HLT_PAHcalNZS_v1_Prescl);
   fChain->SetBranchAddress("HLT_GlobalRunHPDNoise_v8", &HLT_GlobalRunHPDNoise_v8, &b_HLT_GlobalRunHPDNoise_v8);
   fChain->SetBranchAddress("HLT_GlobalRunHPDNoise_v8_Prescl", &HLT_GlobalRunHPDNoise_v8_Prescl, &b_HLT_GlobalRunHPDNoise_v8_Prescl);
   fChain->SetBranchAddress("HLT_Physics_v5", &HLT_Physics_v5, &b_HLT_Physics_v5);
   fChain->SetBranchAddress("HLT_Physics_v5_Prescl", &HLT_Physics_v5_Prescl, &b_HLT_Physics_v5_Prescl);
   fChain->SetBranchAddress("DST_Physics_v5", &DST_Physics_v5, &b_DST_Physics_v5);
   fChain->SetBranchAddress("DST_Physics_v5_Prescl", &DST_Physics_v5_Prescl, &b_DST_Physics_v5_Prescl);
   fChain->SetBranchAddress("HLT_DTCalibration_v2", &HLT_DTCalibration_v2, &b_HLT_DTCalibration_v2);
   fChain->SetBranchAddress("HLT_DTCalibration_v2_Prescl", &HLT_DTCalibration_v2_Prescl, &b_HLT_DTCalibration_v2_Prescl);
   fChain->SetBranchAddress("HLT_EcalCalibration_v3", &HLT_EcalCalibration_v3, &b_HLT_EcalCalibration_v3);
   fChain->SetBranchAddress("HLT_EcalCalibration_v3_Prescl", &HLT_EcalCalibration_v3_Prescl, &b_HLT_EcalCalibration_v3_Prescl);
   fChain->SetBranchAddress("HLT_HcalCalibration_v3", &HLT_HcalCalibration_v3, &b_HLT_HcalCalibration_v3);
   fChain->SetBranchAddress("HLT_HcalCalibration_v3_Prescl", &HLT_HcalCalibration_v3_Prescl, &b_HLT_HcalCalibration_v3_Prescl);
   fChain->SetBranchAddress("HLT_TrackerCalibration_v3", &HLT_TrackerCalibration_v3, &b_HLT_TrackerCalibration_v3);
   fChain->SetBranchAddress("HLT_TrackerCalibration_v3_Prescl", &HLT_TrackerCalibration_v3_Prescl, &b_HLT_TrackerCalibration_v3_Prescl);
   fChain->SetBranchAddress("HLT_L1SingleMuOpen_AntiBPTX_v7", &HLT_L1SingleMuOpen_AntiBPTX_v7, &b_HLT_L1SingleMuOpen_AntiBPTX_v7);
   fChain->SetBranchAddress("HLT_L1SingleMuOpen_AntiBPTX_v7_Prescl", &HLT_L1SingleMuOpen_AntiBPTX_v7_Prescl, &b_HLT_L1SingleMuOpen_AntiBPTX_v7_Prescl);
   fChain->SetBranchAddress("HLT_L1TrackerCosmics_v7", &HLT_L1TrackerCosmics_v7, &b_HLT_L1TrackerCosmics_v7);
   fChain->SetBranchAddress("HLT_L1TrackerCosmics_v7_Prescl", &HLT_L1TrackerCosmics_v7_Prescl, &b_HLT_L1TrackerCosmics_v7_Prescl);
   fChain->SetBranchAddress("AlCa_PAEcalPi0EBonly_v1", &AlCa_PAEcalPi0EBonly_v1, &b_AlCa_PAEcalPi0EBonly_v1);
   fChain->SetBranchAddress("AlCa_PAEcalPi0EBonly_v1_Prescl", &AlCa_PAEcalPi0EBonly_v1_Prescl, &b_AlCa_PAEcalPi0EBonly_v1_Prescl);
   fChain->SetBranchAddress("AlCa_PAEcalPi0EEonly_v1", &AlCa_PAEcalPi0EEonly_v1, &b_AlCa_PAEcalPi0EEonly_v1);
   fChain->SetBranchAddress("AlCa_PAEcalPi0EEonly_v1_Prescl", &AlCa_PAEcalPi0EEonly_v1_Prescl, &b_AlCa_PAEcalPi0EEonly_v1_Prescl);
   fChain->SetBranchAddress("AlCa_PAEcalEtaEBonly_v1", &AlCa_PAEcalEtaEBonly_v1, &b_AlCa_PAEcalEtaEBonly_v1);
   fChain->SetBranchAddress("AlCa_PAEcalEtaEBonly_v1_Prescl", &AlCa_PAEcalEtaEBonly_v1_Prescl, &b_AlCa_PAEcalEtaEBonly_v1_Prescl);
   fChain->SetBranchAddress("AlCa_PAEcalEtaEEonly_v1", &AlCa_PAEcalEtaEEonly_v1, &b_AlCa_PAEcalEtaEEonly_v1);
   fChain->SetBranchAddress("AlCa_PAEcalEtaEEonly_v1_Prescl", &AlCa_PAEcalEtaEEonly_v1_Prescl, &b_AlCa_PAEcalEtaEEonly_v1_Prescl);
   fChain->SetBranchAddress("AlCa_EcalPhiSym_v13", &AlCa_EcalPhiSym_v13, &b_AlCa_EcalPhiSym_v13);
   fChain->SetBranchAddress("AlCa_EcalPhiSym_v13_Prescl", &AlCa_EcalPhiSym_v13_Prescl, &b_AlCa_EcalPhiSym_v13_Prescl);
   fChain->SetBranchAddress("AlCa_RPCMuonNoTriggers_v9", &AlCa_RPCMuonNoTriggers_v9, &b_AlCa_RPCMuonNoTriggers_v9);
   fChain->SetBranchAddress("AlCa_RPCMuonNoTriggers_v9_Prescl", &AlCa_RPCMuonNoTriggers_v9_Prescl, &b_AlCa_RPCMuonNoTriggers_v9_Prescl);
   fChain->SetBranchAddress("AlCa_RPCMuonNoHits_v9", &AlCa_RPCMuonNoHits_v9, &b_AlCa_RPCMuonNoHits_v9);
   fChain->SetBranchAddress("AlCa_RPCMuonNoHits_v9_Prescl", &AlCa_RPCMuonNoHits_v9_Prescl, &b_AlCa_RPCMuonNoHits_v9_Prescl);
   fChain->SetBranchAddress("AlCa_RPCMuonNormalisation_v9", &AlCa_RPCMuonNormalisation_v9, &b_AlCa_RPCMuonNormalisation_v9);
   fChain->SetBranchAddress("AlCa_RPCMuonNormalisation_v9_Prescl", &AlCa_RPCMuonNormalisation_v9_Prescl, &b_AlCa_RPCMuonNormalisation_v9_Prescl);
   fChain->SetBranchAddress("AlCa_LumiPixels_v8", &AlCa_LumiPixels_v8, &b_AlCa_LumiPixels_v8);
   fChain->SetBranchAddress("AlCa_LumiPixels_v8_Prescl", &AlCa_LumiPixels_v8_Prescl, &b_AlCa_LumiPixels_v8_Prescl);
   fChain->SetBranchAddress("AlCa_LumiPixels_ZeroBias_v4", &AlCa_LumiPixels_ZeroBias_v4, &b_AlCa_LumiPixels_ZeroBias_v4);
   fChain->SetBranchAddress("AlCa_LumiPixels_ZeroBias_v4_Prescl", &AlCa_LumiPixels_ZeroBias_v4_Prescl, &b_AlCa_LumiPixels_ZeroBias_v4_Prescl);
   fChain->SetBranchAddress("AlCa_LumiPixels_Random_v1", &AlCa_LumiPixels_Random_v1, &b_AlCa_LumiPixels_Random_v1);
   fChain->SetBranchAddress("AlCa_LumiPixels_Random_v1_Prescl", &AlCa_LumiPixels_Random_v1_Prescl, &b_AlCa_LumiPixels_Random_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAL1SingleJet16_v1", &HLT_PAL1SingleJet16_v1, &b_HLT_PAL1SingleJet16_v1);
   fChain->SetBranchAddress("HLT_PAL1SingleJet16_v1_Prescl", &HLT_PAL1SingleJet16_v1_Prescl, &b_HLT_PAL1SingleJet16_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAL1SingleJet36_v1", &HLT_PAL1SingleJet36_v1, &b_HLT_PAL1SingleJet36_v1);
   fChain->SetBranchAddress("HLT_PAL1SingleJet36_v1_Prescl", &HLT_PAL1SingleJet36_v1_Prescl, &b_HLT_PAL1SingleJet36_v1_Prescl);
   fChain->SetBranchAddress("HLT_PASingleForJet15_v1", &HLT_PASingleForJet15_v1, &b_HLT_PASingleForJet15_v1);
   fChain->SetBranchAddress("HLT_PASingleForJet15_v1_Prescl", &HLT_PASingleForJet15_v1_Prescl, &b_HLT_PASingleForJet15_v1_Prescl);
   fChain->SetBranchAddress("HLT_PASingleForJet25_v1", &HLT_PASingleForJet25_v1, &b_HLT_PASingleForJet25_v1);
   fChain->SetBranchAddress("HLT_PASingleForJet25_v1_Prescl", &HLT_PASingleForJet25_v1_Prescl, &b_HLT_PASingleForJet25_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAJet20_NoJetID_v1", &HLT_PAJet20_NoJetID_v1, &b_HLT_PAJet20_NoJetID_v1);
   fChain->SetBranchAddress("HLT_PAJet20_NoJetID_v1_Prescl", &HLT_PAJet20_NoJetID_v1_Prescl, &b_HLT_PAJet20_NoJetID_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAJet40_NoJetID_v1", &HLT_PAJet40_NoJetID_v1, &b_HLT_PAJet40_NoJetID_v1);
   fChain->SetBranchAddress("HLT_PAJet40_NoJetID_v1_Prescl", &HLT_PAJet40_NoJetID_v1_Prescl, &b_HLT_PAJet40_NoJetID_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAJet60_NoJetID_v1", &HLT_PAJet60_NoJetID_v1, &b_HLT_PAJet60_NoJetID_v1);
   fChain->SetBranchAddress("HLT_PAJet60_NoJetID_v1_Prescl", &HLT_PAJet60_NoJetID_v1_Prescl, &b_HLT_PAJet60_NoJetID_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAJet80_NoJetID_v1", &HLT_PAJet80_NoJetID_v1, &b_HLT_PAJet80_NoJetID_v1);
   fChain->SetBranchAddress("HLT_PAJet80_NoJetID_v1_Prescl", &HLT_PAJet80_NoJetID_v1_Prescl, &b_HLT_PAJet80_NoJetID_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAJet100_NoJetID_v1", &HLT_PAJet100_NoJetID_v1, &b_HLT_PAJet100_NoJetID_v1);
   fChain->SetBranchAddress("HLT_PAJet100_NoJetID_v1_Prescl", &HLT_PAJet100_NoJetID_v1_Prescl, &b_HLT_PAJet100_NoJetID_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAJet120_NoJetID_v1", &HLT_PAJet120_NoJetID_v1, &b_HLT_PAJet120_NoJetID_v1);
   fChain->SetBranchAddress("HLT_PAJet120_NoJetID_v1_Prescl", &HLT_PAJet120_NoJetID_v1_Prescl, &b_HLT_PAJet120_NoJetID_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAForJet20Eta2_v1", &HLT_PAForJet20Eta2_v1, &b_HLT_PAForJet20Eta2_v1);
   fChain->SetBranchAddress("HLT_PAForJet20Eta2_v1_Prescl", &HLT_PAForJet20Eta2_v1_Prescl, &b_HLT_PAForJet20Eta2_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAForJet40Eta2_v1", &HLT_PAForJet40Eta2_v1, &b_HLT_PAForJet40Eta2_v1);
   fChain->SetBranchAddress("HLT_PAForJet40Eta2_v1_Prescl", &HLT_PAForJet40Eta2_v1_Prescl, &b_HLT_PAForJet40Eta2_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAForJet60Eta2_v1", &HLT_PAForJet60Eta2_v1, &b_HLT_PAForJet60Eta2_v1);
   fChain->SetBranchAddress("HLT_PAForJet60Eta2_v1_Prescl", &HLT_PAForJet60Eta2_v1_Prescl, &b_HLT_PAForJet60Eta2_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAForJet80Eta2_v1", &HLT_PAForJet80Eta2_v1, &b_HLT_PAForJet80Eta2_v1);
   fChain->SetBranchAddress("HLT_PAForJet80Eta2_v1_Prescl", &HLT_PAForJet80Eta2_v1_Prescl, &b_HLT_PAForJet80Eta2_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAForJet100Eta2_v1", &HLT_PAForJet100Eta2_v1, &b_HLT_PAForJet100Eta2_v1);
   fChain->SetBranchAddress("HLT_PAForJet100Eta2_v1_Prescl", &HLT_PAForJet100Eta2_v1_Prescl, &b_HLT_PAForJet100Eta2_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAForJet20Eta3_v1", &HLT_PAForJet20Eta3_v1, &b_HLT_PAForJet20Eta3_v1);
   fChain->SetBranchAddress("HLT_PAForJet20Eta3_v1_Prescl", &HLT_PAForJet20Eta3_v1_Prescl, &b_HLT_PAForJet20Eta3_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAForJet40Eta3_v1", &HLT_PAForJet40Eta3_v1, &b_HLT_PAForJet40Eta3_v1);
   fChain->SetBranchAddress("HLT_PAForJet40Eta3_v1_Prescl", &HLT_PAForJet40Eta3_v1_Prescl, &b_HLT_PAForJet40Eta3_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAForJet60Eta3_v1", &HLT_PAForJet60Eta3_v1, &b_HLT_PAForJet60Eta3_v1);
   fChain->SetBranchAddress("HLT_PAForJet60Eta3_v1_Prescl", &HLT_PAForJet60Eta3_v1_Prescl, &b_HLT_PAForJet60Eta3_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAForJet80Eta3_v1", &HLT_PAForJet80Eta3_v1, &b_HLT_PAForJet80Eta3_v1);
   fChain->SetBranchAddress("HLT_PAForJet80Eta3_v1_Prescl", &HLT_PAForJet80Eta3_v1_Prescl, &b_HLT_PAForJet80Eta3_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAForJet100Eta3_v1", &HLT_PAForJet100Eta3_v1, &b_HLT_PAForJet100Eta3_v1);
   fChain->SetBranchAddress("HLT_PAForJet100Eta3_v1_Prescl", &HLT_PAForJet100Eta3_v1_Prescl, &b_HLT_PAForJet100Eta3_v1_Prescl);
   fChain->SetBranchAddress("HLT_PATripleJet20_20_20_v1", &HLT_PATripleJet20_20_20_v1, &b_HLT_PATripleJet20_20_20_v1);
   fChain->SetBranchAddress("HLT_PATripleJet20_20_20_v1_Prescl", &HLT_PATripleJet20_20_20_v1_Prescl, &b_HLT_PATripleJet20_20_20_v1_Prescl);
   fChain->SetBranchAddress("HLT_PATripleJet40_20_20_v1", &HLT_PATripleJet40_20_20_v1, &b_HLT_PATripleJet40_20_20_v1);
   fChain->SetBranchAddress("HLT_PATripleJet40_20_20_v1_Prescl", &HLT_PATripleJet40_20_20_v1_Prescl, &b_HLT_PATripleJet40_20_20_v1_Prescl);
   fChain->SetBranchAddress("HLT_PATripleJet60_20_20_v1", &HLT_PATripleJet60_20_20_v1, &b_HLT_PATripleJet60_20_20_v1);
   fChain->SetBranchAddress("HLT_PATripleJet60_20_20_v1_Prescl", &HLT_PATripleJet60_20_20_v1_Prescl, &b_HLT_PATripleJet60_20_20_v1_Prescl);
   fChain->SetBranchAddress("HLT_PATripleJet80_20_20_v1", &HLT_PATripleJet80_20_20_v1, &b_HLT_PATripleJet80_20_20_v1);
   fChain->SetBranchAddress("HLT_PATripleJet80_20_20_v1_Prescl", &HLT_PATripleJet80_20_20_v1_Prescl, &b_HLT_PATripleJet80_20_20_v1_Prescl);
   fChain->SetBranchAddress("HLT_PATripleJet100_20_20_v1", &HLT_PATripleJet100_20_20_v1, &b_HLT_PATripleJet100_20_20_v1);
   fChain->SetBranchAddress("HLT_PATripleJet100_20_20_v1_Prescl", &HLT_PATripleJet100_20_20_v1_Prescl, &b_HLT_PATripleJet100_20_20_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAJet40ETM30_v1", &HLT_PAJet40ETM30_v1, &b_HLT_PAJet40ETM30_v1);
   fChain->SetBranchAddress("HLT_PAJet40ETM30_v1_Prescl", &HLT_PAJet40ETM30_v1_Prescl, &b_HLT_PAJet40ETM30_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAJet60ETM30_v1", &HLT_PAJet60ETM30_v1, &b_HLT_PAJet60ETM30_v1);
   fChain->SetBranchAddress("HLT_PAJet60ETM30_v1_Prescl", &HLT_PAJet60ETM30_v1_Prescl, &b_HLT_PAJet60ETM30_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAL1DoubleMu0_v1", &HLT_PAL1DoubleMu0_v1, &b_HLT_PAL1DoubleMu0_v1);
   fChain->SetBranchAddress("HLT_PAL1DoubleMu0_v1_Prescl", &HLT_PAL1DoubleMu0_v1_Prescl, &b_HLT_PAL1DoubleMu0_v1_Prescl);
   fChain->SetBranchAddress("HLT_PADimuon0_NoVertexing_v1", &HLT_PADimuon0_NoVertexing_v1, &b_HLT_PADimuon0_NoVertexing_v1);
   fChain->SetBranchAddress("HLT_PADimuon0_NoVertexing_v1_Prescl", &HLT_PADimuon0_NoVertexing_v1_Prescl, &b_HLT_PADimuon0_NoVertexing_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAL1DoubleMu0_HighQ_v1", &HLT_PAL1DoubleMu0_HighQ_v1, &b_HLT_PAL1DoubleMu0_HighQ_v1);
   fChain->SetBranchAddress("HLT_PAL1DoubleMu0_HighQ_v1_Prescl", &HLT_PAL1DoubleMu0_HighQ_v1_Prescl, &b_HLT_PAL1DoubleMu0_HighQ_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAL1DoubleMuOpen_v1", &HLT_PAL1DoubleMuOpen_v1, &b_HLT_PAL1DoubleMuOpen_v1);
   fChain->SetBranchAddress("HLT_PAL1DoubleMuOpen_v1_Prescl", &HLT_PAL1DoubleMuOpen_v1_Prescl, &b_HLT_PAL1DoubleMuOpen_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAL2DoubleMu3_v1", &HLT_PAL2DoubleMu3_v1, &b_HLT_PAL2DoubleMu3_v1);
   fChain->SetBranchAddress("HLT_PAL2DoubleMu3_v1_Prescl", &HLT_PAL2DoubleMu3_v1_Prescl, &b_HLT_PAL2DoubleMu3_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAMu3_v1", &HLT_PAMu3_v1, &b_HLT_PAMu3_v1);
   fChain->SetBranchAddress("HLT_PAMu3_v1_Prescl", &HLT_PAMu3_v1_Prescl, &b_HLT_PAMu3_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAMu7_v1", &HLT_PAMu7_v1, &b_HLT_PAMu7_v1);
   fChain->SetBranchAddress("HLT_PAMu7_v1_Prescl", &HLT_PAMu7_v1_Prescl, &b_HLT_PAMu7_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAMu12_v1", &HLT_PAMu12_v1, &b_HLT_PAMu12_v1);
   fChain->SetBranchAddress("HLT_PAMu12_v1_Prescl", &HLT_PAMu12_v1_Prescl, &b_HLT_PAMu12_v1_Prescl);
   fChain->SetBranchAddress("HLT_Mu15_eta2p1_v5", &HLT_Mu15_eta2p1_v5, &b_HLT_Mu15_eta2p1_v5);
   fChain->SetBranchAddress("HLT_Mu15_eta2p1_v5_Prescl", &HLT_Mu15_eta2p1_v5_Prescl, &b_HLT_Mu15_eta2p1_v5_Prescl);
   fChain->SetBranchAddress("HLT_PABTagMu_Jet20_Mu4_v1", &HLT_PABTagMu_Jet20_Mu4_v1, &b_HLT_PABTagMu_Jet20_Mu4_v1);
   fChain->SetBranchAddress("HLT_PABTagMu_Jet20_Mu4_v1_Prescl", &HLT_PABTagMu_Jet20_Mu4_v1_Prescl, &b_HLT_PABTagMu_Jet20_Mu4_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAMu3PFJet20_v1", &HLT_PAMu3PFJet20_v1, &b_HLT_PAMu3PFJet20_v1);
   fChain->SetBranchAddress("HLT_PAMu3PFJet20_v1_Prescl", &HLT_PAMu3PFJet20_v1_Prescl, &b_HLT_PAMu3PFJet20_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAMu3PFJet40_v1", &HLT_PAMu3PFJet40_v1, &b_HLT_PAMu3PFJet40_v1);
   fChain->SetBranchAddress("HLT_PAMu3PFJet40_v1_Prescl", &HLT_PAMu3PFJet40_v1_Prescl, &b_HLT_PAMu3PFJet40_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAMu7PFJet20_v1", &HLT_PAMu7PFJet20_v1, &b_HLT_PAMu7PFJet20_v1);
   fChain->SetBranchAddress("HLT_PAMu7PFJet20_v1_Prescl", &HLT_PAMu7PFJet20_v1_Prescl, &b_HLT_PAMu7PFJet20_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAPhoton10_NoCaloIdVL_v1", &HLT_PAPhoton10_NoCaloIdVL_v1, &b_HLT_PAPhoton10_NoCaloIdVL_v1);
   fChain->SetBranchAddress("HLT_PAPhoton10_NoCaloIdVL_v1_Prescl", &HLT_PAPhoton10_NoCaloIdVL_v1_Prescl, &b_HLT_PAPhoton10_NoCaloIdVL_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAPhoton15_NoCaloIdVL_v1", &HLT_PAPhoton15_NoCaloIdVL_v1, &b_HLT_PAPhoton15_NoCaloIdVL_v1);
   fChain->SetBranchAddress("HLT_PAPhoton15_NoCaloIdVL_v1_Prescl", &HLT_PAPhoton15_NoCaloIdVL_v1_Prescl, &b_HLT_PAPhoton15_NoCaloIdVL_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAPhoton20_NoCaloIdVL_v1", &HLT_PAPhoton20_NoCaloIdVL_v1, &b_HLT_PAPhoton20_NoCaloIdVL_v1);
   fChain->SetBranchAddress("HLT_PAPhoton20_NoCaloIdVL_v1_Prescl", &HLT_PAPhoton20_NoCaloIdVL_v1_Prescl, &b_HLT_PAPhoton20_NoCaloIdVL_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAPhoton30_NoCaloIdVL_v1", &HLT_PAPhoton30_NoCaloIdVL_v1, &b_HLT_PAPhoton30_NoCaloIdVL_v1);
   fChain->SetBranchAddress("HLT_PAPhoton30_NoCaloIdVL_v1_Prescl", &HLT_PAPhoton30_NoCaloIdVL_v1_Prescl, &b_HLT_PAPhoton30_NoCaloIdVL_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAPhoton40_NoCaloIdVL_v1", &HLT_PAPhoton40_NoCaloIdVL_v1, &b_HLT_PAPhoton40_NoCaloIdVL_v1);
   fChain->SetBranchAddress("HLT_PAPhoton40_NoCaloIdVL_v1_Prescl", &HLT_PAPhoton40_NoCaloIdVL_v1_Prescl, &b_HLT_PAPhoton40_NoCaloIdVL_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAPhoton60_NoCaloIdVL_v1", &HLT_PAPhoton60_NoCaloIdVL_v1, &b_HLT_PAPhoton60_NoCaloIdVL_v1);
   fChain->SetBranchAddress("HLT_PAPhoton60_NoCaloIdVL_v1_Prescl", &HLT_PAPhoton60_NoCaloIdVL_v1_Prescl, &b_HLT_PAPhoton60_NoCaloIdVL_v1_Prescl);
   fChain->SetBranchAddress("HLT_Ele22_CaloIdL_CaloIsoVL_v6", &HLT_Ele22_CaloIdL_CaloIsoVL_v6, &b_HLT_Ele22_CaloIdL_CaloIsoVL_v6);
   fChain->SetBranchAddress("HLT_Ele22_CaloIdL_CaloIsoVL_v6_Prescl", &HLT_Ele22_CaloIdL_CaloIsoVL_v6_Prescl, &b_HLT_Ele22_CaloIdL_CaloIsoVL_v6_Prescl);
   fChain->SetBranchAddress("HLT_PAPhoton10_TightCaloIdVL_v1", &HLT_PAPhoton10_TightCaloIdVL_v1, &b_HLT_PAPhoton10_TightCaloIdVL_v1);
   fChain->SetBranchAddress("HLT_PAPhoton10_TightCaloIdVL_v1_Prescl", &HLT_PAPhoton10_TightCaloIdVL_v1_Prescl, &b_HLT_PAPhoton10_TightCaloIdVL_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAPhoton15_TightCaloIdVL_v1", &HLT_PAPhoton15_TightCaloIdVL_v1, &b_HLT_PAPhoton15_TightCaloIdVL_v1);
   fChain->SetBranchAddress("HLT_PAPhoton15_TightCaloIdVL_v1_Prescl", &HLT_PAPhoton15_TightCaloIdVL_v1_Prescl, &b_HLT_PAPhoton15_TightCaloIdVL_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAPhoton20_TightCaloIdVL_v1", &HLT_PAPhoton20_TightCaloIdVL_v1, &b_HLT_PAPhoton20_TightCaloIdVL_v1);
   fChain->SetBranchAddress("HLT_PAPhoton20_TightCaloIdVL_v1_Prescl", &HLT_PAPhoton20_TightCaloIdVL_v1_Prescl, &b_HLT_PAPhoton20_TightCaloIdVL_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAPhoton30_TightCaloIdVL_v1", &HLT_PAPhoton30_TightCaloIdVL_v1, &b_HLT_PAPhoton30_TightCaloIdVL_v1);
   fChain->SetBranchAddress("HLT_PAPhoton30_TightCaloIdVL_v1_Prescl", &HLT_PAPhoton30_TightCaloIdVL_v1_Prescl, &b_HLT_PAPhoton30_TightCaloIdVL_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAPhoton40_TightCaloIdVL_v1", &HLT_PAPhoton40_TightCaloIdVL_v1, &b_HLT_PAPhoton40_TightCaloIdVL_v1);
   fChain->SetBranchAddress("HLT_PAPhoton40_TightCaloIdVL_v1_Prescl", &HLT_PAPhoton40_TightCaloIdVL_v1_Prescl, &b_HLT_PAPhoton40_TightCaloIdVL_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAPhoton10_Photon10_NoCaloIdVL_v1", &HLT_PAPhoton10_Photon10_NoCaloIdVL_v1, &b_HLT_PAPhoton10_Photon10_NoCaloIdVL_v1);
   fChain->SetBranchAddress("HLT_PAPhoton10_Photon10_NoCaloIdVL_v1_Prescl", &HLT_PAPhoton10_Photon10_NoCaloIdVL_v1_Prescl, &b_HLT_PAPhoton10_Photon10_NoCaloIdVL_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAPhoton15_Photon10_NoCaloIdVL_v1", &HLT_PAPhoton15_Photon10_NoCaloIdVL_v1, &b_HLT_PAPhoton15_Photon10_NoCaloIdVL_v1);
   fChain->SetBranchAddress("HLT_PAPhoton15_Photon10_NoCaloIdVL_v1_Prescl", &HLT_PAPhoton15_Photon10_NoCaloIdVL_v1_Prescl, &b_HLT_PAPhoton15_Photon10_NoCaloIdVL_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAPhoton20_Photon15_NoCaloIdVL_v1", &HLT_PAPhoton20_Photon15_NoCaloIdVL_v1, &b_HLT_PAPhoton20_Photon15_NoCaloIdVL_v1);
   fChain->SetBranchAddress("HLT_PAPhoton20_Photon15_NoCaloIdVL_v1_Prescl", &HLT_PAPhoton20_Photon15_NoCaloIdVL_v1_Prescl, &b_HLT_PAPhoton20_Photon15_NoCaloIdVL_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAPhoton20_Photon20_NoCaloIdVL_v1", &HLT_PAPhoton20_Photon20_NoCaloIdVL_v1, &b_HLT_PAPhoton20_Photon20_NoCaloIdVL_v1);
   fChain->SetBranchAddress("HLT_PAPhoton20_Photon20_NoCaloIdVL_v1_Prescl", &HLT_PAPhoton20_Photon20_NoCaloIdVL_v1_Prescl, &b_HLT_PAPhoton20_Photon20_NoCaloIdVL_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAPhoton30_Photon30_NoCaloIdVL_v1", &HLT_PAPhoton30_Photon30_NoCaloIdVL_v1, &b_HLT_PAPhoton30_Photon30_NoCaloIdVL_v1);
   fChain->SetBranchAddress("HLT_PAPhoton30_Photon30_NoCaloIdVL_v1_Prescl", &HLT_PAPhoton30_Photon30_NoCaloIdVL_v1_Prescl, &b_HLT_PAPhoton30_Photon30_NoCaloIdVL_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAPhoton10_Photon10_TightCaloIdVL_v1", &HLT_PAPhoton10_Photon10_TightCaloIdVL_v1, &b_HLT_PAPhoton10_Photon10_TightCaloIdVL_v1);
   fChain->SetBranchAddress("HLT_PAPhoton10_Photon10_TightCaloIdVL_v1_Prescl", &HLT_PAPhoton10_Photon10_TightCaloIdVL_v1_Prescl, &b_HLT_PAPhoton10_Photon10_TightCaloIdVL_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAPhoton15_Photon10_TightCaloIdVL_v1", &HLT_PAPhoton15_Photon10_TightCaloIdVL_v1, &b_HLT_PAPhoton15_Photon10_TightCaloIdVL_v1);
   fChain->SetBranchAddress("HLT_PAPhoton15_Photon10_TightCaloIdVL_v1_Prescl", &HLT_PAPhoton15_Photon10_TightCaloIdVL_v1_Prescl, &b_HLT_PAPhoton15_Photon10_TightCaloIdVL_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAPhoton20_Photon15_TightCaloIdVL_v1", &HLT_PAPhoton20_Photon15_TightCaloIdVL_v1, &b_HLT_PAPhoton20_Photon15_TightCaloIdVL_v1);
   fChain->SetBranchAddress("HLT_PAPhoton20_Photon15_TightCaloIdVL_v1_Prescl", &HLT_PAPhoton20_Photon15_TightCaloIdVL_v1_Prescl, &b_HLT_PAPhoton20_Photon15_TightCaloIdVL_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAPixelTracks_Multiplicity100_v2", &HLT_PAPixelTracks_Multiplicity100_v2, &b_HLT_PAPixelTracks_Multiplicity100_v2);
   fChain->SetBranchAddress("HLT_PAPixelTracks_Multiplicity100_v2_Prescl", &HLT_PAPixelTracks_Multiplicity100_v2_Prescl, &b_HLT_PAPixelTracks_Multiplicity100_v2_Prescl);
   fChain->SetBranchAddress("HLT_PAPixelTracks_Multiplicity130_v2", &HLT_PAPixelTracks_Multiplicity130_v2, &b_HLT_PAPixelTracks_Multiplicity130_v2);
   fChain->SetBranchAddress("HLT_PAPixelTracks_Multiplicity130_v2_Prescl", &HLT_PAPixelTracks_Multiplicity130_v2_Prescl, &b_HLT_PAPixelTracks_Multiplicity130_v2_Prescl);
   fChain->SetBranchAddress("HLT_PAPixelTracks_Multiplicity160_v2", &HLT_PAPixelTracks_Multiplicity160_v2, &b_HLT_PAPixelTracks_Multiplicity160_v2);
   fChain->SetBranchAddress("HLT_PAPixelTracks_Multiplicity160_v2_Prescl", &HLT_PAPixelTracks_Multiplicity160_v2_Prescl, &b_HLT_PAPixelTracks_Multiplicity160_v2_Prescl);
   fChain->SetBranchAddress("HLT_PAPixelTracks_Multiplicity190_v2", &HLT_PAPixelTracks_Multiplicity190_v2, &b_HLT_PAPixelTracks_Multiplicity190_v2);
   fChain->SetBranchAddress("HLT_PAPixelTracks_Multiplicity190_v2_Prescl", &HLT_PAPixelTracks_Multiplicity190_v2_Prescl, &b_HLT_PAPixelTracks_Multiplicity190_v2_Prescl);
   fChain->SetBranchAddress("HLT_PAPixelTracks_Multiplicity220_v2", &HLT_PAPixelTracks_Multiplicity220_v2, &b_HLT_PAPixelTracks_Multiplicity220_v2);
   fChain->SetBranchAddress("HLT_PAPixelTracks_Multiplicity220_v2_Prescl", &HLT_PAPixelTracks_Multiplicity220_v2_Prescl, &b_HLT_PAPixelTracks_Multiplicity220_v2_Prescl);
   fChain->SetBranchAddress("HLT_PAPixelTrackMultiplicity100_FullTrack12_v2", &HLT_PAPixelTrackMultiplicity100_FullTrack12_v2, &b_HLT_PAPixelTrackMultiplicity100_FullTrack12_v2);
   fChain->SetBranchAddress("HLT_PAPixelTrackMultiplicity100_FullTrack12_v2_Prescl", &HLT_PAPixelTrackMultiplicity100_FullTrack12_v2_Prescl, &b_HLT_PAPixelTrackMultiplicity100_FullTrack12_v2_Prescl);
   fChain->SetBranchAddress("HLT_PAPixelTrackMultiplicity130_FullTrack12_v2", &HLT_PAPixelTrackMultiplicity130_FullTrack12_v2, &b_HLT_PAPixelTrackMultiplicity130_FullTrack12_v2);
   fChain->SetBranchAddress("HLT_PAPixelTrackMultiplicity130_FullTrack12_v2_Prescl", &HLT_PAPixelTrackMultiplicity130_FullTrack12_v2_Prescl, &b_HLT_PAPixelTrackMultiplicity130_FullTrack12_v2_Prescl);
   fChain->SetBranchAddress("HLT_PAPixelTrackMultiplicity160_FullTrack12_v2", &HLT_PAPixelTrackMultiplicity160_FullTrack12_v2, &b_HLT_PAPixelTrackMultiplicity160_FullTrack12_v2);
   fChain->SetBranchAddress("HLT_PAPixelTrackMultiplicity160_FullTrack12_v2_Prescl", &HLT_PAPixelTrackMultiplicity160_FullTrack12_v2_Prescl, &b_HLT_PAPixelTrackMultiplicity160_FullTrack12_v2_Prescl);
   fChain->SetBranchAddress("HLT_PAFullTrack12_v2", &HLT_PAFullTrack12_v2, &b_HLT_PAFullTrack12_v2);
   fChain->SetBranchAddress("HLT_PAFullTrack12_v2_Prescl", &HLT_PAFullTrack12_v2_Prescl, &b_HLT_PAFullTrack12_v2_Prescl);
   fChain->SetBranchAddress("HLT_PAFullTrack20_v2", &HLT_PAFullTrack20_v2, &b_HLT_PAFullTrack20_v2);
   fChain->SetBranchAddress("HLT_PAFullTrack20_v2_Prescl", &HLT_PAFullTrack20_v2_Prescl, &b_HLT_PAFullTrack20_v2_Prescl);
   fChain->SetBranchAddress("HLT_PAFullTrack30_v2", &HLT_PAFullTrack30_v2, &b_HLT_PAFullTrack30_v2);
   fChain->SetBranchAddress("HLT_PAFullTrack30_v2_Prescl", &HLT_PAFullTrack30_v2_Prescl, &b_HLT_PAFullTrack30_v2_Prescl);
   fChain->SetBranchAddress("HLT_PAFullTrack50_v2", &HLT_PAFullTrack50_v2, &b_HLT_PAFullTrack50_v2);
   fChain->SetBranchAddress("HLT_PAFullTrack50_v2_Prescl", &HLT_PAFullTrack50_v2_Prescl, &b_HLT_PAFullTrack50_v2_Prescl);
   fChain->SetBranchAddress("HLT_PAPixelTrackMultiplicity140_Jet80_NoJetID_v2", &HLT_PAPixelTrackMultiplicity140_Jet80_NoJetID_v2, &b_HLT_PAPixelTrackMultiplicity140_Jet80_NoJetID_v2);
   fChain->SetBranchAddress("HLT_PAPixelTrackMultiplicity140_Jet80_NoJetID_v2_Prescl", &HLT_PAPixelTrackMultiplicity140_Jet80_NoJetID_v2_Prescl, &b_HLT_PAPixelTrackMultiplicity140_Jet80_NoJetID_v2_Prescl);
   fChain->SetBranchAddress("HLT_PAPixelTrackMultiplicity100_L2DoubleMu3_v1", &HLT_PAPixelTrackMultiplicity100_L2DoubleMu3_v1, &b_HLT_PAPixelTrackMultiplicity100_L2DoubleMu3_v1);
   fChain->SetBranchAddress("HLT_PAPixelTrackMultiplicity100_L2DoubleMu3_v1_Prescl", &HLT_PAPixelTrackMultiplicity100_L2DoubleMu3_v1_Prescl, &b_HLT_PAPixelTrackMultiplicity100_L2DoubleMu3_v1_Prescl);
   fChain->SetBranchAddress("HLT_PPPixelTracks_Multiplicity55_v1", &HLT_PPPixelTracks_Multiplicity55_v1, &b_HLT_PPPixelTracks_Multiplicity55_v1);
   fChain->SetBranchAddress("HLT_PPPixelTracks_Multiplicity55_v1_Prescl", &HLT_PPPixelTracks_Multiplicity55_v1_Prescl, &b_HLT_PPPixelTracks_Multiplicity55_v1_Prescl);
   fChain->SetBranchAddress("HLT_PPPixelTracks_Multiplicity70_v1", &HLT_PPPixelTracks_Multiplicity70_v1, &b_HLT_PPPixelTracks_Multiplicity70_v1);
   fChain->SetBranchAddress("HLT_PPPixelTracks_Multiplicity70_v1_Prescl", &HLT_PPPixelTracks_Multiplicity70_v1_Prescl, &b_HLT_PPPixelTracks_Multiplicity70_v1_Prescl);
   fChain->SetBranchAddress("HLT_PPPixelTracks_Multiplicity85_v1", &HLT_PPPixelTracks_Multiplicity85_v1, &b_HLT_PPPixelTracks_Multiplicity85_v1);
   fChain->SetBranchAddress("HLT_PPPixelTracks_Multiplicity85_v1_Prescl", &HLT_PPPixelTracks_Multiplicity85_v1_Prescl, &b_HLT_PPPixelTracks_Multiplicity85_v1_Prescl);
   fChain->SetBranchAddress("HLT_PPPixelTrackMultiplicity55_FullTrack12_v1", &HLT_PPPixelTrackMultiplicity55_FullTrack12_v1, &b_HLT_PPPixelTrackMultiplicity55_FullTrack12_v1);
   fChain->SetBranchAddress("HLT_PPPixelTrackMultiplicity55_FullTrack12_v1_Prescl", &HLT_PPPixelTrackMultiplicity55_FullTrack12_v1_Prescl, &b_HLT_PPPixelTrackMultiplicity55_FullTrack12_v1_Prescl);
   fChain->SetBranchAddress("HLT_PPPixelTrackMultiplicity70_FullTrack12_v1", &HLT_PPPixelTrackMultiplicity70_FullTrack12_v1, &b_HLT_PPPixelTrackMultiplicity70_FullTrack12_v1);
   fChain->SetBranchAddress("HLT_PPPixelTrackMultiplicity70_FullTrack12_v1_Prescl", &HLT_PPPixelTrackMultiplicity70_FullTrack12_v1_Prescl, &b_HLT_PPPixelTrackMultiplicity70_FullTrack12_v1_Prescl);
   fChain->SetBranchAddress("HLT_PPL1DoubleJetC36_v1", &HLT_PPL1DoubleJetC36_v1, &b_HLT_PPL1DoubleJetC36_v1);
   fChain->SetBranchAddress("HLT_PPL1DoubleJetC36_v1_Prescl", &HLT_PPL1DoubleJetC36_v1_Prescl, &b_HLT_PPL1DoubleJetC36_v1_Prescl);
   fChain->SetBranchAddress("HLT_PATech35_v1", &HLT_PATech35_v1, &b_HLT_PATech35_v1);
   fChain->SetBranchAddress("HLT_PATech35_v1_Prescl", &HLT_PATech35_v1_Prescl, &b_HLT_PATech35_v1_Prescl);
   fChain->SetBranchAddress("HLT_PATech35_HFSumET100_v3", &HLT_PATech35_HFSumET100_v3, &b_HLT_PATech35_HFSumET100_v3);
   fChain->SetBranchAddress("HLT_PATech35_HFSumET100_v3_Prescl", &HLT_PATech35_HFSumET100_v3_Prescl, &b_HLT_PATech35_HFSumET100_v3_Prescl);
   fChain->SetBranchAddress("HLT_PAHFSumET100_v3", &HLT_PAHFSumET100_v3, &b_HLT_PAHFSumET100_v3);
   fChain->SetBranchAddress("HLT_PAHFSumET100_v3_Prescl", &HLT_PAHFSumET100_v3_Prescl, &b_HLT_PAHFSumET100_v3_Prescl);
   fChain->SetBranchAddress("HLT_PAHFSumET140_v3", &HLT_PAHFSumET140_v3, &b_HLT_PAHFSumET140_v3);
   fChain->SetBranchAddress("HLT_PAHFSumET140_v3_Prescl", &HLT_PAHFSumET140_v3_Prescl, &b_HLT_PAHFSumET140_v3_Prescl);
   fChain->SetBranchAddress("HLT_PAHFSumET170_v3", &HLT_PAHFSumET170_v3, &b_HLT_PAHFSumET170_v3);
   fChain->SetBranchAddress("HLT_PAHFSumET170_v3_Prescl", &HLT_PAHFSumET170_v3_Prescl, &b_HLT_PAHFSumET170_v3_Prescl);
   fChain->SetBranchAddress("HLT_PAHFSumET210_v3", &HLT_PAHFSumET210_v3, &b_HLT_PAHFSumET210_v3);
   fChain->SetBranchAddress("HLT_PAHFSumET210_v3_Prescl", &HLT_PAHFSumET210_v3_Prescl, &b_HLT_PAHFSumET210_v3_Prescl);
   fChain->SetBranchAddress("HLT_PARomanPots_Tech52_v1", &HLT_PARomanPots_Tech52_v1, &b_HLT_PARomanPots_Tech52_v1);
   fChain->SetBranchAddress("HLT_PARomanPots_Tech52_v1_Prescl", &HLT_PARomanPots_Tech52_v1_Prescl, &b_HLT_PARomanPots_Tech52_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAL1Tech53_MB_v1", &HLT_PAL1Tech53_MB_v1, &b_HLT_PAL1Tech53_MB_v1);
   fChain->SetBranchAddress("HLT_PAL1Tech53_MB_v1_Prescl", &HLT_PAL1Tech53_MB_v1_Prescl, &b_HLT_PAL1Tech53_MB_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAL1Tech53_MB_SingleTrack_v1", &HLT_PAL1Tech53_MB_SingleTrack_v1, &b_HLT_PAL1Tech53_MB_SingleTrack_v1);
   fChain->SetBranchAddress("HLT_PAL1Tech53_MB_SingleTrack_v1_Prescl", &HLT_PAL1Tech53_MB_SingleTrack_v1_Prescl, &b_HLT_PAL1Tech53_MB_SingleTrack_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAL1Tech54_ZeroBias_v1", &HLT_PAL1Tech54_ZeroBias_v1, &b_HLT_PAL1Tech54_ZeroBias_v1);
   fChain->SetBranchAddress("HLT_PAL1Tech54_ZeroBias_v1_Prescl", &HLT_PAL1Tech54_ZeroBias_v1_Prescl, &b_HLT_PAL1Tech54_ZeroBias_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAT1minbias_Tech55_v1", &HLT_PAT1minbias_Tech55_v1, &b_HLT_PAT1minbias_Tech55_v1);
   fChain->SetBranchAddress("HLT_PAT1minbias_Tech55_v1_Prescl", &HLT_PAT1minbias_Tech55_v1_Prescl, &b_HLT_PAT1minbias_Tech55_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAL1Tech_HBHEHO_totalOR_v1", &HLT_PAL1Tech_HBHEHO_totalOR_v1, &b_HLT_PAL1Tech_HBHEHO_totalOR_v1);
   fChain->SetBranchAddress("HLT_PAL1Tech_HBHEHO_totalOR_v1_Prescl", &HLT_PAL1Tech_HBHEHO_totalOR_v1_Prescl, &b_HLT_PAL1Tech_HBHEHO_totalOR_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAL1Tech63_CASTORHaloMuon_v1", &HLT_PAL1Tech63_CASTORHaloMuon_v1, &b_HLT_PAL1Tech63_CASTORHaloMuon_v1);
   fChain->SetBranchAddress("HLT_PAL1Tech63_CASTORHaloMuon_v1_Prescl", &HLT_PAL1Tech63_CASTORHaloMuon_v1_Prescl, &b_HLT_PAL1Tech63_CASTORHaloMuon_v1_Prescl);
   fChain->SetBranchAddress("HLT_PACastorEmTotemLowMultiplicity_v1", &HLT_PACastorEmTotemLowMultiplicity_v1, &b_HLT_PACastorEmTotemLowMultiplicity_v1);
   fChain->SetBranchAddress("HLT_PACastorEmTotemLowMultiplicity_v1_Prescl", &HLT_PACastorEmTotemLowMultiplicity_v1_Prescl, &b_HLT_PACastorEmTotemLowMultiplicity_v1_Prescl);
   fChain->SetBranchAddress("HLT_PACastorEmNotHfCoincidencePm_v1", &HLT_PACastorEmNotHfCoincidencePm_v1, &b_HLT_PACastorEmNotHfCoincidencePm_v1);
   fChain->SetBranchAddress("HLT_PACastorEmNotHfCoincidencePm_v1_Prescl", &HLT_PACastorEmNotHfCoincidencePm_v1_Prescl, &b_HLT_PACastorEmNotHfCoincidencePm_v1_Prescl);
   fChain->SetBranchAddress("HLT_PACastorEmNotHfSingleChannel_v1", &HLT_PACastorEmNotHfSingleChannel_v1, &b_HLT_PACastorEmNotHfSingleChannel_v1);
   fChain->SetBranchAddress("HLT_PACastorEmNotHfSingleChannel_v1_Prescl", &HLT_PACastorEmNotHfSingleChannel_v1_Prescl, &b_HLT_PACastorEmNotHfSingleChannel_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAL1CastorTotalTotemLowMultiplicity_v1", &HLT_PAL1CastorTotalTotemLowMultiplicity_v1, &b_HLT_PAL1CastorTotalTotemLowMultiplicity_v1);
   fChain->SetBranchAddress("HLT_PAL1CastorTotalTotemLowMultiplicity_v1_Prescl", &HLT_PAL1CastorTotalTotemLowMultiplicity_v1_Prescl, &b_HLT_PAL1CastorTotalTotemLowMultiplicity_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAMinBiasHF_v1", &HLT_PAMinBiasHF_v1, &b_HLT_PAMinBiasHF_v1);
   fChain->SetBranchAddress("HLT_PAMinBiasHF_v1_Prescl", &HLT_PAMinBiasHF_v1_Prescl, &b_HLT_PAMinBiasHF_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAMinBiasHF_OR_v1", &HLT_PAMinBiasHF_OR_v1, &b_HLT_PAMinBiasHF_OR_v1);
   fChain->SetBranchAddress("HLT_PAMinBiasHF_OR_v1_Prescl", &HLT_PAMinBiasHF_OR_v1_Prescl, &b_HLT_PAMinBiasHF_OR_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAMinBiasBHC_v1", &HLT_PAMinBiasBHC_v1, &b_HLT_PAMinBiasBHC_v1);
   fChain->SetBranchAddress("HLT_PAMinBiasBHC_v1_Prescl", &HLT_PAMinBiasBHC_v1_Prescl, &b_HLT_PAMinBiasBHC_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAMinBiasBHC_OR_v1", &HLT_PAMinBiasBHC_OR_v1, &b_HLT_PAMinBiasBHC_OR_v1);
   fChain->SetBranchAddress("HLT_PAMinBiasBHC_OR_v1_Prescl", &HLT_PAMinBiasBHC_OR_v1_Prescl, &b_HLT_PAMinBiasBHC_OR_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAMinBiasHfOrBHC_v1", &HLT_PAMinBiasHfOrBHC_v1, &b_HLT_PAMinBiasHfOrBHC_v1);
   fChain->SetBranchAddress("HLT_PAMinBiasHfOrBHC_v1_Prescl", &HLT_PAMinBiasHfOrBHC_v1_Prescl, &b_HLT_PAMinBiasHfOrBHC_v1_Prescl);
   fChain->SetBranchAddress("HLT_PABptxPlusNotBptxMinus_v1", &HLT_PABptxPlusNotBptxMinus_v1, &b_HLT_PABptxPlusNotBptxMinus_v1);
   fChain->SetBranchAddress("HLT_PABptxPlusNotBptxMinus_v1_Prescl", &HLT_PABptxPlusNotBptxMinus_v1_Prescl, &b_HLT_PABptxPlusNotBptxMinus_v1_Prescl);
   fChain->SetBranchAddress("HLT_PABptxMinusNotBptxPlus_v1", &HLT_PABptxMinusNotBptxPlus_v1, &b_HLT_PABptxMinusNotBptxPlus_v1);
   fChain->SetBranchAddress("HLT_PABptxMinusNotBptxPlus_v1_Prescl", &HLT_PABptxMinusNotBptxPlus_v1_Prescl, &b_HLT_PABptxMinusNotBptxPlus_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAZeroBias_v1", &HLT_PAZeroBias_v1, &b_HLT_PAZeroBias_v1);
   fChain->SetBranchAddress("HLT_PAZeroBias_v1_Prescl", &HLT_PAZeroBias_v1_Prescl, &b_HLT_PAZeroBias_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1", &HLT_PAZeroBiasPixel_SingleTrack_v1, &b_HLT_PAZeroBiasPixel_SingleTrack_v1);
   fChain->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl", &HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl, &b_HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAHFOR_SingleTrack_v1", &HLT_PAHFOR_SingleTrack_v1, &b_HLT_PAHFOR_SingleTrack_v1);
   fChain->SetBranchAddress("HLT_PAHFOR_SingleTrack_v1_Prescl", &HLT_PAHFOR_SingleTrack_v1_Prescl, &b_HLT_PAHFOR_SingleTrack_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAZeroBiasPixel_DoubleTrack_v1", &HLT_PAZeroBiasPixel_DoubleTrack_v1, &b_HLT_PAZeroBiasPixel_DoubleTrack_v1);
   fChain->SetBranchAddress("HLT_PAZeroBiasPixel_DoubleTrack_v1_Prescl", &HLT_PAZeroBiasPixel_DoubleTrack_v1_Prescl, &b_HLT_PAZeroBiasPixel_DoubleTrack_v1_Prescl);
   fChain->SetBranchAddress("HLT_PADoubleMu4_Acoplanarity03_v1", &HLT_PADoubleMu4_Acoplanarity03_v1, &b_HLT_PADoubleMu4_Acoplanarity03_v1);
   fChain->SetBranchAddress("HLT_PADoubleMu4_Acoplanarity03_v1_Prescl", &HLT_PADoubleMu4_Acoplanarity03_v1_Prescl, &b_HLT_PADoubleMu4_Acoplanarity03_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAExclDijet35_HFOR_v1", &HLT_PAExclDijet35_HFOR_v1, &b_HLT_PAExclDijet35_HFOR_v1);
   fChain->SetBranchAddress("HLT_PAExclDijet35_HFOR_v1_Prescl", &HLT_PAExclDijet35_HFOR_v1_Prescl, &b_HLT_PAExclDijet35_HFOR_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAExclDijet35_HFAND_v1", &HLT_PAExclDijet35_HFAND_v1, &b_HLT_PAExclDijet35_HFAND_v1);
   fChain->SetBranchAddress("HLT_PAExclDijet35_HFAND_v1_Prescl", &HLT_PAExclDijet35_HFAND_v1_Prescl, &b_HLT_PAExclDijet35_HFAND_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAL1DoubleEG3_FwdVeto_v1", &HLT_PAL1DoubleEG3_FwdVeto_v1, &b_HLT_PAL1DoubleEG3_FwdVeto_v1);
   fChain->SetBranchAddress("HLT_PAL1DoubleEG3_FwdVeto_v1_Prescl", &HLT_PAL1DoubleEG3_FwdVeto_v1_Prescl, &b_HLT_PAL1DoubleEG3_FwdVeto_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAL1SingleJet52_TotemDiffractive_v1", &HLT_PAL1SingleJet52_TotemDiffractive_v1, &b_HLT_PAL1SingleJet52_TotemDiffractive_v1);
   fChain->SetBranchAddress("HLT_PAL1SingleJet52_TotemDiffractive_v1_Prescl", &HLT_PAL1SingleJet52_TotemDiffractive_v1_Prescl, &b_HLT_PAL1SingleJet52_TotemDiffractive_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAL1SingleMu20_TotemDiffractive_v1", &HLT_PAL1SingleMu20_TotemDiffractive_v1, &b_HLT_PAL1SingleMu20_TotemDiffractive_v1);
   fChain->SetBranchAddress("HLT_PAL1SingleMu20_TotemDiffractive_v1_Prescl", &HLT_PAL1SingleMu20_TotemDiffractive_v1_Prescl, &b_HLT_PAL1SingleMu20_TotemDiffractive_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAL1SingleEG20_TotemDiffractive_v1", &HLT_PAL1SingleEG20_TotemDiffractive_v1, &b_HLT_PAL1SingleEG20_TotemDiffractive_v1);
   fChain->SetBranchAddress("HLT_PAL1SingleEG20_TotemDiffractive_v1_Prescl", &HLT_PAL1SingleEG20_TotemDiffractive_v1_Prescl, &b_HLT_PAL1SingleEG20_TotemDiffractive_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAL1DoubleJet20_TotemDiffractive_v1", &HLT_PAL1DoubleJet20_TotemDiffractive_v1, &b_HLT_PAL1DoubleJet20_TotemDiffractive_v1);
   fChain->SetBranchAddress("HLT_PAL1DoubleJet20_TotemDiffractive_v1_Prescl", &HLT_PAL1DoubleJet20_TotemDiffractive_v1_Prescl, &b_HLT_PAL1DoubleJet20_TotemDiffractive_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAL1DoubleJetC36_TotemDiffractive_v1", &HLT_PAL1DoubleJetC36_TotemDiffractive_v1, &b_HLT_PAL1DoubleJetC36_TotemDiffractive_v1);
   fChain->SetBranchAddress("HLT_PAL1DoubleJetC36_TotemDiffractive_v1_Prescl", &HLT_PAL1DoubleJetC36_TotemDiffractive_v1_Prescl, &b_HLT_PAL1DoubleJetC36_TotemDiffractive_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAL1DoubleMu5_TotemDiffractive_v1", &HLT_PAL1DoubleMu5_TotemDiffractive_v1, &b_HLT_PAL1DoubleMu5_TotemDiffractive_v1);
   fChain->SetBranchAddress("HLT_PAL1DoubleMu5_TotemDiffractive_v1_Prescl", &HLT_PAL1DoubleMu5_TotemDiffractive_v1_Prescl, &b_HLT_PAL1DoubleMu5_TotemDiffractive_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAL1DoubleEG5_TotemDiffractive_v1", &HLT_PAL1DoubleEG5_TotemDiffractive_v1, &b_HLT_PAL1DoubleEG5_TotemDiffractive_v1);
   fChain->SetBranchAddress("HLT_PAL1DoubleEG5_TotemDiffractive_v1_Prescl", &HLT_PAL1DoubleEG5_TotemDiffractive_v1_Prescl, &b_HLT_PAL1DoubleEG5_TotemDiffractive_v1_Prescl);
   fChain->SetBranchAddress("HLT_PADoubleJet20_ForwardBackward_v1", &HLT_PADoubleJet20_ForwardBackward_v1, &b_HLT_PADoubleJet20_ForwardBackward_v1);
   fChain->SetBranchAddress("HLT_PADoubleJet20_ForwardBackward_v1_Prescl", &HLT_PADoubleJet20_ForwardBackward_v1_Prescl, &b_HLT_PADoubleJet20_ForwardBackward_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAMu7_Ele7_CaloIdT_CaloIsoVL_v1", &HLT_PAMu7_Ele7_CaloIdT_CaloIsoVL_v1, &b_HLT_PAMu7_Ele7_CaloIdT_CaloIsoVL_v1);
   fChain->SetBranchAddress("HLT_PAMu7_Ele7_CaloIdT_CaloIsoVL_v1_Prescl", &HLT_PAMu7_Ele7_CaloIdT_CaloIsoVL_v1_Prescl, &b_HLT_PAMu7_Ele7_CaloIdT_CaloIsoVL_v1_Prescl);
   fChain->SetBranchAddress("HLT_PARandom_v1", &HLT_PARandom_v1, &b_HLT_PARandom_v1);
   fChain->SetBranchAddress("HLT_PARandom_v1_Prescl", &HLT_PARandom_v1_Prescl, &b_HLT_PARandom_v1_Prescl);
   fChain->SetBranchAddress("DQM_FEDIntegrity_v11", &DQM_FEDIntegrity_v11, &b_DQM_FEDIntegrity_v11);
   fChain->SetBranchAddress("DQM_FEDIntegrity_v11_Prescl", &DQM_FEDIntegrity_v11_Prescl, &b_DQM_FEDIntegrity_v11_Prescl);
   fChain->SetBranchAddress("HLT_LogMonitor_v4", &HLT_LogMonitor_v4, &b_HLT_LogMonitor_v4);
   fChain->SetBranchAddress("HLT_LogMonitor_v4_Prescl", &HLT_LogMonitor_v4_Prescl, &b_HLT_LogMonitor_v4_Prescl);
   fChain->SetBranchAddress("HLTriggerFinalPath", &HLTriggerFinalPath, &b_HLTriggerFinalPath);
   fChain->SetBranchAddress("HLTriggerFinalPath_Prescl", &HLTriggerFinalPath_Prescl, &b_HLTriggerFinalPath_Prescl);
   fChain->SetBranchAddress("L1_AlwaysTrue", &L1_AlwaysTrue, &b_L1_AlwaysTrue);
   fChain->SetBranchAddress("L1_AlwaysTrue_Prescl", &L1_AlwaysTrue_Prescl, &b_L1_AlwaysTrue_Prescl);
   fChain->SetBranchAddress("L1_BeamGas_Hf_BptxMinusPostQuiet", &L1_BeamGas_Hf_BptxMinusPostQuiet, &b_L1_BeamGas_Hf_BptxMinusPostQuiet);
   fChain->SetBranchAddress("L1_BeamGas_Hf_BptxMinusPostQuiet_Prescl", &L1_BeamGas_Hf_BptxMinusPostQuiet_Prescl, &b_L1_BeamGas_Hf_BptxMinusPostQuiet_Prescl);
   fChain->SetBranchAddress("L1_BeamGas_Hf_BptxPlusPostQuiet", &L1_BeamGas_Hf_BptxPlusPostQuiet, &b_L1_BeamGas_Hf_BptxPlusPostQuiet);
   fChain->SetBranchAddress("L1_BeamGas_Hf_BptxPlusPostQuiet_Prescl", &L1_BeamGas_Hf_BptxPlusPostQuiet_Prescl, &b_L1_BeamGas_Hf_BptxPlusPostQuiet_Prescl);
   fChain->SetBranchAddress("L1_BeamHalo", &L1_BeamHalo, &b_L1_BeamHalo);
   fChain->SetBranchAddress("L1_BeamHalo_Prescl", &L1_BeamHalo_Prescl, &b_L1_BeamHalo_Prescl);
   fChain->SetBranchAddress("L1_BptxMinus_NotBptxPlus", &L1_BptxMinus_NotBptxPlus, &b_L1_BptxMinus_NotBptxPlus);
   fChain->SetBranchAddress("L1_BptxMinus_NotBptxPlus_Prescl", &L1_BptxMinus_NotBptxPlus_Prescl, &b_L1_BptxMinus_NotBptxPlus_Prescl);
   fChain->SetBranchAddress("L1_BptxPlus_NotBptxMinus", &L1_BptxPlus_NotBptxMinus, &b_L1_BptxPlus_NotBptxMinus);
   fChain->SetBranchAddress("L1_BptxPlus_NotBptxMinus_Prescl", &L1_BptxPlus_NotBptxMinus_Prescl, &b_L1_BptxPlus_NotBptxMinus_Prescl);
   fChain->SetBranchAddress("L1_BscHighMultiplicity_BptxAND", &L1_BscHighMultiplicity_BptxAND, &b_L1_BscHighMultiplicity_BptxAND);
   fChain->SetBranchAddress("L1_BscHighMultiplicity_BptxAND_Prescl", &L1_BscHighMultiplicity_BptxAND_Prescl, &b_L1_BscHighMultiplicity_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_BscMinBiasOR_BptxAND", &L1_BscMinBiasOR_BptxAND, &b_L1_BscMinBiasOR_BptxAND);
   fChain->SetBranchAddress("L1_BscMinBiasOR_BptxAND_Prescl", &L1_BscMinBiasOR_BptxAND_Prescl, &b_L1_BscMinBiasOR_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_BscMinBiasThreshold1_BptxAND", &L1_BscMinBiasThreshold1_BptxAND, &b_L1_BscMinBiasThreshold1_BptxAND);
   fChain->SetBranchAddress("L1_BscMinBiasThreshold1_BptxAND_Prescl", &L1_BscMinBiasThreshold1_BptxAND_Prescl, &b_L1_BscMinBiasThreshold1_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_BscMinBiasThreshold2_BptxAND", &L1_BscMinBiasThreshold2_BptxAND, &b_L1_BscMinBiasThreshold2_BptxAND);
   fChain->SetBranchAddress("L1_BscMinBiasThreshold2_BptxAND_Prescl", &L1_BscMinBiasThreshold2_BptxAND_Prescl, &b_L1_BscMinBiasThreshold2_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_CastorEm_NotHcalHfCoincidencePm", &L1_CastorEm_NotHcalHfCoincidencePm, &b_L1_CastorEm_NotHcalHfCoincidencePm);
   fChain->SetBranchAddress("L1_CastorEm_NotHcalHfCoincidencePm_Prescl", &L1_CastorEm_NotHcalHfCoincidencePm_Prescl, &b_L1_CastorEm_NotHcalHfCoincidencePm_Prescl);
   fChain->SetBranchAddress("L1_CastorEm_NotHcalHfSingleChannel", &L1_CastorEm_NotHcalHfSingleChannel, &b_L1_CastorEm_NotHcalHfSingleChannel);
   fChain->SetBranchAddress("L1_CastorEm_NotHcalHfSingleChannel_Prescl", &L1_CastorEm_NotHcalHfSingleChannel_Prescl, &b_L1_CastorEm_NotHcalHfSingleChannel_Prescl);
   fChain->SetBranchAddress("L1_CastorEm_TotemLowMultiplicity", &L1_CastorEm_TotemLowMultiplicity, &b_L1_CastorEm_TotemLowMultiplicity);
   fChain->SetBranchAddress("L1_CastorEm_TotemLowMultiplicity_Prescl", &L1_CastorEm_TotemLowMultiplicity_Prescl, &b_L1_CastorEm_TotemLowMultiplicity_Prescl);
   fChain->SetBranchAddress("L1_CastorTotalEnergy_TotemLowMultiplicity", &L1_CastorTotalEnergy_TotemLowMultiplicity, &b_L1_CastorTotalEnergy_TotemLowMultiplicity);
   fChain->SetBranchAddress("L1_CastorTotalEnergy_TotemLowMultiplicity_Prescl", &L1_CastorTotalEnergy_TotemLowMultiplicity_Prescl, &b_L1_CastorTotalEnergy_TotemLowMultiplicity_Prescl);
   fChain->SetBranchAddress("L1_DoubleEG3_FwdVeto", &L1_DoubleEG3_FwdVeto, &b_L1_DoubleEG3_FwdVeto);
   fChain->SetBranchAddress("L1_DoubleEG3_FwdVeto_Prescl", &L1_DoubleEG3_FwdVeto_Prescl, &b_L1_DoubleEG3_FwdVeto_Prescl);
   fChain->SetBranchAddress("L1_DoubleEG5", &L1_DoubleEG5, &b_L1_DoubleEG5);
   fChain->SetBranchAddress("L1_DoubleEG5_Prescl", &L1_DoubleEG5_Prescl, &b_L1_DoubleEG5_Prescl);
   fChain->SetBranchAddress("L1_DoubleEG5_TotemDiffractive", &L1_DoubleEG5_TotemDiffractive, &b_L1_DoubleEG5_TotemDiffractive);
   fChain->SetBranchAddress("L1_DoubleEG5_TotemDiffractive_Prescl", &L1_DoubleEG5_TotemDiffractive_Prescl, &b_L1_DoubleEG5_TotemDiffractive_Prescl);
   fChain->SetBranchAddress("L1_DoubleEG6_HTT100", &L1_DoubleEG6_HTT100, &b_L1_DoubleEG6_HTT100);
   fChain->SetBranchAddress("L1_DoubleEG6_HTT100_Prescl", &L1_DoubleEG6_HTT100_Prescl, &b_L1_DoubleEG6_HTT100_Prescl);
   fChain->SetBranchAddress("L1_DoubleEG6_HTT125", &L1_DoubleEG6_HTT125, &b_L1_DoubleEG6_HTT125);
   fChain->SetBranchAddress("L1_DoubleEG6_HTT125_Prescl", &L1_DoubleEG6_HTT125_Prescl, &b_L1_DoubleEG6_HTT125_Prescl);
   fChain->SetBranchAddress("L1_DoubleEG_13_7", &L1_DoubleEG_13_7, &b_L1_DoubleEG_13_7);
   fChain->SetBranchAddress("L1_DoubleEG_13_7_Prescl", &L1_DoubleEG_13_7_Prescl, &b_L1_DoubleEG_13_7_Prescl);
   fChain->SetBranchAddress("L1_DoubleForJet16_EtaOpp", &L1_DoubleForJet16_EtaOpp, &b_L1_DoubleForJet16_EtaOpp);
   fChain->SetBranchAddress("L1_DoubleForJet16_EtaOpp_Prescl", &L1_DoubleForJet16_EtaOpp_Prescl, &b_L1_DoubleForJet16_EtaOpp_Prescl);
   fChain->SetBranchAddress("L1_DoubleJet20", &L1_DoubleJet20, &b_L1_DoubleJet20);
   fChain->SetBranchAddress("L1_DoubleJet20_Prescl", &L1_DoubleJet20_Prescl, &b_L1_DoubleJet20_Prescl);
   fChain->SetBranchAddress("L1_DoubleJet20_TotemDiffractive", &L1_DoubleJet20_TotemDiffractive, &b_L1_DoubleJet20_TotemDiffractive);
   fChain->SetBranchAddress("L1_DoubleJet20_TotemDiffractive_Prescl", &L1_DoubleJet20_TotemDiffractive_Prescl, &b_L1_DoubleJet20_TotemDiffractive_Prescl);
   fChain->SetBranchAddress("L1_DoubleJet24_v1", &L1_DoubleJet24_v1, &b_L1_DoubleJet24_v1);
   fChain->SetBranchAddress("L1_DoubleJet24_v1_Prescl", &L1_DoubleJet24_v1_Prescl, &b_L1_DoubleJet24_v1_Prescl);
   fChain->SetBranchAddress("L1_DoubleJet36_Central", &L1_DoubleJet36_Central, &b_L1_DoubleJet36_Central);
   fChain->SetBranchAddress("L1_DoubleJet36_Central_Prescl", &L1_DoubleJet36_Central_Prescl, &b_L1_DoubleJet36_Central_Prescl);
   fChain->SetBranchAddress("L1_DoubleJet52_Central", &L1_DoubleJet52_Central, &b_L1_DoubleJet52_Central);
   fChain->SetBranchAddress("L1_DoubleJet52_Central_Prescl", &L1_DoubleJet52_Central_Prescl, &b_L1_DoubleJet52_Central_Prescl);
   fChain->SetBranchAddress("L1_DoubleJetC36_TotemDiffractive", &L1_DoubleJetC36_TotemDiffractive, &b_L1_DoubleJetC36_TotemDiffractive);
   fChain->SetBranchAddress("L1_DoubleJetC36_TotemDiffractive_Prescl", &L1_DoubleJetC36_TotemDiffractive_Prescl, &b_L1_DoubleJetC36_TotemDiffractive_Prescl);
   fChain->SetBranchAddress("L1_DoubleJetC44_ETM30", &L1_DoubleJetC44_ETM30, &b_L1_DoubleJetC44_ETM30);
   fChain->SetBranchAddress("L1_DoubleJetC44_ETM30_Prescl", &L1_DoubleJetC44_ETM30_Prescl, &b_L1_DoubleJetC44_ETM30_Prescl);
   fChain->SetBranchAddress("L1_DoubleJetC56", &L1_DoubleJetC56, &b_L1_DoubleJetC56);
   fChain->SetBranchAddress("L1_DoubleJetC56_Prescl", &L1_DoubleJetC56_Prescl, &b_L1_DoubleJetC56_Prescl);
   fChain->SetBranchAddress("L1_DoubleJetC56_Eta1p74_WdEta4", &L1_DoubleJetC56_Eta1p74_WdEta4, &b_L1_DoubleJetC56_Eta1p74_WdEta4);
   fChain->SetBranchAddress("L1_DoubleJetC56_Eta1p74_WdEta4_Prescl", &L1_DoubleJetC56_Eta1p74_WdEta4_Prescl, &b_L1_DoubleJetC56_Eta1p74_WdEta4_Prescl);
   fChain->SetBranchAddress("L1_DoubleMu0", &L1_DoubleMu0, &b_L1_DoubleMu0);
   fChain->SetBranchAddress("L1_DoubleMu0_Prescl", &L1_DoubleMu0_Prescl, &b_L1_DoubleMu0_Prescl);
   fChain->SetBranchAddress("L1_DoubleMu0_HighQ_EtaCuts", &L1_DoubleMu0_HighQ_EtaCuts, &b_L1_DoubleMu0_HighQ_EtaCuts);
   fChain->SetBranchAddress("L1_DoubleMu0_HighQ_EtaCuts_Prescl", &L1_DoubleMu0_HighQ_EtaCuts_Prescl, &b_L1_DoubleMu0_HighQ_EtaCuts_Prescl);
   fChain->SetBranchAddress("L1_DoubleMu3p5_EG5", &L1_DoubleMu3p5_EG5, &b_L1_DoubleMu3p5_EG5);
   fChain->SetBranchAddress("L1_DoubleMu3p5_EG5_Prescl", &L1_DoubleMu3p5_EG5_Prescl, &b_L1_DoubleMu3p5_EG5_Prescl);
   fChain->SetBranchAddress("L1_DoubleMu5_EG5", &L1_DoubleMu5_EG5, &b_L1_DoubleMu5_EG5);
   fChain->SetBranchAddress("L1_DoubleMu5_EG5_Prescl", &L1_DoubleMu5_EG5_Prescl, &b_L1_DoubleMu5_EG5_Prescl);
   fChain->SetBranchAddress("L1_DoubleMu5_TotemDiffractive", &L1_DoubleMu5_TotemDiffractive, &b_L1_DoubleMu5_TotemDiffractive);
   fChain->SetBranchAddress("L1_DoubleMu5_TotemDiffractive_Prescl", &L1_DoubleMu5_TotemDiffractive_Prescl, &b_L1_DoubleMu5_TotemDiffractive_Prescl);
   fChain->SetBranchAddress("L1_DoubleMu5_v1", &L1_DoubleMu5_v1, &b_L1_DoubleMu5_v1);
   fChain->SetBranchAddress("L1_DoubleMu5_v1_Prescl", &L1_DoubleMu5_v1_Prescl, &b_L1_DoubleMu5_v1_Prescl);
   fChain->SetBranchAddress("L1_DoubleMuOpen_BptxAND", &L1_DoubleMuOpen_BptxAND, &b_L1_DoubleMuOpen_BptxAND);
   fChain->SetBranchAddress("L1_DoubleMuOpen_BptxAND_Prescl", &L1_DoubleMuOpen_BptxAND_Prescl, &b_L1_DoubleMuOpen_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_DoubleMu_10_3p5", &L1_DoubleMu_10_3p5, &b_L1_DoubleMu_10_3p5);
   fChain->SetBranchAddress("L1_DoubleMu_10_3p5_Prescl", &L1_DoubleMu_10_3p5_Prescl, &b_L1_DoubleMu_10_3p5_Prescl);
   fChain->SetBranchAddress("L1_DoubleMu_10_Open", &L1_DoubleMu_10_Open, &b_L1_DoubleMu_10_Open);
   fChain->SetBranchAddress("L1_DoubleMu_10_Open_Prescl", &L1_DoubleMu_10_Open_Prescl, &b_L1_DoubleMu_10_Open_Prescl);
   fChain->SetBranchAddress("L1_DoubleMu_12_5", &L1_DoubleMu_12_5, &b_L1_DoubleMu_12_5);
   fChain->SetBranchAddress("L1_DoubleMu_12_5_Prescl", &L1_DoubleMu_12_5_Prescl, &b_L1_DoubleMu_12_5_Prescl);
   fChain->SetBranchAddress("L1_DoubleMu_3er_0er_HighQ_WdEta22", &L1_DoubleMu_3er_0er_HighQ_WdEta22, &b_L1_DoubleMu_3er_0er_HighQ_WdEta22);
   fChain->SetBranchAddress("L1_DoubleMu_3er_0er_HighQ_WdEta22_Prescl", &L1_DoubleMu_3er_0er_HighQ_WdEta22_Prescl, &b_L1_DoubleMu_3er_0er_HighQ_WdEta22_Prescl);
   fChain->SetBranchAddress("L1_DoubleMu_5er_0er_HighQ_WdEta22", &L1_DoubleMu_5er_0er_HighQ_WdEta22, &b_L1_DoubleMu_5er_0er_HighQ_WdEta22);
   fChain->SetBranchAddress("L1_DoubleMu_5er_0er_HighQ_WdEta22_Prescl", &L1_DoubleMu_5er_0er_HighQ_WdEta22_Prescl, &b_L1_DoubleMu_5er_0er_HighQ_WdEta22_Prescl);
   fChain->SetBranchAddress("L1_EG8_DoubleJetC20", &L1_EG8_DoubleJetC20, &b_L1_EG8_DoubleJetC20);
   fChain->SetBranchAddress("L1_EG8_DoubleJetC20_Prescl", &L1_EG8_DoubleJetC20_Prescl, &b_L1_EG8_DoubleJetC20_Prescl);
   fChain->SetBranchAddress("L1_ETM100", &L1_ETM100, &b_L1_ETM100);
   fChain->SetBranchAddress("L1_ETM100_Prescl", &L1_ETM100_Prescl, &b_L1_ETM100_Prescl);
   fChain->SetBranchAddress("L1_ETM30", &L1_ETM30, &b_L1_ETM30);
   fChain->SetBranchAddress("L1_ETM30_Prescl", &L1_ETM30_Prescl, &b_L1_ETM30_Prescl);
   fChain->SetBranchAddress("L1_ETM36", &L1_ETM36, &b_L1_ETM36);
   fChain->SetBranchAddress("L1_ETM36_Prescl", &L1_ETM36_Prescl, &b_L1_ETM36_Prescl);
   fChain->SetBranchAddress("L1_ETM40", &L1_ETM40, &b_L1_ETM40);
   fChain->SetBranchAddress("L1_ETM40_Prescl", &L1_ETM40_Prescl, &b_L1_ETM40_Prescl);
   fChain->SetBranchAddress("L1_ETM50", &L1_ETM50, &b_L1_ETM50);
   fChain->SetBranchAddress("L1_ETM50_Prescl", &L1_ETM50_Prescl, &b_L1_ETM50_Prescl);
   fChain->SetBranchAddress("L1_ETM70", &L1_ETM70, &b_L1_ETM70);
   fChain->SetBranchAddress("L1_ETM70_Prescl", &L1_ETM70_Prescl, &b_L1_ETM70_Prescl);
   fChain->SetBranchAddress("L1_ETT140", &L1_ETT140, &b_L1_ETT140);
   fChain->SetBranchAddress("L1_ETT140_Prescl", &L1_ETT140_Prescl, &b_L1_ETT140_Prescl);
   fChain->SetBranchAddress("L1_ETT20_BptxAND", &L1_ETT20_BptxAND, &b_L1_ETT20_BptxAND);
   fChain->SetBranchAddress("L1_ETT20_BptxAND_Prescl", &L1_ETT20_BptxAND_Prescl, &b_L1_ETT20_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ETT300", &L1_ETT300, &b_L1_ETT300);
   fChain->SetBranchAddress("L1_ETT300_Prescl", &L1_ETT300_Prescl, &b_L1_ETT300_Prescl);
   fChain->SetBranchAddress("L1_ETT40", &L1_ETT40, &b_L1_ETT40);
   fChain->SetBranchAddress("L1_ETT40_Prescl", &L1_ETT40_Prescl, &b_L1_ETT40_Prescl);
   fChain->SetBranchAddress("L1_ETT60", &L1_ETT60, &b_L1_ETT60);
   fChain->SetBranchAddress("L1_ETT60_Prescl", &L1_ETT60_Prescl, &b_L1_ETT60_Prescl);
   fChain->SetBranchAddress("L1_ETT80", &L1_ETT80, &b_L1_ETT80);
   fChain->SetBranchAddress("L1_ETT80_Prescl", &L1_ETT80_Prescl, &b_L1_ETT80_Prescl);
   fChain->SetBranchAddress("L1_HTT100", &L1_HTT100, &b_L1_HTT100);
   fChain->SetBranchAddress("L1_HTT100_Prescl", &L1_HTT100_Prescl, &b_L1_HTT100_Prescl);
   fChain->SetBranchAddress("L1_HTT125", &L1_HTT125, &b_L1_HTT125);
   fChain->SetBranchAddress("L1_HTT125_Prescl", &L1_HTT125_Prescl, &b_L1_HTT125_Prescl);
   fChain->SetBranchAddress("L1_HTT150", &L1_HTT150, &b_L1_HTT150);
   fChain->SetBranchAddress("L1_HTT150_Prescl", &L1_HTT150_Prescl, &b_L1_HTT150_Prescl);
   fChain->SetBranchAddress("L1_HTT175", &L1_HTT175, &b_L1_HTT175);
   fChain->SetBranchAddress("L1_HTT175_Prescl", &L1_HTT175_Prescl, &b_L1_HTT175_Prescl);
   fChain->SetBranchAddress("L1_HTT200", &L1_HTT200, &b_L1_HTT200);
   fChain->SetBranchAddress("L1_HTT200_Prescl", &L1_HTT200_Prescl, &b_L1_HTT200_Prescl);
   fChain->SetBranchAddress("L1_HTT50", &L1_HTT50, &b_L1_HTT50);
   fChain->SetBranchAddress("L1_HTT50_Prescl", &L1_HTT50_Prescl, &b_L1_HTT50_Prescl);
   fChain->SetBranchAddress("L1_HTT75", &L1_HTT75, &b_L1_HTT75);
   fChain->SetBranchAddress("L1_HTT75_Prescl", &L1_HTT75_Prescl, &b_L1_HTT75_Prescl);
   fChain->SetBranchAddress("L1_HcalHfCoincidencePm_BptxAND_v1", &L1_HcalHfCoincidencePm_BptxAND_v1, &b_L1_HcalHfCoincidencePm_BptxAND_v1);
   fChain->SetBranchAddress("L1_HcalHfCoincidencePm_BptxAND_v1_Prescl", &L1_HcalHfCoincidencePm_BptxAND_v1_Prescl, &b_L1_HcalHfCoincidencePm_BptxAND_v1_Prescl);
   fChain->SetBranchAddress("L1_HcalHfSingleChannel_BptxAND", &L1_HcalHfSingleChannel_BptxAND, &b_L1_HcalHfSingleChannel_BptxAND);
   fChain->SetBranchAddress("L1_HcalHfSingleChannel_BptxAND_Prescl", &L1_HcalHfSingleChannel_BptxAND_Prescl, &b_L1_HcalHfSingleChannel_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_HcalHfSingleChannel_BptxAND_Instance1", &L1_HcalHfSingleChannel_BptxAND_Instance1, &b_L1_HcalHfSingleChannel_BptxAND_Instance1);
   fChain->SetBranchAddress("L1_HcalHfSingleChannel_BptxAND_Instance1_Prescl", &L1_HcalHfSingleChannel_BptxAND_Instance1_Prescl, &b_L1_HcalHfSingleChannel_BptxAND_Instance1_Prescl);
   fChain->SetBranchAddress("L1_Mu10er_JetC32", &L1_Mu10er_JetC32, &b_L1_Mu10er_JetC32);
   fChain->SetBranchAddress("L1_Mu10er_JetC32_Prescl", &L1_Mu10er_JetC32_Prescl, &b_L1_Mu10er_JetC32_Prescl);
   fChain->SetBranchAddress("L1_Mu12_EG7", &L1_Mu12_EG7, &b_L1_Mu12_EG7);
   fChain->SetBranchAddress("L1_Mu12_EG7_Prescl", &L1_Mu12_EG7_Prescl, &b_L1_Mu12_EG7_Prescl);
   fChain->SetBranchAddress("L1_Mu3_Jet16", &L1_Mu3_Jet16, &b_L1_Mu3_Jet16);
   fChain->SetBranchAddress("L1_Mu3_Jet16_Prescl", &L1_Mu3_Jet16_Prescl, &b_L1_Mu3_Jet16_Prescl);
   fChain->SetBranchAddress("L1_Mu3_Jet36", &L1_Mu3_Jet36, &b_L1_Mu3_Jet36);
   fChain->SetBranchAddress("L1_Mu3_Jet36_Prescl", &L1_Mu3_Jet36_Prescl, &b_L1_Mu3_Jet36_Prescl);
   fChain->SetBranchAddress("L1_Mu3_JetC16_WdEtaPhi2", &L1_Mu3_JetC16_WdEtaPhi2, &b_L1_Mu3_JetC16_WdEtaPhi2);
   fChain->SetBranchAddress("L1_Mu3_JetC16_WdEtaPhi2_Prescl", &L1_Mu3_JetC16_WdEtaPhi2_Prescl, &b_L1_Mu3_JetC16_WdEtaPhi2_Prescl);
   fChain->SetBranchAddress("L1_Mu3_JetC52_WdEtaPhi2", &L1_Mu3_JetC52_WdEtaPhi2, &b_L1_Mu3_JetC52_WdEtaPhi2);
   fChain->SetBranchAddress("L1_Mu3_JetC52_WdEtaPhi2_Prescl", &L1_Mu3_JetC52_WdEtaPhi2_Prescl, &b_L1_Mu3_JetC52_WdEtaPhi2_Prescl);
   fChain->SetBranchAddress("L1_Mu5_DoubleEG5", &L1_Mu5_DoubleEG5, &b_L1_Mu5_DoubleEG5);
   fChain->SetBranchAddress("L1_Mu5_DoubleEG5_Prescl", &L1_Mu5_DoubleEG5_Prescl, &b_L1_Mu5_DoubleEG5_Prescl);
   fChain->SetBranchAddress("L1_Mu5_DoubleEG6", &L1_Mu5_DoubleEG6, &b_L1_Mu5_DoubleEG6);
   fChain->SetBranchAddress("L1_Mu5_DoubleEG6_Prescl", &L1_Mu5_DoubleEG6_Prescl, &b_L1_Mu5_DoubleEG6_Prescl);
   fChain->SetBranchAddress("L1_Mu7_Jet16", &L1_Mu7_Jet16, &b_L1_Mu7_Jet16);
   fChain->SetBranchAddress("L1_Mu7_Jet16_Prescl", &L1_Mu7_Jet16_Prescl, &b_L1_Mu7_Jet16_Prescl);
   fChain->SetBranchAddress("L1_Mu8_DoubleJetC20", &L1_Mu8_DoubleJetC20, &b_L1_Mu8_DoubleJetC20);
   fChain->SetBranchAddress("L1_Mu8_DoubleJetC20_Prescl", &L1_Mu8_DoubleJetC20_Prescl, &b_L1_Mu8_DoubleJetC20_Prescl);
   fChain->SetBranchAddress("L1_MuOpen_EG12", &L1_MuOpen_EG12, &b_L1_MuOpen_EG12);
   fChain->SetBranchAddress("L1_MuOpen_EG12_Prescl", &L1_MuOpen_EG12_Prescl, &b_L1_MuOpen_EG12_Prescl);
   fChain->SetBranchAddress("L1_MuOpen_EG5", &L1_MuOpen_EG5, &b_L1_MuOpen_EG5);
   fChain->SetBranchAddress("L1_MuOpen_EG5_Prescl", &L1_MuOpen_EG5_Prescl, &b_L1_MuOpen_EG5_Prescl);
   fChain->SetBranchAddress("L1_QuadJetC32", &L1_QuadJetC32, &b_L1_QuadJetC32);
   fChain->SetBranchAddress("L1_QuadJetC32_Prescl", &L1_QuadJetC32_Prescl, &b_L1_QuadJetC32_Prescl);
   fChain->SetBranchAddress("L1_QuadJetC36", &L1_QuadJetC36, &b_L1_QuadJetC36);
   fChain->SetBranchAddress("L1_QuadJetC36_Prescl", &L1_QuadJetC36_Prescl, &b_L1_QuadJetC36_Prescl);
   fChain->SetBranchAddress("L1_QuadJetC40", &L1_QuadJetC40, &b_L1_QuadJetC40);
   fChain->SetBranchAddress("L1_QuadJetC40_Prescl", &L1_QuadJetC40_Prescl, &b_L1_QuadJetC40_Prescl);
   fChain->SetBranchAddress("L1_SingleEG12", &L1_SingleEG12, &b_L1_SingleEG12);
   fChain->SetBranchAddress("L1_SingleEG12_Prescl", &L1_SingleEG12_Prescl, &b_L1_SingleEG12_Prescl);
   fChain->SetBranchAddress("L1_SingleEG18er", &L1_SingleEG18er, &b_L1_SingleEG18er);
   fChain->SetBranchAddress("L1_SingleEG18er_Prescl", &L1_SingleEG18er_Prescl, &b_L1_SingleEG18er_Prescl);
   fChain->SetBranchAddress("L1_SingleEG20", &L1_SingleEG20, &b_L1_SingleEG20);
   fChain->SetBranchAddress("L1_SingleEG20_Prescl", &L1_SingleEG20_Prescl, &b_L1_SingleEG20_Prescl);
   fChain->SetBranchAddress("L1_SingleEG20_TotemDiffractive", &L1_SingleEG20_TotemDiffractive, &b_L1_SingleEG20_TotemDiffractive);
   fChain->SetBranchAddress("L1_SingleEG20_TotemDiffractive_Prescl", &L1_SingleEG20_TotemDiffractive_Prescl, &b_L1_SingleEG20_TotemDiffractive_Prescl);
   fChain->SetBranchAddress("L1_SingleEG22", &L1_SingleEG22, &b_L1_SingleEG22);
   fChain->SetBranchAddress("L1_SingleEG22_Prescl", &L1_SingleEG22_Prescl, &b_L1_SingleEG22_Prescl);
   fChain->SetBranchAddress("L1_SingleEG24", &L1_SingleEG24, &b_L1_SingleEG24);
   fChain->SetBranchAddress("L1_SingleEG24_Prescl", &L1_SingleEG24_Prescl, &b_L1_SingleEG24_Prescl);
   fChain->SetBranchAddress("L1_SingleEG30", &L1_SingleEG30, &b_L1_SingleEG30);
   fChain->SetBranchAddress("L1_SingleEG30_Prescl", &L1_SingleEG30_Prescl, &b_L1_SingleEG30_Prescl);
   fChain->SetBranchAddress("L1_SingleEG5_BptxAND", &L1_SingleEG5_BptxAND, &b_L1_SingleEG5_BptxAND);
   fChain->SetBranchAddress("L1_SingleEG5_BptxAND_Prescl", &L1_SingleEG5_BptxAND_Prescl, &b_L1_SingleEG5_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleEG7", &L1_SingleEG7, &b_L1_SingleEG7);
   fChain->SetBranchAddress("L1_SingleEG7_Prescl", &L1_SingleEG7_Prescl, &b_L1_SingleEG7_Prescl);
   fChain->SetBranchAddress("L1_SingleForJet16", &L1_SingleForJet16, &b_L1_SingleForJet16);
   fChain->SetBranchAddress("L1_SingleForJet16_Prescl", &L1_SingleForJet16_Prescl, &b_L1_SingleForJet16_Prescl);
   fChain->SetBranchAddress("L1_SingleIsoEG18er", &L1_SingleIsoEG18er, &b_L1_SingleIsoEG18er);
   fChain->SetBranchAddress("L1_SingleIsoEG18er_Prescl", &L1_SingleIsoEG18er_Prescl, &b_L1_SingleIsoEG18er_Prescl);
   fChain->SetBranchAddress("L1_SingleIsoEG20er", &L1_SingleIsoEG20er, &b_L1_SingleIsoEG20er);
   fChain->SetBranchAddress("L1_SingleIsoEG20er_Prescl", &L1_SingleIsoEG20er_Prescl, &b_L1_SingleIsoEG20er_Prescl);
   fChain->SetBranchAddress("L1_SingleJet128", &L1_SingleJet128, &b_L1_SingleJet128);
   fChain->SetBranchAddress("L1_SingleJet128_Prescl", &L1_SingleJet128_Prescl, &b_L1_SingleJet128_Prescl);
   fChain->SetBranchAddress("L1_SingleJet12_BptxAND", &L1_SingleJet12_BptxAND, &b_L1_SingleJet12_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet12_BptxAND_Prescl", &L1_SingleJet12_BptxAND_Prescl, &b_L1_SingleJet12_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet16_BptxAND", &L1_SingleJet16_BptxAND, &b_L1_SingleJet16_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet16_BptxAND_Prescl", &L1_SingleJet16_BptxAND_Prescl, &b_L1_SingleJet16_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet16_FwdVeto5", &L1_SingleJet16_FwdVeto5, &b_L1_SingleJet16_FwdVeto5);
   fChain->SetBranchAddress("L1_SingleJet16_FwdVeto5_Prescl", &L1_SingleJet16_FwdVeto5_Prescl, &b_L1_SingleJet16_FwdVeto5_Prescl);
   fChain->SetBranchAddress("L1_SingleJet20_Central_NotBptxOR", &L1_SingleJet20_Central_NotBptxOR, &b_L1_SingleJet20_Central_NotBptxOR);
   fChain->SetBranchAddress("L1_SingleJet20_Central_NotBptxOR_Prescl", &L1_SingleJet20_Central_NotBptxOR_Prescl, &b_L1_SingleJet20_Central_NotBptxOR_Prescl);
   fChain->SetBranchAddress("L1_SingleJet36", &L1_SingleJet36, &b_L1_SingleJet36);
   fChain->SetBranchAddress("L1_SingleJet36_Prescl", &L1_SingleJet36_Prescl, &b_L1_SingleJet36_Prescl);
   fChain->SetBranchAddress("L1_SingleJet36_FwdVeto5", &L1_SingleJet36_FwdVeto5, &b_L1_SingleJet36_FwdVeto5);
   fChain->SetBranchAddress("L1_SingleJet36_FwdVeto5_Prescl", &L1_SingleJet36_FwdVeto5_Prescl, &b_L1_SingleJet36_FwdVeto5_Prescl);
   fChain->SetBranchAddress("L1_SingleJet52", &L1_SingleJet52, &b_L1_SingleJet52);
   fChain->SetBranchAddress("L1_SingleJet52_Prescl", &L1_SingleJet52_Prescl, &b_L1_SingleJet52_Prescl);
   fChain->SetBranchAddress("L1_SingleJet52_TotemDiffractive", &L1_SingleJet52_TotemDiffractive, &b_L1_SingleJet52_TotemDiffractive);
   fChain->SetBranchAddress("L1_SingleJet52_TotemDiffractive_Prescl", &L1_SingleJet52_TotemDiffractive_Prescl, &b_L1_SingleJet52_TotemDiffractive_Prescl);
   fChain->SetBranchAddress("L1_SingleJet68", &L1_SingleJet68, &b_L1_SingleJet68);
   fChain->SetBranchAddress("L1_SingleJet68_Prescl", &L1_SingleJet68_Prescl, &b_L1_SingleJet68_Prescl);
   fChain->SetBranchAddress("L1_SingleJet92", &L1_SingleJet92, &b_L1_SingleJet92);
   fChain->SetBranchAddress("L1_SingleJet92_Prescl", &L1_SingleJet92_Prescl, &b_L1_SingleJet92_Prescl);
   fChain->SetBranchAddress("L1_SingleJetC32_NotBptxOR", &L1_SingleJetC32_NotBptxOR, &b_L1_SingleJetC32_NotBptxOR);
   fChain->SetBranchAddress("L1_SingleJetC32_NotBptxOR_Prescl", &L1_SingleJetC32_NotBptxOR_Prescl, &b_L1_SingleJetC32_NotBptxOR_Prescl);
   fChain->SetBranchAddress("L1_SingleMu12", &L1_SingleMu12, &b_L1_SingleMu12);
   fChain->SetBranchAddress("L1_SingleMu12_Prescl", &L1_SingleMu12_Prescl, &b_L1_SingleMu12_Prescl);
   fChain->SetBranchAddress("L1_SingleMu12er", &L1_SingleMu12er, &b_L1_SingleMu12er);
   fChain->SetBranchAddress("L1_SingleMu12er_Prescl", &L1_SingleMu12er_Prescl, &b_L1_SingleMu12er_Prescl);
   fChain->SetBranchAddress("L1_SingleMu14_Eta2p1", &L1_SingleMu14_Eta2p1, &b_L1_SingleMu14_Eta2p1);
   fChain->SetBranchAddress("L1_SingleMu14_Eta2p1_Prescl", &L1_SingleMu14_Eta2p1_Prescl, &b_L1_SingleMu14_Eta2p1_Prescl);
   fChain->SetBranchAddress("L1_SingleMu16", &L1_SingleMu16, &b_L1_SingleMu16);
   fChain->SetBranchAddress("L1_SingleMu16_Prescl", &L1_SingleMu16_Prescl, &b_L1_SingleMu16_Prescl);
   fChain->SetBranchAddress("L1_SingleMu16_Eta2p1", &L1_SingleMu16_Eta2p1, &b_L1_SingleMu16_Eta2p1);
   fChain->SetBranchAddress("L1_SingleMu16_Eta2p1_Prescl", &L1_SingleMu16_Eta2p1_Prescl, &b_L1_SingleMu16_Eta2p1_Prescl);
   fChain->SetBranchAddress("L1_SingleMu18er", &L1_SingleMu18er, &b_L1_SingleMu18er);
   fChain->SetBranchAddress("L1_SingleMu18er_Prescl", &L1_SingleMu18er_Prescl, &b_L1_SingleMu18er_Prescl);
   fChain->SetBranchAddress("L1_SingleMu20", &L1_SingleMu20, &b_L1_SingleMu20);
   fChain->SetBranchAddress("L1_SingleMu20_Prescl", &L1_SingleMu20_Prescl, &b_L1_SingleMu20_Prescl);
   fChain->SetBranchAddress("L1_SingleMu20_TotemDiffractive", &L1_SingleMu20_TotemDiffractive, &b_L1_SingleMu20_TotemDiffractive);
   fChain->SetBranchAddress("L1_SingleMu20_TotemDiffractive_Prescl", &L1_SingleMu20_TotemDiffractive_Prescl, &b_L1_SingleMu20_TotemDiffractive_Prescl);
   fChain->SetBranchAddress("L1_SingleMu20er", &L1_SingleMu20er, &b_L1_SingleMu20er);
   fChain->SetBranchAddress("L1_SingleMu20er_Prescl", &L1_SingleMu20er_Prescl, &b_L1_SingleMu20er_Prescl);
   fChain->SetBranchAddress("L1_SingleMu25er", &L1_SingleMu25er, &b_L1_SingleMu25er);
   fChain->SetBranchAddress("L1_SingleMu25er_Prescl", &L1_SingleMu25er_Prescl, &b_L1_SingleMu25er_Prescl);
   fChain->SetBranchAddress("L1_SingleMu3", &L1_SingleMu3, &b_L1_SingleMu3);
   fChain->SetBranchAddress("L1_SingleMu3_Prescl", &L1_SingleMu3_Prescl, &b_L1_SingleMu3_Prescl);
   fChain->SetBranchAddress("L1_SingleMu6_NotBptxOR", &L1_SingleMu6_NotBptxOR, &b_L1_SingleMu6_NotBptxOR);
   fChain->SetBranchAddress("L1_SingleMu6_NotBptxOR_Prescl", &L1_SingleMu6_NotBptxOR_Prescl, &b_L1_SingleMu6_NotBptxOR_Prescl);
   fChain->SetBranchAddress("L1_SingleMu7", &L1_SingleMu7, &b_L1_SingleMu7);
   fChain->SetBranchAddress("L1_SingleMu7_Prescl", &L1_SingleMu7_Prescl, &b_L1_SingleMu7_Prescl);
   fChain->SetBranchAddress("L1_SingleMuBeamHalo", &L1_SingleMuBeamHalo, &b_L1_SingleMuBeamHalo);
   fChain->SetBranchAddress("L1_SingleMuBeamHalo_Prescl", &L1_SingleMuBeamHalo_Prescl, &b_L1_SingleMuBeamHalo_Prescl);
   fChain->SetBranchAddress("L1_SingleMuOpen", &L1_SingleMuOpen, &b_L1_SingleMuOpen);
   fChain->SetBranchAddress("L1_SingleMuOpen_Prescl", &L1_SingleMuOpen_Prescl, &b_L1_SingleMuOpen_Prescl);
   fChain->SetBranchAddress("L1_TotemMinBias", &L1_TotemMinBias, &b_L1_TotemMinBias);
   fChain->SetBranchAddress("L1_TotemMinBias_Prescl", &L1_TotemMinBias_Prescl, &b_L1_TotemMinBias_Prescl);
   fChain->SetBranchAddress("L1_TripleEG7", &L1_TripleEG7, &b_L1_TripleEG7);
   fChain->SetBranchAddress("L1_TripleEG7_Prescl", &L1_TripleEG7_Prescl, &b_L1_TripleEG7_Prescl);
   fChain->SetBranchAddress("L1_TripleEG_12_7_5", &L1_TripleEG_12_7_5, &b_L1_TripleEG_12_7_5);
   fChain->SetBranchAddress("L1_TripleEG_12_7_5_Prescl", &L1_TripleEG_12_7_5_Prescl, &b_L1_TripleEG_12_7_5_Prescl);
   fChain->SetBranchAddress("L1_TripleJetC_52_28_28", &L1_TripleJetC_52_28_28, &b_L1_TripleJetC_52_28_28);
   fChain->SetBranchAddress("L1_TripleJetC_52_28_28_Prescl", &L1_TripleJetC_52_28_28_Prescl, &b_L1_TripleJetC_52_28_28_Prescl);
   fChain->SetBranchAddress("L1_TripleJet_64_44_24_VBF", &L1_TripleJet_64_44_24_VBF, &b_L1_TripleJet_64_44_24_VBF);
   fChain->SetBranchAddress("L1_TripleJet_64_44_24_VBF_Prescl", &L1_TripleJet_64_44_24_VBF_Prescl, &b_L1_TripleJet_64_44_24_VBF_Prescl);
   fChain->SetBranchAddress("L1_TripleJet_64_48_28_VBF", &L1_TripleJet_64_48_28_VBF, &b_L1_TripleJet_64_48_28_VBF);
   fChain->SetBranchAddress("L1_TripleJet_64_48_28_VBF_Prescl", &L1_TripleJet_64_48_28_VBF_Prescl, &b_L1_TripleJet_64_48_28_VBF_Prescl);
   fChain->SetBranchAddress("L1_ZdcCaloPlus_TotemDiffractive_QElastic", &L1_ZdcCaloPlus_TotemDiffractive_QElastic, &b_L1_ZdcCaloPlus_TotemDiffractive_QElastic);
   fChain->SetBranchAddress("L1_ZdcCaloPlus_TotemDiffractive_QElastic_Prescl", &L1_ZdcCaloPlus_TotemDiffractive_QElastic_Prescl, &b_L1_ZdcCaloPlus_TotemDiffractive_QElastic_Prescl);
   fChain->SetBranchAddress("L1_ZeroBias", &L1_ZeroBias, &b_L1_ZeroBias);
   fChain->SetBranchAddress("L1_ZeroBias_Prescl", &L1_ZeroBias_Prescl, &b_L1_ZeroBias_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_PreBPTX.v0", &L1Tech_BPTX_PreBPTX_v0, &b_L1Tech_BPTX_PreBPTX_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_PreBPTX.v0_Prescl", &L1Tech_BPTX_PreBPTX_v0_Prescl, &b_L1Tech_BPTX_PreBPTX_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_minus.v0", &L1Tech_BPTX_minus_v0, &b_L1Tech_BPTX_minus_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_minus.v0_Prescl", &L1Tech_BPTX_minus_v0_Prescl, &b_L1Tech_BPTX_minus_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_minus_AND_not_plus.v0", &L1Tech_BPTX_minus_AND_not_plus_v0, &b_L1Tech_BPTX_minus_AND_not_plus_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_minus_AND_not_plus.v0_Prescl", &L1Tech_BPTX_minus_AND_not_plus_v0_Prescl, &b_L1Tech_BPTX_minus_AND_not_plus_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_plus.v0", &L1Tech_BPTX_plus_v0, &b_L1Tech_BPTX_plus_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_plus.v0_Prescl", &L1Tech_BPTX_plus_v0_Prescl, &b_L1Tech_BPTX_plus_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_AND_NOT_minus.v0", &L1Tech_BPTX_plus_AND_NOT_minus_v0, &b_L1Tech_BPTX_plus_AND_NOT_minus_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_AND_NOT_minus.v0_Prescl", &L1Tech_BPTX_plus_AND_NOT_minus_v0_Prescl, &b_L1Tech_BPTX_plus_AND_NOT_minus_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_AND_minus.v0", &L1Tech_BPTX_plus_AND_minus_v0, &b_L1Tech_BPTX_plus_AND_minus_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_AND_minus.v0_Prescl", &L1Tech_BPTX_plus_AND_minus_v0_Prescl, &b_L1Tech_BPTX_plus_AND_minus_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_AND_minus_instance1.v0", &L1Tech_BPTX_plus_AND_minus_instance1_v0, &b_L1Tech_BPTX_plus_AND_minus_instance1_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_AND_minus_instance1.v0_Prescl", &L1Tech_BPTX_plus_AND_minus_instance1_v0_Prescl, &b_L1Tech_BPTX_plus_AND_minus_instance1_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_OR_minus.v0", &L1Tech_BPTX_plus_OR_minus_v0, &b_L1Tech_BPTX_plus_OR_minus_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_OR_minus.v0_Prescl", &L1Tech_BPTX_plus_OR_minus_v0_Prescl, &b_L1Tech_BPTX_plus_OR_minus_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_quiet.v0", &L1Tech_BPTX_quiet_v0, &b_L1Tech_BPTX_quiet_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_quiet.v0_Prescl", &L1Tech_BPTX_quiet_v0_Prescl, &b_L1Tech_BPTX_quiet_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_HighMultiplicity.v0", &L1Tech_BSC_HighMultiplicity_v0, &b_L1Tech_BSC_HighMultiplicity_v0);
   fChain->SetBranchAddress("L1Tech_BSC_HighMultiplicity.v0_Prescl", &L1Tech_BSC_HighMultiplicity_v0_Prescl, &b_L1Tech_BSC_HighMultiplicity_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_halo_beam1_inner.v0", &L1Tech_BSC_halo_beam1_inner_v0, &b_L1Tech_BSC_halo_beam1_inner_v0);
   fChain->SetBranchAddress("L1Tech_BSC_halo_beam1_inner.v0_Prescl", &L1Tech_BSC_halo_beam1_inner_v0_Prescl, &b_L1Tech_BSC_halo_beam1_inner_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_halo_beam1_outer.v0", &L1Tech_BSC_halo_beam1_outer_v0, &b_L1Tech_BSC_halo_beam1_outer_v0);
   fChain->SetBranchAddress("L1Tech_BSC_halo_beam1_outer.v0_Prescl", &L1Tech_BSC_halo_beam1_outer_v0_Prescl, &b_L1Tech_BSC_halo_beam1_outer_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_halo_beam2_inner.v0", &L1Tech_BSC_halo_beam2_inner_v0, &b_L1Tech_BSC_halo_beam2_inner_v0);
   fChain->SetBranchAddress("L1Tech_BSC_halo_beam2_inner.v0_Prescl", &L1Tech_BSC_halo_beam2_inner_v0_Prescl, &b_L1Tech_BSC_halo_beam2_inner_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_halo_beam2_outer.v0", &L1Tech_BSC_halo_beam2_outer_v0, &b_L1Tech_BSC_halo_beam2_outer_v0);
   fChain->SetBranchAddress("L1Tech_BSC_halo_beam2_outer.v0_Prescl", &L1Tech_BSC_halo_beam2_outer_v0_Prescl, &b_L1Tech_BSC_halo_beam2_outer_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_OR.v0", &L1Tech_BSC_minBias_OR_v0, &b_L1Tech_BSC_minBias_OR_v0);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_OR.v0_Prescl", &L1Tech_BSC_minBias_OR_v0_Prescl, &b_L1Tech_BSC_minBias_OR_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_inner_threshold1.v0", &L1Tech_BSC_minBias_inner_threshold1_v0, &b_L1Tech_BSC_minBias_inner_threshold1_v0);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_inner_threshold1.v0_Prescl", &L1Tech_BSC_minBias_inner_threshold1_v0_Prescl, &b_L1Tech_BSC_minBias_inner_threshold1_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_inner_threshold2.v0", &L1Tech_BSC_minBias_inner_threshold2_v0, &b_L1Tech_BSC_minBias_inner_threshold2_v0);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_inner_threshold2.v0_Prescl", &L1Tech_BSC_minBias_inner_threshold2_v0_Prescl, &b_L1Tech_BSC_minBias_inner_threshold2_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_threshold1.v0", &L1Tech_BSC_minBias_threshold1_v0, &b_L1Tech_BSC_minBias_threshold1_v0);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_threshold1.v0_Prescl", &L1Tech_BSC_minBias_threshold1_v0_Prescl, &b_L1Tech_BSC_minBias_threshold1_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_threshold2.v0", &L1Tech_BSC_minBias_threshold2_v0, &b_L1Tech_BSC_minBias_threshold2_v0);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_threshold2.v0_Prescl", &L1Tech_BSC_minBias_threshold2_v0_Prescl, &b_L1Tech_BSC_minBias_threshold2_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_splash_beam1.v0", &L1Tech_BSC_splash_beam1_v0, &b_L1Tech_BSC_splash_beam1_v0);
   fChain->SetBranchAddress("L1Tech_BSC_splash_beam1.v0_Prescl", &L1Tech_BSC_splash_beam1_v0_Prescl, &b_L1Tech_BSC_splash_beam1_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_splash_beam2.v0", &L1Tech_BSC_splash_beam2_v0, &b_L1Tech_BSC_splash_beam2_v0);
   fChain->SetBranchAddress("L1Tech_BSC_splash_beam2.v0_Prescl", &L1Tech_BSC_splash_beam2_v0_Prescl, &b_L1Tech_BSC_splash_beam2_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_CASTOR_0.v0", &L1Tech_CASTOR_0_v0, &b_L1Tech_CASTOR_0_v0);
   fChain->SetBranchAddress("L1Tech_CASTOR_0.v0_Prescl", &L1Tech_CASTOR_0_v0_Prescl, &b_L1Tech_CASTOR_0_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_CASTOR_EM.v0", &L1Tech_CASTOR_EM_v0, &b_L1Tech_CASTOR_EM_v0);
   fChain->SetBranchAddress("L1Tech_CASTOR_EM.v0_Prescl", &L1Tech_CASTOR_EM_v0_Prescl, &b_L1Tech_CASTOR_EM_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_CASTOR_HaloMuon.v0", &L1Tech_CASTOR_HaloMuon_v0, &b_L1Tech_CASTOR_HaloMuon_v0);
   fChain->SetBranchAddress("L1Tech_CASTOR_HaloMuon.v0_Prescl", &L1Tech_CASTOR_HaloMuon_v0_Prescl, &b_L1Tech_CASTOR_HaloMuon_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_CASTOR_TotalEnergy.v0", &L1Tech_CASTOR_TotalEnergy_v0, &b_L1Tech_CASTOR_TotalEnergy_v0);
   fChain->SetBranchAddress("L1Tech_CASTOR_TotalEnergy.v0_Prescl", &L1Tech_CASTOR_TotalEnergy_v0_Prescl, &b_L1Tech_CASTOR_TotalEnergy_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_DT_GlobalOR.v0", &L1Tech_DT_GlobalOR_v0, &b_L1Tech_DT_GlobalOR_v0);
   fChain->SetBranchAddress("L1Tech_DT_GlobalOR.v0_Prescl", &L1Tech_DT_GlobalOR_v0_Prescl, &b_L1Tech_DT_GlobalOR_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_FSC_St3Sect45_downLeft.v0", &L1Tech_FSC_St3Sect45_downLeft_v0, &b_L1Tech_FSC_St3Sect45_downLeft_v0);
   fChain->SetBranchAddress("L1Tech_FSC_St3Sect45_downLeft.v0_Prescl", &L1Tech_FSC_St3Sect45_downLeft_v0_Prescl, &b_L1Tech_FSC_St3Sect45_downLeft_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_FSC_St3Sect45_downRight.v0", &L1Tech_FSC_St3Sect45_downRight_v0, &b_L1Tech_FSC_St3Sect45_downRight_v0);
   fChain->SetBranchAddress("L1Tech_FSC_St3Sect45_downRight.v0_Prescl", &L1Tech_FSC_St3Sect45_downRight_v0_Prescl, &b_L1Tech_FSC_St3Sect45_downRight_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_FSC_St3Sect45_uppLeft.v0", &L1Tech_FSC_St3Sect45_uppLeft_v0, &b_L1Tech_FSC_St3Sect45_uppLeft_v0);
   fChain->SetBranchAddress("L1Tech_FSC_St3Sect45_uppLeft.v0_Prescl", &L1Tech_FSC_St3Sect45_uppLeft_v0_Prescl, &b_L1Tech_FSC_St3Sect45_uppLeft_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_FSC_St3Sect45_uppRight.v0", &L1Tech_FSC_St3Sect45_uppRight_v0, &b_L1Tech_FSC_St3Sect45_uppRight_v0);
   fChain->SetBranchAddress("L1Tech_FSC_St3Sect45_uppRight.v0_Prescl", &L1Tech_FSC_St3Sect45_uppRight_v0_Prescl, &b_L1Tech_FSC_St3Sect45_uppRight_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_HCAL_HBHE_totalOR.v0", &L1Tech_HCAL_HBHE_totalOR_v0, &b_L1Tech_HCAL_HBHE_totalOR_v0);
   fChain->SetBranchAddress("L1Tech_HCAL_HBHE_totalOR.v0_Prescl", &L1Tech_HCAL_HBHE_totalOR_v0_Prescl, &b_L1Tech_HCAL_HBHE_totalOR_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_HCAL_HF_MMP_or_MPP.v1", &L1Tech_HCAL_HF_MMP_or_MPP_v1, &b_L1Tech_HCAL_HF_MMP_or_MPP_v1);
   fChain->SetBranchAddress("L1Tech_HCAL_HF_MMP_or_MPP.v1_Prescl", &L1Tech_HCAL_HF_MMP_or_MPP_v1_Prescl, &b_L1Tech_HCAL_HF_MMP_or_MPP_v1_Prescl);
   fChain->SetBranchAddress("L1Tech_HCAL_HF_coincidence_PM.v2", &L1Tech_HCAL_HF_coincidence_PM_v2, &b_L1Tech_HCAL_HF_coincidence_PM_v2);
   fChain->SetBranchAddress("L1Tech_HCAL_HF_coincidence_PM.v2_Prescl", &L1Tech_HCAL_HF_coincidence_PM_v2_Prescl, &b_L1Tech_HCAL_HF_coincidence_PM_v2_Prescl);
   fChain->SetBranchAddress("L1Tech_HCAL_HF_single_channel.v0", &L1Tech_HCAL_HF_single_channel_v0, &b_L1Tech_HCAL_HF_single_channel_v0);
   fChain->SetBranchAddress("L1Tech_HCAL_HF_single_channel.v0_Prescl", &L1Tech_HCAL_HF_single_channel_v0_Prescl, &b_L1Tech_HCAL_HF_single_channel_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_HCAL_HO_totalOR.v0", &L1Tech_HCAL_HO_totalOR_v0, &b_L1Tech_HCAL_HO_totalOR_v0);
   fChain->SetBranchAddress("L1Tech_HCAL_HO_totalOR.v0_Prescl", &L1Tech_HCAL_HO_totalOR_v0_Prescl, &b_L1Tech_HCAL_HO_totalOR_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RB0_Cosmics.v0", &L1Tech_RPC_TTU_RB0_Cosmics_v0, &b_L1Tech_RPC_TTU_RB0_Cosmics_v0);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RB0_Cosmics.v0_Prescl", &L1Tech_RPC_TTU_RB0_Cosmics_v0_Prescl, &b_L1Tech_RPC_TTU_RB0_Cosmics_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBminus1_Cosmics.v0", &L1Tech_RPC_TTU_RBminus1_Cosmics_v0, &b_L1Tech_RPC_TTU_RBminus1_Cosmics_v0);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBminus1_Cosmics.v0_Prescl", &L1Tech_RPC_TTU_RBminus1_Cosmics_v0_Prescl, &b_L1Tech_RPC_TTU_RBminus1_Cosmics_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBminus2_Cosmics.v0", &L1Tech_RPC_TTU_RBminus2_Cosmics_v0, &b_L1Tech_RPC_TTU_RBminus2_Cosmics_v0);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBminus2_Cosmics.v0_Prescl", &L1Tech_RPC_TTU_RBminus2_Cosmics_v0_Prescl, &b_L1Tech_RPC_TTU_RBminus2_Cosmics_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBplus1_Cosmics.v0", &L1Tech_RPC_TTU_RBplus1_Cosmics_v0, &b_L1Tech_RPC_TTU_RBplus1_Cosmics_v0);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBplus1_Cosmics.v0_Prescl", &L1Tech_RPC_TTU_RBplus1_Cosmics_v0_Prescl, &b_L1Tech_RPC_TTU_RBplus1_Cosmics_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBplus2_Cosmics.v0", &L1Tech_RPC_TTU_RBplus2_Cosmics_v0, &b_L1Tech_RPC_TTU_RBplus2_Cosmics_v0);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBplus2_Cosmics.v0_Prescl", &L1Tech_RPC_TTU_RBplus2_Cosmics_v0_Prescl, &b_L1Tech_RPC_TTU_RBplus2_Cosmics_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBst1_collisions.v0", &L1Tech_RPC_TTU_RBst1_collisions_v0, &b_L1Tech_RPC_TTU_RBst1_collisions_v0);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBst1_collisions.v0_Prescl", &L1Tech_RPC_TTU_RBst1_collisions_v0_Prescl, &b_L1Tech_RPC_TTU_RBst1_collisions_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_barrel_Cosmics.v0", &L1Tech_RPC_TTU_barrel_Cosmics_v0, &b_L1Tech_RPC_TTU_barrel_Cosmics_v0);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_barrel_Cosmics.v0_Prescl", &L1Tech_RPC_TTU_barrel_Cosmics_v0_Prescl, &b_L1Tech_RPC_TTU_barrel_Cosmics_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_pointing_Cosmics.v0", &L1Tech_RPC_TTU_pointing_Cosmics_v0, &b_L1Tech_RPC_TTU_pointing_Cosmics_v0);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_pointing_Cosmics.v0_Prescl", &L1Tech_RPC_TTU_pointing_Cosmics_v0_Prescl, &b_L1Tech_RPC_TTU_pointing_Cosmics_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_TOTEM_Diffractive.v0", &L1Tech_TOTEM_Diffractive_v0, &b_L1Tech_TOTEM_Diffractive_v0);
   fChain->SetBranchAddress("L1Tech_TOTEM_Diffractive.v0_Prescl", &L1Tech_TOTEM_Diffractive_v0_Prescl, &b_L1Tech_TOTEM_Diffractive_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_TOTEM_LowMultiplicity.v0", &L1Tech_TOTEM_LowMultiplicity_v0, &b_L1Tech_TOTEM_LowMultiplicity_v0);
   fChain->SetBranchAddress("L1Tech_TOTEM_LowMultiplicity.v0_Prescl", &L1Tech_TOTEM_LowMultiplicity_v0_Prescl, &b_L1Tech_TOTEM_LowMultiplicity_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_TOTEM_MinBias.v0", &L1Tech_TOTEM_MinBias_v0, &b_L1Tech_TOTEM_MinBias_v0);
   fChain->SetBranchAddress("L1Tech_TOTEM_MinBias.v0_Prescl", &L1Tech_TOTEM_MinBias_v0_Prescl, &b_L1Tech_TOTEM_MinBias_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_TOTEM_ZeroBias.v0", &L1Tech_TOTEM_ZeroBias_v0, &b_L1Tech_TOTEM_ZeroBias_v0);
   fChain->SetBranchAddress("L1Tech_TOTEM_ZeroBias.v0_Prescl", &L1Tech_TOTEM_ZeroBias_v0_Prescl, &b_L1Tech_TOTEM_ZeroBias_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_ZDC_Scint_loose_vertex.v0", &L1Tech_ZDC_Scint_loose_vertex_v0, &b_L1Tech_ZDC_Scint_loose_vertex_v0);
   fChain->SetBranchAddress("L1Tech_ZDC_Scint_loose_vertex.v0_Prescl", &L1Tech_ZDC_Scint_loose_vertex_v0_Prescl, &b_L1Tech_ZDC_Scint_loose_vertex_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_ZDC_Scint_minus.v0", &L1Tech_ZDC_Scint_minus_v0, &b_L1Tech_ZDC_Scint_minus_v0);
   fChain->SetBranchAddress("L1Tech_ZDC_Scint_minus.v0_Prescl", &L1Tech_ZDC_Scint_minus_v0_Prescl, &b_L1Tech_ZDC_Scint_minus_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_ZDC_Scint_plus.v0", &L1Tech_ZDC_Scint_plus_v0, &b_L1Tech_ZDC_Scint_plus_v0);
   fChain->SetBranchAddress("L1Tech_ZDC_Scint_plus.v0_Prescl", &L1Tech_ZDC_Scint_plus_v0_Prescl, &b_L1Tech_ZDC_Scint_plus_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_ZDC_Scint_tight_vertex.v0", &L1Tech_ZDC_Scint_tight_vertex_v0, &b_L1Tech_ZDC_Scint_tight_vertex_v0);
   fChain->SetBranchAddress("L1Tech_ZDC_Scint_tight_vertex.v0_Prescl", &L1Tech_ZDC_Scint_tight_vertex_v0_Prescl, &b_L1Tech_ZDC_Scint_tight_vertex_v0_Prescl);
   Notify();
}

Bool_t HltTree_pp::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void HltTree_pp::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t HltTree_pp::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
void HltTree_pp::Close()
{
  f->Close();
}
#endif // #ifdef HltTree_pp_cxx
