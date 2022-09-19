#ifndef PDNtupleData_H
#define PDNtupleData_H

#include <string>
#include <vector>

#include "NtuAnalysis/Common/interface/NtuData.h"

class PDNtupleData: public NtuData<float> {

 public:

  PDNtupleData();
  virtual ~PDNtupleData();

 protected:

  // Declaration of leaf types

  // hlt path names
  bool use_hltlist;
  int nHLTPaths;
  std::vector<int        >* hltCode;
  std::vector<std::string>* hltName;

  // hlt status
  bool use_hlts;
  int nHLTStatus;
  std::vector<int >* hltPath;
  std::vector<int >* hltVersion;
  std::vector<int >* hltScale;
  std::vector<bool>* hltRun;
  std::vector<bool>* hltAccept;

  // hlt objects
  bool use_hlto;
  bool use_hlto_sphe;
  bool use_hlto_cart;
  int nHLTObjects;
  std::vector<int   >* hltObjType;
  std::vector<number>* hltPt;
  std::vector<number>* hltEta;
  std::vector<number>* hltPhi;
  std::vector<number>* hltPx;
  std::vector<number>* hltPy;
  std::vector<number>* hltPz;
  std::vector<number>* hltE;

  // hlt match
  bool use_hltm;
  int nHLTMatches;
  std::vector<int>*  hltPathIndex;
  std::vector<int>*  hltObjectIndex;
  std::vector<char>* hltFilter;

  // beam spot
  bool use_bspot;
  number bwX;         // beam width
  number bwY;
  number bwXError;
  number bwYError;
  number bsX;         // beam spot
  number bsY;
  number bsZ;
  number bsXError;
  number bsYError;
  number bsZError;
//  number bsSxy;
//  number bsSxz;
//  number bsSyz;
  number bsdXdZ;      // beam slope
  number bsdYdZ;
  number bsdXdZError;
  number bsdYdZError;
  number bsSigmaZ;    // sigma Z
  number bsSigmaZError;

  // met
  bool use_met;
  number mEt;
  number mEx;
  number mEy;
  number mEsign;

  // muons
  bool use_muons;
  bool use_muons_sphe;
  bool use_muons_cart;
  int nMuons;
  std::vector<number>* muoPt;
  std::vector<number>* muoEta;
  std::vector<number>* muoPhi;
  std::vector<number>* muoPx;
  std::vector<number>* muoPy;
  std::vector<number>* muoPz;
  std::vector<number>* muoE;
  std::vector<int   >* muoCharge;
  std::vector<int   >* muoTrk;        // track index
  std::vector<int   >* muoTrg;        // trigger object index
  std::vector<number>* muoChaIso;     //.chargedHadronIso()
  std::vector<number>* muoNeuIso;     //.neutralHadronIso()
  std::vector<number>* muoPhoIso;     //.photonIso()
  std::vector<number>* muoSumCHpt;    //.pfIsolationR04().sumChargedHadronPt
  std::vector<number>* muoSumCPpt;    //.pfIsolationR04().sumChargedParticlePt
  std::vector<number>* muoSumNHet;    //.pfIsolationR04().sumNeutralHadronEt
  std::vector<number>* muoSumPHet;    //.pfIsolationR04().sumPhotonEt
  std::vector<number>* muoSumPUpt;    //.pfIsolationR04().sumPUPt
  std::vector<int   >* muoNumMatches; //.numberOfMatchedStations()
  std::vector<int   >* muoType;       // bitted, see PDEnumString::muonType
  std::vector<number>* muoDb;         //.dB()
  std::vector<number>* muoChi2LP;     //.combinedQuality().chi2LocalPosition
  std::vector<number>* muoTrkKink;    //.combinedQuality().trkKink
  std::vector<number>* muoSegmComp;   //muon::segmentCompatibility
//  std::vector<int   >* muoNumValidHits; //. innerTrack()->numberOfValidHits()
//  std::vector<number>* muoNormChi2;     //.globalTrack()->normalizedChi2()
//  std::vector<int   >* muoNumMuHits;    //.globalTrack()->hitPattern().numberOfValidMuonHits()
//  std::vector<int   >* muoNumPixHits;   //.globalTrack()->hitPattern().numberOfValidPixelHits()
//  std::vector<int   >* muoNumTkHits;    //.globalTrack()->hitPattern().numberOfValidTrackerHits()
  std::vector<number>* muoChi2LM;
  std::vector<number>* muoTrkRelChi2;
  std::vector<number>* muoGlbTrackTailProb;
  std::vector<number>* muoGlbKinkFinderLOG;  
  std::vector<number>* muoGlbDeltaEtaPhi;
  std::vector<number>* muoStaRelChi2;
  std::vector<number>* muoTimeAtIpInOut;
  std::vector<number>* muoTimeAtIpInOutErr;
//  std::vector<number>* muoInnerChi2;
  std::vector<number>* muoIValFrac;
//  std::vector<number>* muoValPixHits;
//  std::vector<number>* muoNTrkVHits;
//  std::vector<number>* muoOuterChi2;
//  std::vector<number>* muoGNchi2;
//  std::vector<number>* muoVMuHits;
//  std::vector<number>* muoQprod;
//  std::vector<number>* muoLWH;
  std::vector<number>* muoVMuonHitComb;

  // electrons
  bool use_electrons;
  bool use_electrons_sphe;
  bool use_electrons_cart;
  int nElectrons;
  std::vector<number>* elePt;
  std::vector<number>* eleEta;
  std::vector<number>* elePhi;
  std::vector<number>* elePx;
  std::vector<number>* elePy;
  std::vector<number>* elePz;
  std::vector<number>* eleE;
  std::vector<int   >* eleCharge;
  std::vector<int   >* eleTrk;     // track index
  std::vector<int   >* eleTrg;     // trigger object index
  std::vector<number>* eleChaIso;  //.chargedHadronIso()
  std::vector<number>* eleNeuIso;  //.neutralHadronIso()
  std::vector<number>* elePhoIso;  //.photonIso()
  std::vector<number>* elePCHIso;  //.puChargedHadronIso()
  std::vector<number>* eleSumCHpt; //.pfIsolationVariables().sumChargedHadronPt
  std::vector<number>* eleSumCPpt; //.pfIsolationVariables().sumChargedParticlePt
  std::vector<number>* eleSumNHet; //.pfIsolationVariables().sumNeutralHadronEt
  std::vector<number>* eleSumPHet; //.pfIsolationVariables().sumPhotonEt
  std::vector<number>* eleSumPUpt; //.pfIsolationVariables().sumPUPt
  std::vector<number>* eleAbsEta;  //fabs(.superCluster()->eta())
  std::vector<number>* eleAEff;
    //ElectronEffectiveArea::GetElectronEffectiveArea(
    //ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, absEta,
    //ElectronEffectiveArea::kEleEAData2011 )
  std::vector<int   >* eleID;      //.electronID( "eidMediumMC" );
  std::vector<bool  >* eleConsCha; //.isGsfCtfScPixChargeConsistent();
  std::vector<bool  >* eleEBEEGap; //.isEBEEGap();
  std::vector<number>* eleDb;      //.dB();
  std::vector<number>* eleGsfPt;
  std::vector<number>* eleGsfEta;
  std::vector<number>* eleGsfPhi;
  std::vector<number>* eleGsfPx;
  std::vector<number>* eleGsfPy;
  std::vector<number>* eleGsfPz;
  std::vector<number>* eleGsfPtAtVtx;
  std::vector<number>* eleGsfEtaAtVtx;
  std::vector<number>* eleGsfPhiAtVtx;
  std::vector<number>* eleGsfPxAtVtx;
  std::vector<number>* eleGsfPyAtVtx;
  std::vector<number>* eleGsfPzAtVtx;
  std::vector<number>* eleGsfChi2;
  std::vector<number>* eleGsfNdof;
  std::vector<int   >* eleGsfLostHits;
  std::vector<int   >* eleGsfValidHits;
  std::vector<int   >* eleGsfPVtx;
  std::vector<int   >* eleGsfQuality;
  std::vector<int   >* eleGsfHitPattern;
  std::vector<int   >* eleGsfLayPattern;
  std::vector<number>* eleGsfPtError;
  std::vector<number>* eleGsfDxy;
  std::vector<number>* eleGsfDxyPV;
  std::vector<number>* eleGsfDz;
  std::vector<number>* eleGsfDzPV;
  std::vector<number>* eleGsfExy;
//   std::vector<number>* eleGsfExyPV;
  std::vector<number>* eleGsfEz;
//   std::vector<number>* eleGsfEzPV;
  std::vector<number>* eleCtfPt;
  std::vector<number>* eleCtfEta;
  std::vector<number>* eleCtfPhi;
  std::vector<number>* eleCtfPx;
  std::vector<number>* eleCtfPy;
  std::vector<number>* eleCtfPz;
  std::vector<int   >* eleCtfTrk;
  std::vector<number>* eleCtfChi2;
  std::vector<number>* eleCtfNdof;
  std::vector<int   >* eleCtfLostHits;
  std::vector<int   >* eleCtfValidHits;
  std::vector<int   >* eleCtfPVtx;
  std::vector<int   >* eleCtfQuality;
  std::vector<int   >* eleCtfHitPattern;
  std::vector<int   >* eleCtfLayPattern;
  std::vector<number>* eleCtfPtError;
  std::vector<number>* eleCtfDxy;
  std::vector<number>* eleCtfDxyPV;
  std::vector<number>* eleCtfDz;
  std::vector<number>* eleCtfDzPV;
  std::vector<number>* eleCtfExy;
//   std::vector<number>* eleCtfExyPV;
  std::vector<number>* eleCtfEz;
//   std::vector<number>* eleCtfEzPV;

  std::vector<int   >* eleProps;   // bitmask containing the above bools
  std::vector<number>* elePartIso;  //.particleIso()
  std::vector<number>* eleEcalIso; //.ecalIso()
  std::vector<number>* eleHcalIso; //.hcalIso()
  std::vector<number>* eleTrackIso;  //.trackIso()
  std::vector<int   >* eleClass;   //.classification()
  std::vector<number>* eleCtfGsfOv; //.ctfGsfOverlap()
  std::vector<number>* eleDEtaElClTkAtCalo; //.deltaEtaEleClusterTrackAtCalo()
  std::vector<number>* eleDEtaSeSlTkAtCalo; //.deltaEtaSeedClusterTrackAtCalo()
  std::vector<number>* eleDEtaSeClTkAtVtx;  //.deltaEtaSeedClusterTrackAtVtx()
  std::vector<number>* eleDEtaSuClTkAtVtx;  //.deltaEtaSuperClusterTrackAtVtx()
  std::vector<number>* eleDPhiElClTkAtCalo; //.deltaPhiEleClusterTrackAtCalo()
  std::vector<number>* eleDPhiSeClTkAtCalo; //.deltaPhiSeedClusterTrackAtCalo()
  std::vector<number>* eleDPhiSuClTkAtVtx;  //.deltaPhiSuperClusterTrackAtVtx()
  std::vector<number>* eleEcalE;   //.ecalEnergy()
  std::vector<number>* eleEcalEE;  //.ecalEnergyError()
  std::vector<number>* eleEElClOvPOut;  //.eEleClusterOverPout()
  std::vector<number>* eleESeClOvP;     //.eSeedClusterOverP()
  std::vector<number>* eleESeClOvPOut;  //.eSeedClusterOverPout()
  std::vector<number>* eleESuClOvP;     //.eSuperClusterOverP()
  std::vector<number>* eleFBrem;   //.fbrem()
  std::vector<number>* eleHcalOvEcal;   //.hcalOverEcal()
  std::vector<number>* eleF5x5E1x5;  //.full5x5_e1x5()
  std::vector<number>* eleF5x5E2x5B; //.full5x5_e2x5Bottom()
  std::vector<number>* eleF5x5E2x5L; //.full5x5_e2x5Left()
  std::vector<number>* eleF5x5E2x5M; //.full5x5_e2x5Max()
  std::vector<number>* eleF5x5E2x5R; //.full5x5_e2x5Right()
  std::vector<number>* eleF5x5E2x5T; //.full5x5_e2x5Top()
  std::vector<number>* eleF5x5E5x5;  //.full5x5_e5x5()
  std::vector<number>* eleF5x5EB;   //.full5x5_eBottom()
  std::vector<number>* eleF5x5EL;   //.full5x5_eLeft()
  std::vector<number>* eleF5x5ER;   //.full5x5_eRight()
  std::vector<number>* eleF5x5ET;   //.full5x5_eTop()
  std::vector<number>* eleF5x5HcalOvEcal;  //.full5x5_hcalOverEcal()
  std::vector<number>* eleF5x5HcalOvEcalV; //.full5x5_hcalOverEcalValid() // MG: it is actually a bool returned as a float...
  std::vector<number>* eleF5x5R9;   //.full5x5_r9()
  std::vector<number>* eleF5x5SEtaEta; //.full5x5_sigmaEtaEta()
  std::vector<number>* eleF5x5SIetaIeta; //.full5x5_sigmaIetaIeta()
  std::vector<number>* eleF5x5SIphiIphi; //.full5x5_sigmaIphiIphi()
  std::vector<number>* eleSCEtaWidth; //.superCluster()->etaWidth()
  std::vector<number>* eleSCPhiWidth; //.superCluster()->phiWidth()
  std::vector<number>* eleSCPreshowerE; //.superCluster()->preshowerEnergy()
  std::vector<number>* eleSCRawE; //.superCluster()->rawEnergy()
  std::vector<number>* eleSCEta; //.superCluster()->eta()
  std::vector<number>* eleConvDcot; //.convDcot()
  std::vector<number>* eleConvDist; //.convDist()
  std::vector<int   >* eleConvFlags; //.convFlags()
  std::vector<number>* eleConvRadius; //.convRadius()
  std::vector<number>* eleConvVtxProb; // vertex fit probability of electron conversion (code to fill it taken from MVAHZZ)

  // taus
  bool use_taus;
  bool use_taus_sphe;
  bool use_taus_cart;
  int nTaus;
  std::vector<number>* tauPt;
  std::vector<number>* tauEta;
  std::vector<number>* tauPhi;
  std::vector<number>* tauPx;
  std::vector<number>* tauPy;
  std::vector<number>* tauPz;
  std::vector<number>* tauE;
  std::vector<int   >* tauCharge;
  std::vector<int   >* tauTrg;     // trigger object index

  // jets
  bool use_jets;
  bool use_jets_sphe;
  bool use_jets_cart;
  int nJets;
  std::vector<number>* jetPt;
  std::vector<number>* jetEta;
  std::vector<number>* jetPhi;
  std::vector<number>* jetPx;
  std::vector<number>* jetPy;
  std::vector<number>* jetPz;
  std::vector<number>* jetE;
  std::vector<number>* jetCSV;  //.bDiscriminator( "combinedSecondaryVertexBJetTags" )
  std::vector<number>* jetTCHE; //.bDiscriminator( "trackCountingHighEffBJetTags" )
  std::vector<int   >* jetTrg;  // trigger object index
  std::vector<bool  >* jetPF;   //.isPFJet()
  std::vector<int   >* jetNDau; //.numberOfDaughters()
  std::vector<number>* jetNHF;  //.neutralHadronEnergyFraction()
  std::vector<number>* jetNEF;  //.neutralEmEnergyFraction()
  std::vector<number>* jetCHF;  //.chargedHadronEnergyFraction()
  std::vector<number>* jetCEF;  //.chargedEmEnergyFraction()
  std::vector<number>* jetNCH;  //.chargedMultiplicity()

  // jet tags
  bool use_tags;
  int nTags;
  std::vector<int   >* tagJet;   // jet index
  std::vector<int   >* tagType;  // jet tag type ( see PDEnumString::csvType )
  std::vector<number>* tagProb;  //.bDiscriminator( label )

  // user info
  bool use_info;
  int nUserInfo;
  std::vector<int   >* useObjType;
  std::vector<int   >* useObjIndex;
  std::vector<int   >* useInfoType;
  std::vector<number>* useInfoValue;

  // particle flow
  bool use_pflow;
  bool use_pflow_sphe;
  bool use_pflow_cart;
  int nPF;
  std::vector<number>* pfcPt;
  std::vector<number>* pfcEta;
  std::vector<number>* pfcPhi;
  std::vector<number>* pfcPx;
  std::vector<number>* pfcPy;
  std::vector<number>* pfcPz;
  std::vector<number>* pfcE;
  std::vector<int   >* pfcCharge;
  std::vector<int   >* pfcJet;      //   jet index
  std::vector<int   >* pfcTrk;      // track index
  std::vector<int   >* pfcHasTrackDetails;

  // tracks
  bool use_tracks;
  bool use_tracks_sphe;
  bool use_tracks_cart;
  int nTracks;
  std::vector<number>* trkPt;
  std::vector<number>* trkEta;
  std::vector<number>* trkPhi;
  std::vector<number>* trkPx;
  std::vector<number>* trkPy;
  std::vector<number>* trkPz;
  std::vector<int   >* trkCharge;
  std::vector<int   >* trkType;         // bitted, see PDEnumString::trackType
  std::vector<int   >* trkNext;         // next track linked to the same cand
  std::vector<int   >* trkPFC;          // PFCandidate index
  std::vector<int   >* trkJet;          //         jet index
  std::vector<int   >* trkPVtx;         //   primary vertex index
//  std::vector<int   >* trkSVtx;         // secondary vertex index
  std::vector<int   >* trkQuality;      // qualityMask()
  std::vector<int   >* trkHitPattern;   // MMPPTTVV (Muons,Pixel,Tracker,Valid)
  std::vector<int   >* trkMissingInnerHits;
  std::vector<int   >* trkLayPattern;   // 
  std::vector<number>* trkNormChi2;     // normalizedChi2()
  std::vector<number>* trkPtError;      // pt error
  std::vector<number>* trkDxy;          // dxy()
  std::vector<number>* trkDz;           // dz()
  std::vector<number>* trkExy;          // dxy()
  std::vector<number>* trkEz;           // dz()
  std::vector<number>* trkVtxPx;        // momentum at nearest point to BS
  std::vector<number>* trkVtxPy;
  std::vector<number>* trkVtxPz;
  std::vector<number>* trkQOverPError;
  std::vector<number>* trkLambdaError;
  std::vector<number>* trkPhiError;
  std::vector<number>* trkDszError;
  std::vector<number>* trkEtaError;
  std::vector<number>* trkDxyDszCov;
  std::vector<number>* trkLambdaDszCov;
  std::vector<number>* trkPhiDxyCov;

  // track perigee parameters
  bool use_trkper;
  int nTrkPer;
  std::vector<int   >* tppTrk;
  std::vector<number>* tppQop;
  std::vector<number>* tppLam;
  std::vector<number>* tppPhi;
  std::vector<number>* tppDxy;
  std::vector<number>* tppDsz;
  std::vector<number>* tppSQopQop;
  std::vector<number>* tppSQopLam;
  std::vector<number>* tppSQopPhi;
  std::vector<number>* tppSQopDxy;
  std::vector<number>* tppSQopDsz;
  std::vector<number>* tppSLamLam;
  std::vector<number>* tppSLamPhi;
  std::vector<number>* tppSLamDxy;
  std::vector<number>* tppSLamDsz;
  std::vector<number>* tppSPhiPhi;
  std::vector<number>* tppSPhiDxy;
  std::vector<number>* tppSPhiDsz;
  std::vector<number>* tppSDxyDxy;
  std::vector<number>* tppSDxyDsz;
  std::vector<number>* tppSDszDsz;

  // primary vertices
  bool use_pvts;
  int nPVTotal;
  int nPVertices;
  std::vector<number>* pvtX;
  std::vector<number>* pvtY;
  std::vector<number>* pvtZ;
  std::vector<number>* pvtSxx;
  std::vector<number>* pvtSyy;
  std::vector<number>* pvtSzz;
  std::vector<number>* pvtSxy;
  std::vector<number>* pvtSxz;
  std::vector<number>* pvtSyz;
//  std::vector<
//  std::vector<number>
//                    >* pvtCovariance;
  std::vector<int   >* pvtNTracks;
//  std::vector<number>* pvtNormChi2;
  std::vector<number>* pvtChi2;
  std::vector<number>* pvtNDOF;
  std::vector<int   >* pvtBadQuality; // 0=good, 1=fake, 2=invalid

  // secondary vertices
  bool use_svts;
  int nSVertices;
  std::vector<number>* svtX;
  std::vector<number>* svtY;
  std::vector<number>* svtZ;
  std::vector<number>* svtSxx;
  std::vector<number>* svtSyy;
  std::vector<number>* svtSzz;
  std::vector<number>* svtSxy;
  std::vector<number>* svtSxz;
  std::vector<number>* svtSyz;
  std::vector<int   >* svtPVtx;
  std::vector<number>* svtDirX;
  std::vector<number>* svtDirY;
  std::vector<number>* svtDirZ;
  std::vector<int   >* svtType;
  std::vector<int   >* svtNTracks;
  std::vector<number>* svtChi2;
  std::vector<number>* svtNDOF;
  std::vector<number>* svtMass;
  std::vector<number>* svtDist2D;
  std::vector<number>* svtSigma2D;
  std::vector<number>* svtDist3D;
  std::vector<number>* svtSigma3D;
  std::vector<int   >* svtJet;
  std::vector<int   >* svtBadQuality; // 0=good, 1=fake, 2=invalid

  // composite particle vertices
  bool use_vsub;
  int nCompVts;
  std::vector<int   >* subPart;     // daughter vertex index
  std::vector<int   >* subSVtx;     // global   vertex index
  std::vector<number>* subMass;

  // impact parameters
  bool use_tkips;
  int nTkIPs;
  std::vector<int   >* tipTrk;      // track index
  std::vector<int   >* tipSVtx;     // secondary vertex index
  std::vector<number>* tipMass;     // mass hypothesis
  std::vector<number>* tipDxy;      // dxy()
  std::vector<number>* tipDz;       // dz()
  std::vector<number>* tipExy;      // dxy()
  std::vector<number>* tipEz;       // dz()

  // momenta at vertices
  bool use_vtxps;
  bool use_vtxps_sphe;
  bool use_vtxps_cart;
  int nVtxPs;
  std::vector<int   >* tvpTip;      // impact parameter index
  std::vector<number>* tvpPt;
  std::vector<number>* tvpEta;
  std::vector<number>* tvpPhi;
  std::vector<number>* tvpPx;
  std::vector<number>* tvpPy;
  std::vector<number>* tvpPz;

  // PU weight
  bool use_puwgt;
  number puWeight;

  // gen particles
  bool use_gen;
  bool use_gen_sphe;
  bool use_gen_cart;
  int nGenP;
  std::vector<int   >* genId;
  std::vector<int   >* genStatus;
  std::vector<int   >* genMother;  // mother  index
  std::vector<int   >* genPartner; // partner index
  std::vector<number>* genPt;
  std::vector<number>* genEta;
  std::vector<number>* genPhi;
  std::vector<number>* genPx;
  std::vector<number>* genPy;
  std::vector<number>* genPz;
  std::vector<number>* genE;
  std::vector<int   >* genCharge;
  std::vector<number>* genMass;
//  std::vector<int   >* genJet;
  std::vector<number>* genVx;
  std::vector<number>* genVy;
  std::vector<number>* genVz;

  // gen jets
  bool use_gpj;
  bool use_gpj_sphe;
  bool use_gpj_cart;
  int nGenJets;
  std::vector<number>* gpjPt;
  std::vector<number>* gpjEta;
  std::vector<number>* gpjPhi;
  std::vector<number>* gpjPx;
  std::vector<number>* gpjPy;
  std::vector<number>* gpjPz;
  std::vector<number>* gpjE;
  std::vector<int   >* gpjNDau; //.numberOfDaughters()
  std::vector<int   >* gpjReco; // reco jet index

 private:

  // dummy copy constructor and assignment
  PDNtupleData           ( const PDNtupleData& td );
  PDNtupleData& operator=( const PDNtupleData& td );

};

#endif

