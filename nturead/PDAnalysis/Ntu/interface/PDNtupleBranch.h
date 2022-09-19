#ifndef PDNtupleBranch_H
#define PDNtupleBranch_H

#include "NtuAnalysis/Common/interface/NtuEventBranch.h"
#include "PDAnalysis/Ntu/interface/PDNtupleData.h"
#include <string>
#include <vector>

class TBranch;

template <class T>
class PDNtupleBranch: public virtual PDNtupleData, public virtual T,
                      public NtuEventBranch<T> {

 public:

  PDNtupleBranch();
  virtual ~PDNtupleBranch();

  void initTree();

 protected:

  // List of branches

  // hlt path names
  void setBranches_hltlist();
  TBranch* b_nHLTPaths;
  TBranch* b_hltCode;
  TBranch* b_hltName;

  // hlt status
  void setBranches_hlts();
  TBranch* b_nHLTStatus;
  TBranch* b_hltPath;
  TBranch* b_hltVersion;
  TBranch* b_hltScale;
  TBranch* b_hltRun;
  TBranch* b_hltAccept;

  // hlt objects
  void setBranches_hlto();
  TBranch* b_nHLTObjects;
  TBranch* b_hltObjType;
  TBranch* b_hltPt;
  TBranch* b_hltEta;
  TBranch* b_hltPhi;
  TBranch* b_hltPx;
  TBranch* b_hltPy;
  TBranch* b_hltPz;
  TBranch* b_hltE;

  // hlt match
  void setBranches_hltm();
  bool use_hltm;
  TBranch* b_nHLTMatches;
  TBranch* b_hltPathIndex;
  TBranch* b_hltObjectIndex;
  TBranch* b_hltFilter;

  // beam spot
  void setBranches_bspot();
  TBranch* b_bwX;         // beam width
  TBranch* b_bwY;
  TBranch* b_bwXError;
  TBranch* b_bwYError;
  TBranch* b_bsX;         // beam spot
  TBranch* b_bsY;
  TBranch* b_bsZ;
  TBranch* b_bsXError;
  TBranch* b_bsYError;
  TBranch* b_bsZError;
//  TBranch* b_bsSxy;
//  TBranch* b_bsSxz;
//  TBranch* b_bsSyz;
  TBranch* b_bsdXdZ;      // beam slope
  TBranch* b_bsdYdZ;
  TBranch* b_bsdXdZError;
  TBranch* b_bsdYdZError;
  TBranch* b_bsSigmaZ;      // sigma Z
  TBranch* b_bsSigmaZError;

  // met
  void setBranches_met();
  TBranch* b_mEt;
  TBranch* b_mEx;
  TBranch* b_mEy;
  TBranch* b_mEsign;

  // muons
  void setBranches_muons();
  TBranch* b_nMuons;
  TBranch* b_muoPt;
  TBranch* b_muoEta;
  TBranch* b_muoPhi;
  TBranch* b_muoPx;
  TBranch* b_muoPy;
  TBranch* b_muoPz;
  TBranch* b_muoE;
  TBranch* b_muoCharge;
  TBranch* b_muoTrk;
  TBranch* b_muoTrg;
  TBranch* b_muoChaIso;
  TBranch* b_muoNeuIso;
  TBranch* b_muoPhoIso;
  TBranch* b_muoSumCHpt;
  TBranch* b_muoSumCPpt;
  TBranch* b_muoSumNHet;
  TBranch* b_muoSumPHet;
  TBranch* b_muoSumPUpt;
  TBranch* b_muoNumMatches;
  TBranch* b_muoType;
  TBranch* b_muoDb;
  TBranch* b_muoChi2LP;
  TBranch* b_muoTrkKink;
  TBranch* b_muoSegmComp;
//  TBranch* b_muoNumValidHits;
//  TBranch* b_muoNormChi2;
//  TBranch* b_muoNumMuHits;
//  TBranch* b_muoNumPixHits;
//  TBranch* b_muoNumTkHits;
  TBranch* b_muoChi2LM;
  TBranch* b_muoTrkRelChi2;
  TBranch* b_muoGlbTrackTailProb;
  TBranch* b_muoGlbKinkFinderLOG;  
  TBranch* b_muoGlbDeltaEtaPhi;
  TBranch* b_muoStaRelChi2;
  TBranch* b_muoTimeAtIpInOut;
  TBranch* b_muoTimeAtIpInOutErr;
//  TBranch* b_muoInnerChi2;
  TBranch* b_muoIValFrac;
//  TBranch* b_muoValPixHits;
//  TBranch* b_muoNTrkVHits;
//  TBranch* b_muoOuterChi2;
//  TBranch* b_muoGNchi2;
//  TBranch* b_muoVMuHits;
//  TBranch* b_muoQprod;
//  TBranch* b_muoLWH;
  TBranch* b_muoVMuonHitComb;

  // electrons
  void setBranches_electrons();
  TBranch* b_nElectrons;
  TBranch* b_elePt;
  TBranch* b_eleEta;
  TBranch* b_elePhi;
  TBranch* b_elePx;
  TBranch* b_elePy;
  TBranch* b_elePz;
  TBranch* b_eleE;
  TBranch* b_eleCharge;
  TBranch* b_eleTrk;
  TBranch* b_eleTrg;
  TBranch* b_eleChaIso;
  TBranch* b_eleNeuIso;
  TBranch* b_elePhoIso;
  TBranch* b_elePCHIso;
  TBranch* b_eleSumCHpt;
  TBranch* b_eleSumCPpt;
  TBranch* b_eleSumNHet;
  TBranch* b_eleSumPHet;
  TBranch* b_eleSumPUpt;
  TBranch* b_eleAbsEta;
  TBranch* b_eleAEff;
  TBranch* b_eleID;
  TBranch* b_eleConsCha;
  TBranch* b_eleEBEEGap;
  TBranch* b_eleDb;
  TBranch* b_eleGsfPt;
  TBranch* b_eleGsfEta;
  TBranch* b_eleGsfPhi;
  TBranch* b_eleGsfPx;
  TBranch* b_eleGsfPy;
  TBranch* b_eleGsfPz;
  TBranch* b_eleGsfPtAtVtx;
  TBranch* b_eleGsfEtaAtVtx;
  TBranch* b_eleGsfPhiAtVtx;
  TBranch* b_eleGsfPxAtVtx;
  TBranch* b_eleGsfPyAtVtx;
  TBranch* b_eleGsfPzAtVtx;
  TBranch* b_eleGsfChi2;
  TBranch* b_eleGsfNdof;
  TBranch* b_eleGsfLostHits;
  TBranch* b_eleGsfValidHits;
  TBranch* b_eleGsfPVtx;
  TBranch* b_eleGsfQuality;
  TBranch* b_eleGsfHitPattern;
  TBranch* b_eleGsfLayPattern;
  TBranch* b_eleGsfPtError;
  TBranch* b_eleGsfDxy;
  TBranch* b_eleGsfDxyPV;
  TBranch* b_eleGsfDz;
  TBranch* b_eleGsfDzPV;
  TBranch* b_eleGsfExy;
//   TBranch* b_eleGsfExyPV;
  TBranch* b_eleGsfEz;
//   TBranch* b_eleGsfEzPV;
  TBranch* b_eleCtfPt;
  TBranch* b_eleCtfEta;
  TBranch* b_eleCtfPhi;
  TBranch* b_eleCtfPx;
  TBranch* b_eleCtfPy;
  TBranch* b_eleCtfPz;
  TBranch* b_eleCtfTrk;
  TBranch* b_eleCtfChi2;
  TBranch* b_eleCtfNdof;
  TBranch* b_eleCtfLostHits;
  TBranch* b_eleCtfValidHits;
  TBranch* b_eleCtfPVtx;
  TBranch* b_eleCtfQuality;
  TBranch* b_eleCtfHitPattern;
  TBranch* b_eleCtfLayPattern;
  TBranch* b_eleCtfPtError;
  TBranch* b_eleCtfDxy;
  TBranch* b_eleCtfDxyPV;
  TBranch* b_eleCtfDz;
  TBranch* b_eleCtfDzPV;
  TBranch* b_eleCtfExy;
//   TBranch* b_eleCtfExyPV;
  TBranch* b_eleCtfEz;
//   TBranch* b_eleCtfEzPV;
  TBranch* b_eleProps;
  TBranch* b_elePartIso;
  TBranch* b_eleEcalIso;
  TBranch* b_eleHcalIso;
  TBranch* b_eleTrackIso;
  TBranch* b_eleClass;
  TBranch* b_eleCtfGsfOv;
  TBranch* b_eleDEtaElClTkAtCalo;
  TBranch* b_eleDEtaSeSlTkAtCalo;
  TBranch* b_eleDEtaSeClTkAtVtx;
  TBranch* b_eleDEtaSuClTkAtVtx;
  TBranch* b_eleDPhiElClTkAtCalo;
  TBranch* b_eleDPhiSeClTkAtCalo;
  TBranch* b_eleDPhiSuClTkAtVtx;
  TBranch* b_eleEcalE;
  TBranch* b_eleEcalEE;
  TBranch* b_eleEElClOvPOut;
  TBranch* b_eleESeClOvP;
  TBranch* b_eleESeClOvPOut;
  TBranch* b_eleESuClOvP;
  TBranch* b_eleFBrem;
  TBranch* b_eleHcalOvEcal;
  TBranch* b_eleF5x5E1x5;
  TBranch* b_eleF5x5E2x5B;
  TBranch* b_eleF5x5E2x5L;
  TBranch* b_eleF5x5E2x5M;
  TBranch* b_eleF5x5E2x5R;
  TBranch* b_eleF5x5E2x5T;
  TBranch* b_eleF5x5E5x5;
  TBranch* b_eleF5x5EB;
  TBranch* b_eleF5x5EL;
  TBranch* b_eleF5x5ER;
  TBranch* b_eleF5x5ET;
  TBranch* b_eleF5x5HcalOvEcal;
  TBranch* b_eleF5x5HcalOvEcalV;
  TBranch* b_eleF5x5R9;
  TBranch* b_eleF5x5SEtaEta;
  TBranch* b_eleF5x5SIetaIeta;
  TBranch* b_eleF5x5SIphiIphi;
  TBranch* b_eleSCEtaWidth;
  TBranch* b_eleSCPhiWidth;
  TBranch* b_eleSCPreshowerE;
  TBranch* b_eleSCRawE;
  TBranch* b_eleSCEta;
  TBranch* b_eleConvDcot;
  TBranch* b_eleConvDist;
  TBranch* b_eleConvFlags;
  TBranch* b_eleConvRadius;
  TBranch* b_eleConvVtxProb;

  // taus
  void setBranches_taus();
  TBranch* b_nTaus;
  TBranch* b_tauPt;
  TBranch* b_tauEta;
  TBranch* b_tauPhi;
  TBranch* b_tauPx;
  TBranch* b_tauPy;
  TBranch* b_tauPz;
  TBranch* b_tauE;
  TBranch* b_tauCharge;
  TBranch* b_tauTrg;

  // jets
  void setBranches_jets();
  TBranch* b_nJets;
  TBranch* b_jetPt;
  TBranch* b_jetEta;
  TBranch* b_jetPhi;
  TBranch* b_jetPx;
  TBranch* b_jetPy;
  TBranch* b_jetPz;
  TBranch* b_jetE;
  TBranch* b_jetCSV;
  TBranch* b_jetTCHE;
  TBranch* b_jetTrg;
  TBranch* b_jetPF;
  TBranch* b_jetNDau;
  TBranch* b_jetNHF;
  TBranch* b_jetNEF;
  TBranch* b_jetCHF;
  TBranch* b_jetCEF;
  TBranch* b_jetNCH;

  // jet tags
  void setBranches_tags();
  TBranch* b_nTags;
  TBranch* b_tagJet;
  TBranch* b_tagType;
  TBranch* b_tagProb;

  // user info
  void setBranches_userinfo();
  TBranch* b_nUserInfo;
  TBranch* b_useObjType;
  TBranch* b_useObjIndex;
  TBranch* b_useInfoType;
  TBranch* b_useInfoValue;

  // particle flow
  void setBranches_pflow();
  TBranch* b_nPF;
  TBranch* b_pfcPt;
  TBranch* b_pfcEta;
  TBranch* b_pfcPhi;
  TBranch* b_pfcPx;
  TBranch* b_pfcPy;
  TBranch* b_pfcPz;
  TBranch* b_pfcE;
  TBranch* b_pfcCharge;
  TBranch* b_pfcJet;
  TBranch* b_pfcTrk;
  TBranch* b_pfcHasTrackDetails;

  // tracks
  void setBranches_tracks();
  TBranch* b_nTracks;
  TBranch* b_trkPt;
  TBranch* b_trkEta;
  TBranch* b_trkPhi;
  TBranch* b_trkPx;
  TBranch* b_trkPy;
  TBranch* b_trkPz;
  TBranch* b_trkCharge;
  TBranch* b_trkType;
  TBranch* b_trkNext;
  TBranch* b_trkPFC;
  TBranch* b_trkJet;
  TBranch* b_trkPVtx;
//  TBranch* b_trkSVtx;
  TBranch* b_trkQuality;
  TBranch* b_trkHitPattern;
  TBranch* b_trkLayPattern;
  TBranch* b_trkMissingInnerHits;
  TBranch* b_trkNormChi2;
  TBranch* b_trkPtError;
  TBranch* b_trkDxy;
  TBranch* b_trkDz;
  TBranch* b_trkExy;
  TBranch* b_trkEz;
  TBranch* b_trkVtxPx;
  TBranch* b_trkVtxPy;
  TBranch* b_trkVtxPz;
  TBranch* b_trkQOverPError;
  TBranch* b_trkLambdaError;
  TBranch* b_trkPhiError;
  TBranch* b_trkDszError;
  TBranch* b_trkEtaError;
  TBranch* b_trkDxyDszCov;
  TBranch* b_trkLambdaDszCov;
  TBranch* b_trkPhiDxyCov;

  // track perigee parameters
  void setBranches_trkper();
  TBranch* b_nTrkPer;
  TBranch* b_tppTrk;
  TBranch* b_tppQop;
  TBranch* b_tppLam;
  TBranch* b_tppPhi;
  TBranch* b_tppDxy;
  TBranch* b_tppDsz;
  TBranch* b_tppSQopQop;
  TBranch* b_tppSQopLam;
  TBranch* b_tppSQopPhi;
  TBranch* b_tppSQopDxy;
  TBranch* b_tppSQopDsz;
  TBranch* b_tppSLamLam;
  TBranch* b_tppSLamPhi;
  TBranch* b_tppSLamDxy;
  TBranch* b_tppSLamDsz;
  TBranch* b_tppSPhiPhi;
  TBranch* b_tppSPhiDxy;
  TBranch* b_tppSPhiDsz;
  TBranch* b_tppSDxyDxy;
  TBranch* b_tppSDxyDsz;
  TBranch* b_tppSDszDsz;

  // primary vertices
  void setBranches_pVertices();
  TBranch* b_nPVTotal;
  TBranch* b_nPVertices;
  TBranch* b_pvtX;
  TBranch* b_pvtY;
  TBranch* b_pvtZ;
  TBranch* b_pvtSxx;
  TBranch* b_pvtSyy;
  TBranch* b_pvtSzz;
  TBranch* b_pvtSxy;
  TBranch* b_pvtSxz;
  TBranch* b_pvtSyz;
//  TBranch* b_pvtCovariance;
  TBranch* b_pvtNTracks;
//  TBranch* b_pvtNormChi2;
  TBranch* b_pvtChi2;
  TBranch* b_pvtNDOF;
  TBranch* b_pvtBadQuality;

  // secondary vertices
  void setBranches_sVertices();
  TBranch* b_nSVertices;
  TBranch* b_svtX;
  TBranch* b_svtY;
  TBranch* b_svtZ;
  TBranch* b_svtSxx;
  TBranch* b_svtSyy;
  TBranch* b_svtSzz;
  TBranch* b_svtSxy;
  TBranch* b_svtSxz;
  TBranch* b_svtSyz;
  TBranch* b_svtPVtx;
  TBranch* b_svtDirX;
  TBranch* b_svtDirY;
  TBranch* b_svtDirZ;
  TBranch* b_svtType;
  TBranch* b_svtNTracks;
  TBranch* b_svtChi2;
  TBranch* b_svtNDOF;
  TBranch* b_svtMass;
  TBranch* b_svtDist2D;
  TBranch* b_svtSigma2D;
  TBranch* b_svtDist3D;
  TBranch* b_svtSigma3D;
  TBranch* b_svtJet;
  TBranch* b_svtBadQuality;

  // composite particle vertices
  void setBranches_vsub();
  TBranch* b_nCompVts;
  TBranch* b_subPart;
  TBranch* b_subSVtx;
  TBranch* b_subMass;

  // impact parameters
  void setBranches_tkips();
  TBranch* b_nTkIPs;
  TBranch* b_tipTrk;
  TBranch* b_tipSVtx;
  TBranch* b_tipMass;
  TBranch* b_tipDxy;
  TBranch* b_tipDz;
  TBranch* b_tipExy;
  TBranch* b_tipEz;

  // momenta at vertices
  void setBranches_vtxps();
  TBranch* b_nVtxPs;
  TBranch* b_tvpTip;
  TBranch* b_tvpPt;
  TBranch* b_tvpEta;
  TBranch* b_tvpPhi;
  TBranch* b_tvpPx;
  TBranch* b_tvpPy;
  TBranch* b_tvpPz;

  // PU weight
  void setBranches_puwgt();
  TBranch* b_puWeight;

  // gen particles
  void setBranches_gen();
  TBranch* b_nGenP;
  TBranch* b_genId;
  TBranch* b_genStatus;
  TBranch* b_genMother;
  TBranch* b_genPartner;
  TBranch* b_genPt;
  TBranch* b_genEta;
  TBranch* b_genPhi;
  TBranch* b_genPx;
  TBranch* b_genPy;
  TBranch* b_genPz;
  TBranch* b_genE;
  TBranch* b_genCharge;
  TBranch* b_genMass;
//  TBranch* b_genJet;
  TBranch* b_genVx;
  TBranch* b_genVy;
  TBranch* b_genVz;

  // gen jets
  void setBranches_gpj();
  TBranch* b_nGenJets;
  TBranch* b_gpjPt;
  TBranch* b_gpjEta;
  TBranch* b_gpjPhi;
  TBranch* b_gpjPx;
  TBranch* b_gpjPy;
  TBranch* b_gpjPz;
  TBranch* b_gpjE;
  TBranch* b_gpjNDau;
  TBranch* b_gpjReco;

 private:

  // dummy copy constructor and assignment
  PDNtupleBranch           ( const PDNtupleBranch& td );
  PDNtupleBranch& operator=( const PDNtupleBranch& td );

  const char* bCat( const std::string& name, const std::string& type );

};

#include "PDNtupleBranch.hpp"

#endif

