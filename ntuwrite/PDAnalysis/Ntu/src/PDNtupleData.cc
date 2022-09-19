#include <iostream>

#include "PDAnalysis/Ntu/interface/PDNtupleData.h"

PDNtupleData::PDNtupleData() {

  // hlt path names
  use_hltlist = false;
  hltCode = new std::vector<int        >;
  hltName = new std::vector<std::string>;

  // hlt status
  use_hlts = false;
  hltPath         = new std::vector<int >;
  hltVersion      = new std::vector<int >;
  hltScale        = new std::vector<int >;
  hltRun          = new std::vector<bool>;
  hltAccept       = new std::vector<bool>;

  // hlt objects
  use_hlto      = false;
  use_hlto_sphe = true;
  use_hlto_cart = false;
  hltObjType      = new std::vector<int   >;
  hltPt           = new std::vector<number>;
  hltEta          = new std::vector<number>;
  hltPhi          = new std::vector<number>;
  hltPx           = new std::vector<number>;
  hltPy           = new std::vector<number>;
  hltPz           = new std::vector<number>;
  hltE            = new std::vector<number>;

  // hlt match
  use_hltm = false;
  hltPathIndex   = new std::vector<int>;
  hltObjectIndex = new std::vector<int>;
  hltFilter      = new std::vector<char>;

  // beam spot
  use_bspot = false;

  // met
  use_met = false;

  // muons
  use_muons      = false;
  use_muons_sphe = true;
  use_muons_cart = false;
  nMuons = 0;
  muoPt           = new std::vector<number>;
  muoEta          = new std::vector<number>;
  muoPhi          = new std::vector<number>;
  muoPx           = new std::vector<number>;
  muoPy           = new std::vector<number>;
  muoPz           = new std::vector<number>;
  muoE            = new std::vector<number>;
  muoCharge       = new std::vector<int   >;
  muoTrk          = new std::vector<int   >;
  muoTrg          = new std::vector<int   >;
  muoChaIso       = new std::vector<number>;
  muoNeuIso       = new std::vector<number>;
  muoPhoIso       = new std::vector<number>;
  muoSumCHpt      = new std::vector<number>;
  muoSumCPpt      = new std::vector<number>;
  muoSumNHet      = new std::vector<number>;
  muoSumPHet      = new std::vector<number>;
  muoSumPUpt      = new std::vector<number>;
  muoNumMatches   = new std::vector<int   >;
  muoType         = new std::vector<int   >;
  muoDb           = new std::vector<number>;
  muoChi2LP       = new std::vector<number>;
  muoTrkKink      = new std::vector<number>;
  muoSegmComp     = new std::vector<number>;
//  muoNumValidHits = new std::vector<int   >;
//  muoNormChi2     = new std::vector<number>;
//  muoNumMuHits    = new std::vector<int   >;
//  muoNumPixHits   = new std::vector<int   >;
//  muoNumTkHits    = new std::vector<int   >;
  muoChi2LM           = new std::vector<number>;
  muoTrkRelChi2       = new std::vector<number>;
  muoGlbTrackTailProb = new std::vector<number>;
  muoGlbKinkFinderLOG = new std::vector<number>; 
  muoGlbDeltaEtaPhi   = new std::vector<number>;
  muoStaRelChi2       = new std::vector<number>;
  muoTimeAtIpInOut    = new std::vector<number>;
  muoTimeAtIpInOutErr = new std::vector<number>;
//  muoInnerChi2        = new std::vector<number>;
  muoIValFrac         = new std::vector<number>;
//  muoValPixHits       = new std::vector<number>;
//  muoNTrkVHits        = new std::vector<number>;
//  muoOuterChi2        = new std::vector<number>;
//  muoGNchi2           = new std::vector<number>;
//  muoVMuHits          = new std::vector<number>;
//  muoQprod            = new std::vector<number>;
//  muoLWH              = new std::vector<number>;
  muoVMuonHitComb     = new std::vector<number>; 

  // electrons
  use_electrons      = false;
  use_electrons_sphe = true;
  use_electrons_cart = false;
  nElectrons = 0;
  elePt           = new std::vector<number>;
  eleEta          = new std::vector<number>;
  elePhi          = new std::vector<number>;
  elePx           = new std::vector<number>;
  elePy           = new std::vector<number>;
  elePz           = new std::vector<number>;
  eleE            = new std::vector<number>;
  eleCharge       = new std::vector<int   >;
  eleTrk          = new std::vector<int   >;
  eleTrg          = new std::vector<int   >;
  eleChaIso       = new std::vector<number>;
  eleNeuIso       = new std::vector<number>;
  elePhoIso       = new std::vector<number>;
  elePCHIso       = new std::vector<number>;
  eleSumCHpt          = new std::vector<number>;
  eleSumCPpt          = new std::vector<number>;
  eleSumNHet          = new std::vector<number>;
  eleSumPHet          = new std::vector<number>;
  eleSumPUpt          = new std::vector<number>;
  eleAbsEta       = new std::vector<number>;
  eleAEff         = new std::vector<number>;
  eleID           = new std::vector<int   >;
  eleConsCha      = new std::vector<bool  >;
  eleEBEEGap      = new std::vector<bool  >;
  eleDb           = new std::vector<number>;
  eleGsfPt            = new std::vector<number>;
  eleGsfEta           = new std::vector<number>;
  eleGsfPhi           = new std::vector<number>;
  eleGsfPx            = new std::vector<number>;
  eleGsfPy            = new std::vector<number>;
  eleGsfPz            = new std::vector<number>;
  eleGsfPtAtVtx       = new std::vector<number>;
  eleGsfEtaAtVtx      = new std::vector<number>;
  eleGsfPhiAtVtx      = new std::vector<number>;
  eleGsfPxAtVtx       = new std::vector<number>;
  eleGsfPyAtVtx       = new std::vector<number>;
  eleGsfPzAtVtx       = new std::vector<number>;
  eleGsfChi2          = new std::vector<number>;
  eleGsfNdof          = new std::vector<number>;
  eleGsfLostHits      = new std::vector<int   >;
  eleGsfValidHits     = new std::vector<int   >;
  eleGsfPVtx          = new std::vector<int   >;
  eleGsfQuality       = new std::vector<int   >;
  eleGsfHitPattern    = new std::vector<int   >;
  eleGsfLayPattern    = new std::vector<int   >;
  eleGsfPtError       = new std::vector<number>;
  eleGsfDxy           = new std::vector<number>;
  eleGsfDxyPV         = new std::vector<number>;
  eleGsfDz            = new std::vector<number>;
  eleGsfDzPV          = new std::vector<number>;
  eleGsfExy           = new std::vector<number>;
//   eleGsfExyPV         = new std::vector<number>;
  eleGsfEz            = new std::vector<number>;
//   eleGsfEzPV          = new std::vector<number>;
  eleCtfPt            = new std::vector<number>;
  eleCtfEta           = new std::vector<number>;
  eleCtfPhi           = new std::vector<number>;
  eleCtfPx            = new std::vector<number>;
  eleCtfPy            = new std::vector<number>;
  eleCtfPz            = new std::vector<number>;
  eleCtfTrk           = new std::vector<int   >;
  eleCtfChi2          = new std::vector<number>;
  eleCtfNdof          = new std::vector<number>;
  eleCtfLostHits      = new std::vector<int   >;
  eleCtfValidHits     = new std::vector<int   >;
  eleCtfPVtx          = new std::vector<int   >;
  eleCtfQuality       = new std::vector<int   >;
  eleCtfHitPattern    = new std::vector<int   >;
  eleCtfLayPattern    = new std::vector<int   >;
  eleCtfPtError       = new std::vector<number>;
  eleCtfDxy           = new std::vector<number>;
  eleCtfDxyPV         = new std::vector<number>;
  eleCtfDz            = new std::vector<number>;
  eleCtfDzPV          = new std::vector<number>;
  eleCtfExy           = new std::vector<number>;
//   eleCtfExyPV         = new std::vector<number>;
  eleCtfEz            = new std::vector<number>;
//   eleCtfEzPV          = new std::vector<number>;
  eleProps            = new std::vector<int   >;
  elePartIso          = new std::vector<number>;
  eleEcalIso          = new std::vector<number>;
  eleHcalIso          = new std::vector<number>;
  eleTrackIso         = new std::vector<number>;
  eleClass            = new std::vector<int   >;
  eleCtfGsfOv         = new std::vector<number>;
  eleDEtaElClTkAtCalo = new std::vector<number>;
  eleDEtaSeSlTkAtCalo = new std::vector<number>;
  eleDEtaSeClTkAtVtx  = new std::vector<number>;
  eleDEtaSuClTkAtVtx  = new std::vector<number>;
  eleDPhiElClTkAtCalo = new std::vector<number>;
  eleDPhiSeClTkAtCalo = new std::vector<number>;
  eleDPhiSuClTkAtVtx  = new std::vector<number>;
  eleEcalE            = new std::vector<number>;
  eleEcalEE           = new std::vector<number>;
  eleEElClOvPOut      = new std::vector<number>;
  eleESeClOvP         = new std::vector<number>;
  eleESeClOvPOut      = new std::vector<number>;
  eleESuClOvP         = new std::vector<number>;
  eleFBrem            = new std::vector<number>;
  eleHcalOvEcal       = new std::vector<number>;
  eleF5x5E1x5         = new std::vector<number>;
  eleF5x5E2x5B        = new std::vector<number>;
  eleF5x5E2x5L        = new std::vector<number>;
  eleF5x5E2x5M        = new std::vector<number>;
  eleF5x5E2x5R        = new std::vector<number>;
  eleF5x5E2x5T        = new std::vector<number>;
  eleF5x5E5x5         = new std::vector<number>;
  eleF5x5EB           = new std::vector<number>;
  eleF5x5EL           = new std::vector<number>;
  eleF5x5ER           = new std::vector<number>;
  eleF5x5ET           = new std::vector<number>;
  eleF5x5HcalOvEcal   = new std::vector<number>;
  eleF5x5HcalOvEcalV  = new std::vector<number>;
  eleF5x5R9           = new std::vector<number>;
  eleF5x5SEtaEta      = new std::vector<number>;
  eleF5x5SIetaIeta    = new std::vector<number>;
  eleF5x5SIphiIphi    = new std::vector<number>;
  eleSCEtaWidth       = new std::vector<number>;
  eleSCPhiWidth       = new std::vector<number>;
  eleSCPreshowerE     = new std::vector<number>;
  eleSCRawE           = new std::vector<number>;
  eleSCEta            = new std::vector<number>;
  eleConvDcot         = new std::vector<number>;
  eleConvDist         = new std::vector<number>;
  eleConvFlags        = new std::vector<int   >;
  eleConvRadius       = new std::vector<number>;
  eleConvVtxProb      = new std::vector<number>;

  // taus
  use_taus      = false;
  use_taus_sphe = true;
  use_taus_cart = false;
  nTaus = 0;
  tauPt           = new std::vector<number>;
  tauEta          = new std::vector<number>;
  tauPhi          = new std::vector<number>;
  tauPx           = new std::vector<number>;
  tauPy           = new std::vector<number>;
  tauPz           = new std::vector<number>;
  tauE            = new std::vector<number>;
  tauCharge       = new std::vector<int   >;
  tauTrg          = new std::vector<int   >;

  // jets
  use_jets      = false;
  use_jets_sphe = true;
  use_jets_cart = false;
  nJets = 0;
  jetPt           = new std::vector<number>;
  jetEta          = new std::vector<number>;
  jetPhi          = new std::vector<number>;
  jetPx           = new std::vector<number>;
  jetPy           = new std::vector<number>;
  jetPz           = new std::vector<number>;
  jetE            = new std::vector<number>;
  jetCSV          = new std::vector<number>;
  jetTCHE         = new std::vector<number>;
  jetTrg          = new std::vector<int   >;
  jetPF           = new std::vector<bool  >;
  jetNDau         = new std::vector<int   >;
  jetNHF          = new std::vector<number>;
  jetNEF          = new std::vector<number>;
  jetCHF          = new std::vector<number>;
  jetCEF          = new std::vector<number>;
  jetNCH          = new std::vector<number>;

  // jet tags
  use_tags = false;
  nTags = 0;
  tagJet  = new std::vector<int   >;
  tagType = new std::vector<int   >;
  tagProb = new std::vector<number>;

  // user info
  use_info = false;
  nUserInfo = 0;
  useObjType   = new std::vector<int   >;
  useObjIndex  = new std::vector<int   >;
  useInfoType  = new std::vector<int   >;
  useInfoValue = new std::vector<number>;

  // particle flow
  use_pflow      = false;
  use_pflow_sphe = true;
  use_pflow_cart = false;
  nPF = 0;
  pfcPt           = new std::vector<number>;
  pfcEta          = new std::vector<number>;
  pfcPhi          = new std::vector<number>;
  pfcPx           = new std::vector<number>;
  pfcPy           = new std::vector<number>;
  pfcPz           = new std::vector<number>;
  pfcE            = new std::vector<number>;
  pfcCharge       = new std::vector<int   >;
  pfcJet          = new std::vector<int   >;
  pfcTrk          = new std::vector<int   >;
  pfcHasTrackDetails = new std::vector<int>;

  // tracks
  use_tracks      = false;
  use_tracks_sphe = true;
  use_tracks_cart = false;
  nTracks = 0;
  trkPt           = new std::vector<number>;
  trkEta          = new std::vector<number>;
  trkPhi          = new std::vector<number>;
  trkPx           = new std::vector<number>;
  trkPy           = new std::vector<number>;
  trkPz           = new std::vector<number>;
  trkCharge       = new std::vector<int   >;
  trkType         = new std::vector<int   >;
  trkNext         = new std::vector<int   >;
  trkPFC          = new std::vector<int   >;
  trkJet          = new std::vector<int   >;
  trkPVtx         = new std::vector<int   >;
//  trkSVtx         = new std::vector<int   >;
  trkQuality      = new std::vector<int   >;
  trkHitPattern   = new std::vector<int   >;
  trkLayPattern   = new std::vector<int   >;
  trkMissingInnerHits = new std::vector<int   >;
  trkNormChi2     = new std::vector<number>;
  trkPtError      = new std::vector<number>;
  trkDxy          = new std::vector<number>;
  trkDz           = new std::vector<number>;
  trkExy          = new std::vector<number>;
  trkEz           = new std::vector<number>;
  trkVtxPx        = new std::vector<number>;
  trkVtxPy        = new std::vector<number>;
  trkVtxPz        = new std::vector<number>;
  trkQOverPError  = new std::vector<number>;
  trkLambdaError  = new std::vector<number>;
  trkPhiError     = new std::vector<number>;
  trkDszError     = new std::vector<number>;
  trkEtaError     = new std::vector<number>;
  trkDxyDszCov    = new std::vector<number>;
  trkLambdaDszCov = new std::vector<number>;
  trkPhiDxyCov    = new std::vector<number>;

  // track perigee parameters
  use_trkper = false;
  nTrkPer = 0;
  tppTrk     = new std::vector<int   >;
  tppQop     = new std::vector<number>;
  tppLam     = new std::vector<number>;
  tppPhi     = new std::vector<number>;
  tppDxy     = new std::vector<number>;
  tppDsz     = new std::vector<number>;
  tppSQopQop = new std::vector<number>;
  tppSQopLam = new std::vector<number>;
  tppSQopPhi = new std::vector<number>;
  tppSQopDxy = new std::vector<number>;
  tppSQopDsz = new std::vector<number>;
  tppSLamLam = new std::vector<number>;
  tppSLamPhi = new std::vector<number>;
  tppSLamDxy = new std::vector<number>;
  tppSLamDsz = new std::vector<number>;
  tppSPhiPhi = new std::vector<number>;
  tppSPhiDxy = new std::vector<number>;
  tppSPhiDsz = new std::vector<number>;
  tppSDxyDxy = new std::vector<number>;
  tppSDxyDsz = new std::vector<number>;
  tppSDszDsz = new std::vector<number>;

  // primary vertices
  use_pvts =  false;
  nPVTotal   = 0;
  nPVertices = 0;
  pvtX            = new std::vector<number>;
  pvtY            = new std::vector<number>;
  pvtZ            = new std::vector<number>;
  pvtSxx          = new std::vector<number>;
  pvtSyy          = new std::vector<number>;
  pvtSzz          = new std::vector<number>;
  pvtSxy          = new std::vector<number>;
  pvtSxz          = new std::vector<number>;
  pvtSyz          = new std::vector<number>;
//  pvtCovariance   = new std::vector<
//                        std::vector<number>
//                                          >;
  pvtNTracks      = new std::vector<int   >;
//  pvtNormChi2     = new std::vector<number>;
  pvtChi2         = new std::vector<number>;
  pvtNDOF         = new std::vector<number>;
  pvtBadQuality   = new std::vector<int   >;

  // secondary vertices
  use_svts = false;
  nSVertices = 0;
  svtX          = new std::vector<number>;
  svtY          = new std::vector<number>;
  svtZ          = new std::vector<number>;
  svtSxx        = new std::vector<number>;
  svtSyy        = new std::vector<number>;
  svtSzz        = new std::vector<number>;
  svtSxy        = new std::vector<number>;
  svtSxz        = new std::vector<number>;
  svtSyz        = new std::vector<number>;
  svtPVtx       = new std::vector<int   >;
  svtDirX       = new std::vector<number>;
  svtDirY       = new std::vector<number>;
  svtDirZ       = new std::vector<number>;
  svtType       = new std::vector<int   >;
  svtNTracks    = new std::vector<int   >;
  svtChi2       = new std::vector<number>;
  svtNDOF       = new std::vector<number>;
  svtMass       = new std::vector<number>;
  svtDist2D     = new std::vector<number>;
  svtSigma2D    = new std::vector<number>;
  svtDist3D     = new std::vector<number>;
  svtSigma3D    = new std::vector<number>;
  svtJet        = new std::vector<int   >;
  svtBadQuality = new std::vector<int   >;

  // composite particle vertices
  use_vsub = false;
  nCompVts = 0;
  subPart = new std::vector<int   >;
  subSVtx = new std::vector<int   >;
  subMass = new std::vector<number>;

  // impact parameters
  use_tkips = false;
  nTkIPs = 0;
  tipTrk          = new std::vector<int   >;
  tipSVtx         = new std::vector<int   >;
  tipMass         = new std::vector<number>;
  tipDxy          = new std::vector<number>;
  tipDz           = new std::vector<number>;
  tipExy          = new std::vector<number>;
  tipEz           = new std::vector<number>;

  // momenta at vertices
  use_vtxps      = false;
  use_vtxps_sphe = true;
  use_vtxps_cart = false;
  nVtxPs = 0;
  tvpTip          = new std::vector<int   >;
  tvpPt           = new std::vector<number>;
  tvpEta          = new std::vector<number>;
  tvpPhi          = new std::vector<number>;
  tvpPx           = new std::vector<number>;
  tvpPy           = new std::vector<number>;
  tvpPz           = new std::vector<number>;

  // PU weight
  use_puwgt = false;

  // gen particles
  use_gen      = false;
  use_gen_sphe = true;
  use_gen_cart = false;
  nGenP = 0;
  genId           = new std::vector<int   >;
  genStatus       = new std::vector<int   >;
  genMother       = new std::vector<int   >;
  genPartner      = new std::vector<int   >;
  genPt           = new std::vector<number>;
  genEta          = new std::vector<number>;
  genPhi          = new std::vector<number>;
  genPx           = new std::vector<number>;
  genPy           = new std::vector<number>;
  genPz           = new std::vector<number>;
  genE            = new std::vector<number>;
  genCharge       = new std::vector<int   >;
  genMass         = new std::vector<number>;
//  genJet          = new std::vector<int   >;
  genVx           = new std::vector<number>;
  genVy           = new std::vector<number>;
  genVz           = new std::vector<number>;

  // gen jets
  use_gpj      = false;
  use_gpj_sphe = true;
  use_gpj_cart = false;
  nGenJets = 0;
  gpjPt           = new std::vector<number>;
  gpjEta          = new std::vector<number>;
  gpjPhi          = new std::vector<number>;
  gpjPx           = new std::vector<number>;
  gpjPy           = new std::vector<number>;
  gpjPz           = new std::vector<number>;
  gpjE            = new std::vector<number>;
  gpjNDau         = new std::vector<int   >;
  gpjReco         = new std::vector<int   >;

}

PDNtupleData::~PDNtupleData() {
}


