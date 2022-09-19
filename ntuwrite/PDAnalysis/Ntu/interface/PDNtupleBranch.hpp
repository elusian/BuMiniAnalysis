#include <iostream>

template <class T>
PDNtupleBranch<T>::PDNtupleBranch() {
  T::treeName = "PDtree";
}


template <class T>
PDNtupleBranch<T>::~PDNtupleBranch() {
}


template <class T>
void PDNtupleBranch<T>::initTree() {

  std::cout << "PDNtupleBranch<T>::initTree()" << std::endl;

  // define branches
  // use local function "setBranch" (defined in the utility)
  // in place of ROOT "Branch" & "SetBranchAddress"

  if ( use_hltlist   ) setBranches_hltlist  ();
  if ( use_hlts      ) setBranches_hlts     ();
  if ( use_hlto      ) setBranches_hlto     ();
  if ( use_hltm      ) setBranches_hltm     ();
  if ( use_bspot     ) setBranches_bspot    ();
  if ( use_met       ) setBranches_met      ();
  if ( use_muons     ) setBranches_muons    ();
  if ( use_electrons ) setBranches_electrons();
  if ( use_taus      ) setBranches_taus     ();
  if ( use_jets      ) setBranches_jets     ();
  if ( use_tags      ) setBranches_tags     ();
  if ( use_info      ) setBranches_userinfo ();
  if ( use_pflow     ) setBranches_pflow    ();
  if ( use_tracks    ) setBranches_tracks   ();
  if ( use_trkper    ) setBranches_trkper   ();
  if ( use_pvts      ) setBranches_pVertices();
  if ( use_svts      ) setBranches_sVertices();
  if ( use_vsub      ) setBranches_vsub     ();
  if ( use_tkips     ) setBranches_tkips    ();
  if ( use_vtxps     ) setBranches_vtxps    ();
  if ( use_puwgt     ) setBranches_puwgt    ();
  if ( use_gen       ) setBranches_gen      ();
  if ( use_gpj       ) setBranches_gpj      ();

  return;

}


template <class T>
void PDNtupleBranch<T>::setBranches_hltlist() {
  // hlt path names
  int size = 8192;
  this->setBranch( "nHLTPaths"                , &  nHLTPaths     ,
                   "nHLTPaths/I"              , &b_nHLTPaths      );
  this->setBranch( "hltCode"                  , &  hltCode      ,
                   size, 99                   , &b_hltCode       );
  this->setBranch( "hltName"                  , &  hltName      ,
                   size, 99                   , &b_hltName       );
  return;
}


template <class T>
void PDNtupleBranch<T>::setBranches_hlts() {
  // hlt status
  int size = 8192;
  this->setBranch( "nHLTStatus"               , &  nHLTStatus     ,
                   "nHLTStatus/I"             , &b_nHLTStatus      );
  this->setBranch( "hltPath"                  , &  hltPath        ,
                   size, 99                   , &b_hltPath         );
  this->setBranch( "hltVersion"               , &  hltVersion     ,
                   size, 99                   , &b_hltVersion      );
  this->setBranch( "hltScale"                 , &  hltScale       ,
                   size, 99                   , &b_hltScale        );
  this->setBranch( "hltRun"                   , &  hltRun         ,
                   size, 99                   , &b_hltRun          );
  this->setBranch( "hltAccept"                , &  hltAccept      ,
                   size, 99                   , &b_hltAccept       );
  return;
}


template <class T>
void PDNtupleBranch<T>::setBranches_hlto() {
  // hlt objects
  int size = 8192;
  this->setBranch( "nHLTObjects"              , &  nHLTObjects    ,
                   "nHLTObjects/I"            , &b_nHLTObjects     );
  this->setBranch( "hltObjType"               , &  hltObjType     ,
                   size, 99                   , &b_hltObjType      );
  if ( use_hlto_sphe ) {
  this->setBranch( "hltPt"                    , &  hltPt          ,
                   size, 99                   , &b_hltPt           );
  this->setBranch( "hltEta"                   , &  hltEta         ,
                   size, 99                   , &b_hltEta          );
  this->setBranch( "hltPhi"                   , &  hltPhi         ,
                   size, 99                   , &b_hltPhi          );
  }
  if ( use_hlto_cart ) {
  this->setBranch( "hltPx"                    , &  hltPx          ,
                   size, 99                   , &b_hltPx           );
  this->setBranch( "hltPy"                    , &  hltPy          ,
                   size, 99                   , &b_hltPy           );
  this->setBranch( "hltPz"                    , &  hltPz          ,
                   size, 99                   , &b_hltPz           );
  }
  this->setBranch( "hltE"                     , &  hltE           ,
                   size, 99                   , &b_hltE            );
  return;
}


template <class T>
void PDNtupleBranch<T>::setBranches_hltm() {
  // hlt matches
  int size = 8192;
  this->setBranch( "nHLTMatches"              , &  nHLTMatches    ,
                   "nHLTMatches/I"            , &b_nHLTMatches     );
  this->setBranch( "hltPathIndex"             , &  hltPathIndex   ,
                   size, 99                   , &b_hltPathIndex    );
  this->setBranch( "hltObjectIndex"           , &  hltObjectIndex ,
                   size, 99                   , &b_hltObjectIndex  );
  this->setBranch( "hltFilter"                , &  hltFilter      ,
                   size, 99                   , &b_hltFilter       );
  return;
}


template <class T>
void PDNtupleBranch<T>::setBranches_bspot() {
  // beam spot
  this->setBranch( "bwX"                        , &  bwX           ,
             bCat( "bwX",           numberType ), &b_bwX            );
  this->setBranch( "bwY"                        , &  bwY           ,
             bCat( "bwY",           numberType ), &b_bwY            );
  this->setBranch( "bwXError"                   , &  bwXError      ,
             bCat( "bwXError",      numberType ), &b_bwXError       );
  this->setBranch( "bwYError"                   , &  bwYError      ,
             bCat( "bwYError",      numberType ), &b_bwYError       );
  this->setBranch( "bsX"                        , &  bsX           ,
             bCat( "bsX",           numberType ), &b_bsX            );
  this->setBranch( "bsY"                        , &  bsY           ,
             bCat( "bsY",           numberType ), &b_bsY            );
  this->setBranch( "bsZ"                        , &  bsZ           ,
             bCat( "bsZ",           numberType ), &b_bsZ            );
  this->setBranch( "bsXError"                   , &  bsXError      ,
             bCat( "bsXError",      numberType ), &b_bsXError       );
  this->setBranch( "bsYError"                   , &  bsYError      ,
             bCat( "bsYError",      numberType ), &b_bsYError       );
  this->setBranch( "bsZError"                   , &  bsZError      ,
             bCat( "bsZError",      numberType ), &b_bsZError       );
//  this->setBranch( "bsSxy"                      , &  bsSxy         ,
//             bCat( "bsSxy",         numberType ), &b_bsSxy          );
//  this->setBranch( "bsSxz"                      , &  bsSxz         ,
//             bCat( "bsSxz",         numberType ), &b_bsSxz          );
//  this->setBranch( "bsSyz"                      , &  bsSyz         ,
//             bCat( "bsSyz",         numberType ), &b_bsSyz          );
  this->setBranch( "bsdXdZ"                     , &  bsdXdZ        ,
             bCat( "bsdXdZ",        numberType ), &b_bsdXdZ         );
  this->setBranch( "bsdYdZ"                     , &  bsdYdZ        ,
             bCat( "bsdYdZ",        numberType ), &b_bsdYdZ         );
  this->setBranch( "bsdXdZError"                , &  bsdXdZError   ,
             bCat( "bsdXdZError",   numberType ), &b_bsdXdZError    );
  this->setBranch( "bsdYdZError"                , &  bsdYdZError   ,
             bCat( "bsdYdZError",   numberType ), &b_bsdYdZError    );
  this->setBranch( "bsSigmaZ"                   , &  bsSigmaZ      ,
             bCat( "bsSigmaZ",      numberType ), &b_bsSigmaZ       );
  this->setBranch( "bsSigmaZError"              , &  bsSigmaZError ,
             bCat( "bsSigmaZError", numberType ), &b_bsSigmaZError  );
  return;
}


template <class T>
void PDNtupleBranch<T>::setBranches_met() {
  // met
  this->setBranch( "mEt"                      , &  mEt            ,
             bCat( "mEt",         numberType ), &b_mEt             );
  this->setBranch( "mEx"                      , &  mEx            ,
             bCat( "mEx",         numberType ), &b_mEx             );
  this->setBranch( "mEy"                      , &  mEy            ,
             bCat( "mEy",         numberType ), &b_mEy             );
  this->setBranch( "mEsign"                   , &  mEsign         ,
             bCat( "mEsign",      numberType ), &b_mEsign          );
  return;
}


template <class T>
void PDNtupleBranch<T>::setBranches_muons() {
  // muons
  int size = 8192;
  this->setBranch( "nMuons"                   , &  nMuons         ,
                   "nMuons/I"                 , &b_nMuons          );
  if ( use_muons_sphe ) {
  this->setBranch( "muoPt"                    , &  muoPt          ,
                   size, 99                   , &b_muoPt           );
  this->setBranch( "muoEta"                   , &  muoEta         ,
                   size, 99                   , &b_muoEta          );
  this->setBranch( "muoPhi"                   , &  muoPhi         ,
                   size, 99                   , &b_muoPhi          );
  }
  if ( use_muons_cart ) {
  this->setBranch( "muoPx"                    , &  muoPx          ,
                   size, 99                   , &b_muoPx           );
  this->setBranch( "muoPy"                    , &  muoPy          ,
                   size, 99                   , &b_muoPy           );
  this->setBranch( "muoPz"                    , &  muoPz          ,
                   size, 99                   , &b_muoPz           );
  }
  this->setBranch( "muoE"                     , &  muoE           ,
                   size, 99                   , &b_muoE            );
  this->setBranch( "muoTrk"                   , &  muoTrk         ,
                   size, 99                   , &b_muoTrk          );
  this->setBranch( "muoTrg"                   , &  muoTrg         ,
                   size, 99                   , &b_muoTrg          );
  this->setBranch( "muoCharge"                , &  muoCharge      ,
                   size, 99                   , &b_muoCharge       );
  this->setBranch( "muoChaIso"                , &  muoChaIso      ,
                   size, 99                   , &b_muoChaIso       );
  this->setBranch( "muoNeuIso"                , &  muoNeuIso      ,
                   size, 99                   , &b_muoNeuIso       );
  this->setBranch( "muoPhoIso"                , &  muoPhoIso      ,
                   size, 99                   , &b_muoPhoIso       );
  this->setBranch( "muoSumCHpt"               , &  muoSumCHpt     ,
                   size, 99                   , &b_muoSumCHpt      );
  this->setBranch( "muoSumCPpt"               , &  muoSumCPpt     ,
                   size, 99                   , &b_muoSumCPpt      );
  this->setBranch( "muoSumNHet"               , &  muoSumNHet     ,
                   size, 99                   , &b_muoSumNHet      );
  this->setBranch( "muoSumPHet"               , &  muoSumPHet     ,
                   size, 99                   , &b_muoSumPHet      );
  this->setBranch( "muoSumPUpt"               , &  muoSumPUpt     ,
                   size, 99                   , &b_muoSumPUpt      );
  this->setBranch( "muoNumMatches"            , &  muoNumMatches  ,
                   size, 99                   , &b_muoNumMatches   );
  this->setBranch( "muoType"                  , &  muoType        ,
                   size, 99                   , &b_muoType         );
  this->setBranch( "muoDb"                    , &  muoDb          ,
                   size, 99                   , &b_muoDb           );
  this->setBranch( "muoChi2LP"                , &  muoChi2LP      ,
                   size, 99                   , &b_muoChi2LP       );
  this->setBranch( "muoTrkKink"               , &  muoTrkKink     ,
                   size, 99                   , &b_muoTrkKink      );
  this->setBranch( "muoSegmComp"              , &  muoSegmComp    ,
                   size, 99                   , &b_muoSegmComp     );
//  this->setBranch( "muoNumValidHits"          , &  muoNumValidHits,
//                   size, 99                   , &b_muoNumValidHits );
//  this->setBranch( "muoNormChi2"              , &  muoNormChi2    ,
//                   size, 99                   , &b_muoNormChi2     );
//  this->setBranch( "muoNumMuHits"             , &  muoNumMuHits   ,
//                   size, 99                   , &b_muoNumMuHits    );
//  this->setBranch( "muoNumPixHits"            , &  muoNumPixHits  ,
//                   size, 99                   , &b_muoNumPixHits   );
//  this->setBranch( "muoNumTkHits"             , &  muoNumTkHits   ,
//                   size, 99                   , &b_muoNumTkHits    );

  this->setBranch( "muoChi2LM"                , &  muoChi2LM           ,
                   size, 99                   , &b_muoChi2LM            );
  this->setBranch( "muoTrkRelChi2"            , &  muoTrkRelChi2       ,
                   size, 99                   , &b_muoTrkRelChi2        );
  this->setBranch( "muoGlbTrackTailProb"      , &  muoGlbTrackTailProb ,
                   size, 99                   , &b_muoGlbTrackTailProb  );
  this->setBranch( "muoGlbKinkFinderLOG"      , &  muoGlbKinkFinderLOG ,
                   size, 99                   , &b_muoGlbKinkFinderLOG  );
  this->setBranch( "muoGlbDeltaEtaPhi"        , &  muoGlbDeltaEtaPhi   ,
                   size, 99                   , &b_muoGlbDeltaEtaPhi    );
  this->setBranch( "muoStaRelChi2"            , &  muoStaRelChi2       ,
                   size, 99                   , &b_muoStaRelChi2        );
  this->setBranch( "muoTimeAtIpInOut"         , &  muoTimeAtIpInOut    ,
                   size, 99                   , &b_muoTimeAtIpInOut     );
  this->setBranch( "muoTimeAtIpInOutErr"      , &  muoTimeAtIpInOutErr ,
                   size, 99                   , &b_muoTimeAtIpInOutErr  );
//  this->setBranch( "muoInnerChi2"             , &  muoInnerChi2        ,
//                   size, 99                   , &b_muoInnerChi2         );
  this->setBranch( "muoIValFrac"              , &  muoIValFrac         ,
                   size, 99                   , &b_muoIValFrac          );
//  this->setBranch( "muoValPixHits"            , &  muoValPixHits       ,
//                   size, 99                   , &b_muoValPixHits        );
//  this->setBranch( "muoNTrkVHits"             , &  muoNTrkVHits        ,
//                   size, 99                   , &b_muoNTrkVHits         );
//  this->setBranch( "muoOuterChi2"             , &  muoOuterChi2        ,
//                   size, 99                   , &b_muoOuterChi2         );
//  this->setBranch( "muoGNchi2"                , &  muoGNchi2           ,
//                   size, 99                   , &b_muoGNchi2            );
//  this->setBranch( "muoVMuHits"               , &  muoVMuHits          ,
//                   size, 99                   , &b_muoVMuHits           );
//  this->setBranch( "muoQprod"                 , &  muoQprod            ,
//                   size, 99                   , &b_muoQprod             );
//  this->setBranch( "muoLWH"                   , &  muoLWH              ,
//                   size, 99                   , &b_muoLWH               );
  this->setBranch( "muoVMuonHitComb"          , &  muoVMuonHitComb     ,
                   size, 99                   , &b_muoVMuonHitComb      );

  return;
}


template <class T>
void PDNtupleBranch<T>::setBranches_electrons() {
  // electrons
  int size = 8192;
  this->setBranch( "nElectrons"               , &  nElectrons     ,
                   "nElectrons/I"             , &b_nElectrons      );
  if ( use_electrons_sphe ) {
  this->setBranch( "elePt"                    , &  elePt          ,
                   size, 99                   , &b_elePt           );
  this->setBranch( "eleEta"                   , &  eleEta         ,
                   size, 99                   , &b_eleEta          );
  this->setBranch( "elePhi"                   , &  elePhi         ,
                   size, 99                   , &b_elePhi          );
  this->setBranch( "eleGsfPt"                 , &  eleGsfPt           ,
                   size, 99                   , &b_eleGsfPt            );
  this->setBranch( "eleGsfEta"                , &  eleGsfEta          ,
                   size, 99                   , &b_eleGsfEta           );
  this->setBranch( "eleGsfPhi"                , &  eleGsfPhi          ,
                   size, 99                   , &b_eleGsfPhi           );
  this->setBranch( "eleGsfPtAtVtx"            , &  eleGsfPtAtVtx      ,
                   size, 99                   , &b_eleGsfPtAtVtx       );
  this->setBranch( "eleGsfEtaAtVtx"           , &  eleGsfEtaAtVtx     ,
                   size, 99                   , &b_eleGsfEtaAtVtx      );
  this->setBranch( "eleGsfPhiAtVtx"           , &  eleGsfPhiAtVtx     ,
                   size, 99                   , &b_eleGsfPhiAtVtx      );
  this->setBranch( "eleCtfPt"                 , &  eleCtfPt           ,
                   size, 99                   , &b_eleCtfPt            );
  this->setBranch( "eleCtfEta"                , &  eleCtfEta          ,
                   size, 99                   , &b_eleCtfEta           );
  this->setBranch( "eleCtfPhi"                , &  eleCtfPhi          ,
                   size, 99                   , &b_eleCtfPhi           );
  }
  if ( use_electrons_cart ) {
  this->setBranch( "elePx"                    , &  elePx          ,
                   size, 99                   , &b_elePx           );
  this->setBranch( "elePy"                    , &  elePy          ,
                   size, 99                   , &b_elePy           );
  this->setBranch( "elePz"                    , &  elePz          ,
                   size, 99                   , &b_elePz           );
  this->setBranch( "eleGsfPx"                 , &  eleGsfPx           ,
                   size, 99                   , &b_eleGsfPx            );
  this->setBranch( "eleGsfPy"                 , &  eleGsfPy           ,
                   size, 99                   , &b_eleGsfPy            );
  this->setBranch( "eleGsfPz"                 , &  eleGsfPz           ,
                   size, 99                   , &b_eleGsfPz            );
  this->setBranch( "eleGsfPxAtVtx"            , &  eleGsfPxAtVtx      ,
                   size, 99                   , &b_eleGsfPxAtVtx       );
  this->setBranch( "eleGsfPyAtVtx"            , &  eleGsfPyAtVtx      ,
                   size, 99                   , &b_eleGsfPyAtVtx       );
  this->setBranch( "eleGsfPzAtVtx"            , &  eleGsfPzAtVtx      ,
                   size, 99                   , &b_eleGsfPzAtVtx       );
  this->setBranch( "eleCtfPx"                 , &  eleCtfPx           ,
                   size, 99                   , &b_eleCtfPx            );
  this->setBranch( "eleCtfPy"                 , &  eleCtfPy           ,
                   size, 99                   , &b_eleCtfPy            );
  this->setBranch( "eleCtfPz"                 , &  eleCtfPz           ,
                   size, 99                   , &b_eleCtfPz            );
  }
  this->setBranch( "eleE"                     , &  eleE           ,
                   size, 99                   , &b_eleE            );
  this->setBranch( "eleTrk"                   , &  eleTrk         ,
                   size, 99                   , &b_eleTrk          );
  this->setBranch( "eleTrg"                   , &  eleTrg         ,
                   size, 99                   , &b_eleTrg          );
  this->setBranch( "eleCharge"                , &  eleCharge      ,
                   size, 99                   , &b_eleCharge       );
  this->setBranch( "eleChaIso"                , &  eleChaIso      ,
                   size, 99                   , &b_eleChaIso       );
  this->setBranch( "eleNeuIso"                , &  eleNeuIso      ,
                   size, 99                   , &b_eleNeuIso       );
  this->setBranch( "elePhoIso"                , &  elePhoIso      ,
                   size, 99                   , &b_elePhoIso       );
  this->setBranch( "elePCHIso"                , &  elePCHIso      ,
                   size, 99                   , &b_elePCHIso       );
  this->setBranch( "eleSumCHpt"               , &  eleSumCHpt         ,
                   size, 99                   , &b_eleSumCHpt          );
  this->setBranch( "eleSumCPpt"               , &  eleSumCPpt         ,
                   size, 99                   , &b_eleSumCPpt          );
  this->setBranch( "eleSumNHet"               , &  eleSumNHet         ,
                   size, 99                   , &b_eleSumNHet          );
  this->setBranch( "eleSumPHet"               , &  eleSumPHet         ,
                   size, 99                   , &b_eleSumPHet          );
  this->setBranch( "eleSumPUpt"               , &  eleSumPUpt         ,
                   size, 99                   , &b_eleSumPUpt          );
  this->setBranch( "eleAbsEta"                , &  eleAbsEta      ,
                   size, 99                   , &b_eleAbsEta       );
  this->setBranch( "eleAEff"                  , &  eleAEff        ,
                   size, 99                   , &b_eleAEff         );
  this->setBranch( "eleID"                    , &  eleID          ,
                   size, 99                   , &b_eleID           );
  this->setBranch( "eleConsCha"               , &  eleConsCha     ,
                   size, 99                   , &b_eleConsCha      );
  this->setBranch( "eleEBEEGap"               , &  eleEBEEGap     ,
                   size, 99                   , &b_eleEBEEGap      );
  this->setBranch( "eleDb"                    , &  eleDb          ,
                   size, 99                   , &b_eleDb           );
  this->setBranch( "eleGsfChi2"               , &  eleGsfChi2         ,
                   size, 99                   , &b_eleGsfChi2          );
  this->setBranch( "eleGsfNdof"               , &  eleGsfNdof         ,
                   size, 99                   , &b_eleGsfNdof          );
  this->setBranch( "eleGsfLostHits"           , &  eleGsfLostHits     ,
                   size, 99                   , &b_eleGsfLostHits      );
  this->setBranch( "eleGsfValidHits"          , &  eleGsfValidHits    ,
                   size, 99                   , &b_eleGsfValidHits     );
  this->setBranch( "eleGsfPVtx"               , &  eleGsfPVtx         ,
                   size, 99                   , &b_eleGsfPVtx          );
  this->setBranch( "eleGsfQuality"            , &  eleGsfQuality      ,
                   size, 99                   , &b_eleGsfQuality       );
  this->setBranch( "eleGsfHitPattern"         , &  eleGsfHitPattern   ,
                   size, 99                   , &b_eleGsfHitPattern    );
  this->setBranch( "eleGsfLayPattern"         , &  eleGsfLayPattern   ,
                   size, 99                   , &b_eleGsfLayPattern    );
  this->setBranch( "eleGsfPtError"            , &  eleGsfPtError      ,
                   size, 99                   , &b_eleGsfPtError       );
  this->setBranch( "eleGsfDxy"                , &  eleGsfDxy          ,
                   size, 99                   , &b_eleGsfDxy           );
  this->setBranch( "eleGsfDxyPV"              , &  eleGsfDxyPV        ,
                   size, 99                   , &b_eleGsfDxyPV         );
  this->setBranch( "eleGsfDz"                 , &  eleGsfDz           ,
                   size, 99                   , &b_eleGsfDz            );
  this->setBranch( "eleGsfDzPV"               , &  eleGsfDzPV         ,
                   size, 99                   , &b_eleGsfDzPV          );
  this->setBranch( "eleGsfExy"                , &  eleGsfExy          ,
                   size, 99                   , &b_eleGsfExy           );
//   this->setBranch( "eleGsfExyPV"              , &  eleGsfExyPV        ,
//                    size, 99                   , &b_eleGsfExyPV         );
  this->setBranch( "eleGsfEz"                 , &  eleGsfEz           ,
                   size, 99                   , &b_eleGsfEz            );
//   this->setBranch( "eleGsfEzPV"               , &  eleGsfEzPV         ,
//                    size, 99                   , &b_eleGsfEzPV          );
  this->setBranch( "eleCtfTrk"                , &  eleCtfTrk          ,
                   size, 99                   , &b_eleCtfTrk           );
  this->setBranch( "eleCtfChi2"               , &  eleCtfChi2         ,
                   size, 99                   , &b_eleCtfChi2          );
  this->setBranch( "eleCtfNdof"               , &  eleCtfNdof         ,
                   size, 99                   , &b_eleCtfNdof          );
  this->setBranch( "eleCtfLostHits"           , &  eleCtfLostHits     ,
                   size, 99                   , &b_eleCtfLostHits      );
  this->setBranch( "eleCtfValidHits"          , &  eleCtfValidHits    ,
                   size, 99                   , &b_eleCtfValidHits     );
  this->setBranch( "eleCtfPVtx"               , &  eleCtfPVtx         ,
                   size, 99                   , &b_eleCtfPVtx          );
  this->setBranch( "eleCtfQuality"            , &  eleCtfQuality      ,
                   size, 99                   , &b_eleCtfQuality       );
  this->setBranch( "eleCtfHitPattern"         , &  eleCtfHitPattern   ,
                   size, 99                   , &b_eleCtfHitPattern    );
  this->setBranch( "eleCtfLayPattern"         , &  eleCtfLayPattern   ,
                   size, 99                   , &b_eleCtfLayPattern    );
  this->setBranch( "eleCtfPtError"            , &  eleCtfPtError      ,
                   size, 99                   , &b_eleCtfPtError       );
  this->setBranch( "eleCtfDxy"                , &  eleCtfDxy          ,
                   size, 99                   , &b_eleCtfDxy           );
  this->setBranch( "eleCtfDxyPV"              , &  eleCtfDxyPV        ,
                   size, 99                   , &b_eleCtfDxyPV         );
  this->setBranch( "eleCtfDz"                 , &  eleCtfDz           ,
                   size, 99                   , &b_eleCtfDz            );
  this->setBranch( "eleCtfDzPV"               , &  eleCtfDzPV         ,
                   size, 99                   , &b_eleCtfDzPV          );
  this->setBranch( "eleCtfExy"                , &  eleCtfExy          ,
                   size, 99                   , &b_eleCtfExy           );
//   this->setBranch( "eleCtfExyPV"              , &  eleCtfExyPV     ,
//                    size, 99                   , &b_eleCtfExyPV      );
  this->setBranch( "eleCtfEz"                 , &  eleCtfEz           ,
                   size, 99                   , &b_eleCtfEz            );
//   this->setBranch( "eleCtfEzPV"               , &  eleCtfEzPV         ,
//                    size, 99                   , &b_eleCtfEzPV          );
  this->setBranch( "eleProps"                 , &  eleProps           ,
                   size, 99                   , &b_eleProps            );
  this->setBranch( "elePartIso"               , &  elePartIso         ,
                   size, 99                   , &b_elePartIso          );
  this->setBranch( "eleEcalIso"               , &  eleEcalIso         ,
                   size, 99                   , &b_eleEcalIso          );
  this->setBranch( "eleHcalIso"               , &  eleHcalIso         ,
                   size, 99                   , &b_eleHcalIso          );
  this->setBranch( "eleTrackIso"              , &  eleTrackIso        ,
                   size, 99                   , &b_eleTrackIso         );
  this->setBranch( "eleClass"                 , &  eleClass           ,
                   size, 99                   , &b_eleClass            );
  this->setBranch( "eleCtfGsfOv"              , &  eleCtfGsfOv        ,
                   size, 99                   , &b_eleCtfGsfOv         );
  this->setBranch( "eleDEtaElClTkAtCalo"      , &  eleDEtaElClTkAtCalo,
                   size, 99                   , &b_eleDEtaElClTkAtCalo );
  this->setBranch( "eleDEtaSeSlTkAtCalo"      , &  eleDEtaSeSlTkAtCalo,
                   size, 99                   , &b_eleDEtaSeSlTkAtCalo );
  this->setBranch( "eleDEtaSeClTkAtVtx"       , &  eleDEtaSeClTkAtVtx ,
                   size, 99                   , &b_eleDEtaSeClTkAtVtx  );
  this->setBranch( "eleDEtaSuClTkAtVtx"       , &  eleDEtaSuClTkAtVtx ,
                   size, 99                   , &b_eleDEtaSuClTkAtVtx  );
  this->setBranch( "eleDPhiElClTkAtCalo"      , &  eleDPhiElClTkAtCalo,
                   size, 99                   , &b_eleDPhiElClTkAtCalo );
  this->setBranch( "eleDPhiSeClTkAtCalo"      , &  eleDPhiSeClTkAtCalo,
                   size, 99                   , &b_eleDPhiSeClTkAtCalo );
  this->setBranch( "eleDPhiSuClTkAtVtx"       , &  eleDPhiSuClTkAtVtx ,
                   size, 99                   , &b_eleDPhiSuClTkAtVtx  );
  this->setBranch( "eleEcalE"                 , &  eleEcalE           ,
                   size, 99                   , &b_eleEcalE            );
  this->setBranch( "eleEcalEE"                , &  eleEcalEE          ,
                   size, 99                   , &b_eleEcalEE           );
  this->setBranch( "eleEElClOvPOut"           , &  eleEElClOvPOut     ,
                   size, 99                   , &b_eleEElClOvPOut      );
  this->setBranch( "eleESeClOvP"              , &  eleESeClOvP        ,
                   size, 99                   , &b_eleESeClOvP         );
  this->setBranch( "eleESeClOvPOut"           , &  eleESeClOvPOut     ,
                   size, 99                   , &b_eleESeClOvPOut      );
  this->setBranch( "eleESuClOvP"              , &  eleESuClOvP        ,
                   size, 99                   , &b_eleESuClOvP         );
  this->setBranch( "eleFBrem"                 , &  eleFBrem           ,
                   size, 99                   , &b_eleFBrem            );
  this->setBranch( "eleHcalOvEcal"            , &  eleHcalOvEcal      ,
                   size, 99                   , &b_eleHcalOvEcal       );
  this->setBranch( "eleF5x5E1x5"              , &  eleF5x5E1x5        ,
                   size, 99                   , &b_eleF5x5E1x5         );
  this->setBranch( "eleF5x5E2x5B"             , &  eleF5x5E2x5B       ,
                   size, 99                   , &b_eleF5x5E2x5B        );
  this->setBranch( "eleF5x5E2x5L"             , &  eleF5x5E2x5L       ,
                   size, 99                   , &b_eleF5x5E2x5L        );
  this->setBranch( "eleF5x5E2x5M"             , &  eleF5x5E2x5M       ,
                   size, 99                   , &b_eleF5x5E2x5M        );
  this->setBranch( "eleF5x5E2x5R"             , &  eleF5x5E2x5R       ,
                   size, 99                   , &b_eleF5x5E2x5R        );
  this->setBranch( "eleF5x5E2x5T"             , &  eleF5x5E2x5T       ,
                   size, 99                   , &b_eleF5x5E2x5T        );
  this->setBranch( "eleF5x5E5x5"              , &  eleF5x5E5x5        ,
                   size, 99                   , &b_eleF5x5E5x5         );
  this->setBranch( "eleF5x5EB"                , &  eleF5x5EB          ,
                   size, 99                   , &b_eleF5x5EB           );
  this->setBranch( "eleF5x5EL"                , &  eleF5x5EL          ,
                   size, 99                   , &b_eleF5x5EL           );
  this->setBranch( "eleF5x5ER"                , &  eleF5x5ER          ,
                   size, 99                   , &b_eleF5x5ER           );
  this->setBranch( "eleF5x5ET"                , &  eleF5x5ET          ,
                   size, 99                   , &b_eleF5x5ET           );
  this->setBranch( "eleF5x5HcalOvEcal"        , &  eleF5x5HcalOvEcal  ,
                   size, 99                   , &b_eleF5x5HcalOvEcal   );
  this->setBranch( "eleF5x5HcalOvEcalV"       , &  eleF5x5HcalOvEcalV ,
                   size, 99                   , &b_eleF5x5HcalOvEcalV  );
  this->setBranch( "eleF5x5R9"                , &  eleF5x5R9          ,
                   size, 99                   , &b_eleF5x5R9           );
  this->setBranch( "eleF5x5SEtaEta"           , &  eleF5x5SEtaEta     ,
                   size, 99                   , &b_eleF5x5SEtaEta      );
  this->setBranch( "eleF5x5SIetaIeta"         , &  eleF5x5SIetaIeta   ,
                   size, 99                   , &b_eleF5x5SIetaIeta    );
  this->setBranch( "eleF5x5SIphiIphi"         , &  eleF5x5SIphiIphi   ,
                   size, 99                   , &b_eleF5x5SIphiIphi    );
  this->setBranch( "eleSCEtaWidth"            , &  eleSCEtaWidth      ,
                   size, 99                   , &b_eleSCEtaWidth       );
  this->setBranch( "eleSCPhiWidth"            , &  eleSCPhiWidth      ,
                   size, 99                   , &b_eleSCPhiWidth       );
  this->setBranch( "eleSCPreshowerE"          , &  eleSCPreshowerE    ,
                   size, 99                   , &b_eleSCPreshowerE     );
  this->setBranch( "eleSCRawE"                , &  eleSCRawE          ,
                   size, 99                   , &b_eleSCRawE           );
  this->setBranch( "eleSCEta"                 , &  eleSCEta           ,
                   size, 99                   , &b_eleSCEta            );
  this->setBranch( "eleConvDcot"              , &  eleConvDcot        ,
                   size, 99                   , &b_eleConvDcot         );
  this->setBranch( "eleConvDist"              , &  eleConvDist        ,
                   size, 99                   , &b_eleConvDist         );
  this->setBranch( "eleConvFlags"             , &  eleConvFlags       ,
                   size, 99                   , &b_eleConvFlags        );
  this->setBranch( "eleConvRadius"            , &  eleConvRadius      ,
                   size, 99                   , &b_eleConvRadius       );
  this->setBranch( "eleConvVtxProb"           , &  eleConvVtxProb     ,
                   size, 99                   , &b_eleConvVtxProb      );
  return;
}


template <class T>
void PDNtupleBranch<T>::setBranches_taus() {
  // taus
  int size = 8192;
  this->setBranch( "nTaus"                    , &  nTaus          ,
                   "nTaus/I"                  , &b_nTaus           );
  if ( use_taus_sphe ) {
  this->setBranch( "tauPt"                    , &  tauPt          ,
                   size, 99                   , &b_tauPt           );
  this->setBranch( "tauEta"                   , &  tauEta         ,
                   size, 99                   , &b_tauEta          );
  this->setBranch( "tauPhi"                   , &  tauPhi         ,
                   size, 99                   , &b_tauPhi          );
  }
  if ( use_taus_cart ) {
  this->setBranch( "tauPx"                    , &  tauPx          ,
                   size, 99                   , &b_tauPx           );
  this->setBranch( "tauPy"                    , &  tauPy          ,
                   size, 99                   , &b_tauPy           );
  this->setBranch( "tauPz"                    , &  tauPz          ,
                   size, 99                   , &b_tauPz           );
  }
  this->setBranch( "tauE"                     , &  tauE           ,
                   size, 99                   , &b_tauE            );
  this->setBranch( "tauCharge"                , &  tauCharge      ,
                   size, 99                   , &b_tauCharge       );
  this->setBranch( "tauTrg"                   , &  tauTrg         ,
                   size, 99                   , &b_tauTrg          );
  return;
}


template <class T>
void PDNtupleBranch<T>::setBranches_jets() {
  // jets
  int size = 8192;
  this->setBranch( "nJets"                    , &  nJets          ,
                   "nJets/I"                  , &b_nJets           );
  if ( use_jets_sphe ) {
  this->setBranch( "jetPt"                    , &  jetPt          ,
                   size, 99                   , &b_jetPt           );
  this->setBranch( "jetEta"                   , &  jetEta         ,
                   size, 99                   , &b_jetEta          );
  this->setBranch( "jetPhi"                   , &  jetPhi         ,
                   size, 99                   , &b_jetPhi          );
  }
  if ( use_jets_cart ) {
  this->setBranch( "jetPx"                    , &  jetPx          ,
                   size, 99                   , &b_jetPx           );
  this->setBranch( "jetPy"                    , &  jetPy          ,
                   size, 99                   , &b_jetPy           );
  this->setBranch( "jetPz"                    , &  jetPz          ,
                   size, 99                   , &b_jetPz           );
  }
  this->setBranch( "jetE"                     , &  jetE           ,
                   size, 99                   , &b_jetE            );
  this->setBranch( "jetCSV"                   , &  jetCSV         ,
                   size, 99                   , &b_jetCSV          );
  this->setBranch( "jetTCHE"                  , &  jetTCHE        ,
                   size, 99                   , &b_jetTCHE         );
  this->setBranch( "jetTrg"                   , &  jetTrg         ,
                   size, 99                   , &b_jetTrg          );
  this->setBranch( "jetPF"                    , &  jetPF          ,
                   size, 99                   , &b_jetPF           );
  this->setBranch( "jetNDau"                  , &  jetNDau        ,
                   size, 99                   , &b_jetNDau         );
  this->setBranch( "jetNHF"                   , &  jetNHF         ,
                   size, 99                   , &b_jetNHF          );
  this->setBranch( "jetNEF"                   , &  jetNEF         ,
                   size, 99                   , &b_jetNEF          );
  this->setBranch( "jetCHF"                   , &  jetCHF         ,
                   size, 99                   , &b_jetCHF          );
  this->setBranch( "jetCEF"                   , &  jetCEF         ,
                   size, 99                   , &b_jetCEF          );
  this->setBranch( "jetNCH"                   , &  jetNCH         ,
                   size, 99                   , &b_jetNCH          );
  return;
}


template <class T>
void PDNtupleBranch<T>::setBranches_tags() {
  // jets
  int size = 8192;
  this->setBranch( "nTags"                    , &  nTags          ,
                   "nTags/I"                  , &b_nTags           );
  this->setBranch( "tagJet"                   , &  tagJet         ,
                   size, 99                   , &b_tagJet          );
  this->setBranch( "tagType"                  , &  tagType        ,
                   size, 99                   , &b_tagType         );
  this->setBranch( "tagProb"                  , &  tagProb        ,
                   size, 99                   , &b_tagProb         );
}


template <class T>
void PDNtupleBranch<T>::setBranches_userinfo() {
  int size = 8192;
  this->setBranch( "nUserInfo"                , &  nUserInfo      ,
                   "nUserInfo/I"              , &b_nUserInfo       );
  this->setBranch( "useObjType"               , &  useObjType     ,
                   size, 99                   , &b_useObjType      );
  this->setBranch( "useObjIndex"              , &  useObjIndex    ,
                   size, 99                   , &b_useObjIndex     );
  this->setBranch( "useInfoType"              , &  useInfoType    ,
                   size, 99                   , &b_useInfoType     );
  this->setBranch( "useInfoValue"             , &  useInfoValue   ,
                   size, 99                   , &b_useInfoValue    );
  return;
}


template <class T>
void PDNtupleBranch<T>::setBranches_pflow() {
  // particle flow
  int size = 32768;
  this->setBranch( "nPF"                      , &  nPF            ,
                   "nPF/I"                    , &b_nPF             );
  if ( use_pflow_sphe ) {
  this->setBranch( "pfcPt"                    , &  pfcPt          ,
                   size, 99                   , &b_pfcPt           );
  this->setBranch( "pfcEta"                   , &  pfcEta         ,
                   size, 99                   , &b_pfcEta          );
  this->setBranch( "pfcPhi"                   , &  pfcPhi         ,
                   size, 99                   , &b_pfcPhi          );
  }
  if ( use_pflow_cart ) {
  this->setBranch( "pfcPx"                    , &  pfcPx          ,
                   size, 99                   , &b_pfcPx           );
  this->setBranch( "pfcPy"                    , &  pfcPy          ,
                   size, 99                   , &b_pfcPy           );
  this->setBranch( "pfcPz"                    , &  pfcPz          ,
                   size, 99                   , &b_pfcPz           );
  }
  this->setBranch( "pfcE"                     , &  pfcE           ,
                   size, 99                   , &b_pfcE            );
  this->setBranch( "pfcCharge"                , &  pfcCharge      ,
                   size, 99                   , &b_pfcCharge       );
  this->setBranch( "pfcJet"                   , &  pfcJet         ,
                   size, 99                   , &b_pfcJet          );
  this->setBranch( "pfcTrk"                   , &  pfcTrk         ,
                   size, 99                   , &b_pfcTrk          );
  this->setBranch( "pfcHasTrackDetais"        , &  pfcHasTrackDetails,
                   size, 99                   , &b_pfcHasTrackDetails);
  return;
}


template <class T>
void PDNtupleBranch<T>::setBranches_tracks() {
  // tracks
  int size = 32768;
  this->setBranch( "nTracks"                  , &  nTracks     ,
                   "nTracks/I"                , &b_nTracks      );
  if ( use_tracks_sphe ) {
  this->setBranch( "trkPt"                    , &  trkPt          ,
                   size, 99                   , &b_trkPt           );
  this->setBranch( "trkEta"                   , &  trkEta         ,
                   size, 99                   , &b_trkEta          );
  this->setBranch( "trkPhi"                   , &  trkPhi         ,
                   size, 99                   , &b_trkPhi          );
  }
  if ( use_tracks_cart ) {
  this->setBranch( "trkPx"                    , &  trkPx          ,
                   size, 99                   , &b_trkPx           );
  this->setBranch( "trkPy"                    , &  trkPy          ,
                   size, 99                   , &b_trkPy           );
  this->setBranch( "trkPz"                    , &  trkPz          ,
                   size, 99                   , &b_trkPz           );
  }
  this->setBranch( "trkCharge"                , &  trkCharge      ,
                   size, 99                   , &b_trkCharge       );
  this->setBranch( "trkType"                  , &  trkType        ,
                   size, 99                   , &b_trkType         );
  this->setBranch( "trkNext"                  , &  trkNext        ,
                   size, 99                   , &b_trkNext         );
  this->setBranch( "trkPFC"                   , &  trkPFC         ,
                   size, 99                   , &b_trkPFC          );
  this->setBranch( "trkJet"                   , &  trkJet         ,
                   size, 99                   , &b_trkJet          );
  this->setBranch( "trkPVtx"                  , &  trkPVtx        ,
                   size, 99                   , &b_trkPVtx         );
//  this->setBranch( "trkSVtx"                  , &  trkSVtx        ,
//                   size, 99                   , &b_trkSVtx         );
  this->setBranch( "trkQuality"               , &  trkQuality     ,
                   size, 99                   , &b_trkQuality      );
  this->setBranch( "trkHitPattern"            , &  trkHitPattern  ,
                   size, 99                   , &b_trkHitPattern   );
  this->setBranch( "trkLayPattern"            , &  trkLayPattern  ,
                   size, 99                   , &b_trkLayPattern   );
  this->setBranch( "trkMissingInnerHits"      , &  trkMissingInnerHits  ,
                   size, 99                   , &b_trkMissingInnerHits   );
  this->setBranch( "trkNormChi2"              , &  trkNormChi2    ,
                   size, 99                   , &b_trkNormChi2     );
  this->setBranch( "trkPtError"               , &  trkPtError     ,
                   size, 99                   , &b_trkPtError      );
  this->setBranch( "trkDxy"                   , &  trkDxy         ,
                   size, 99                   , &b_trkDxy          );
  this->setBranch( "trkDz"                    , &  trkDz          ,
                   size, 99                   , &b_trkDz           );
  this->setBranch( "trkExy"                   , &  trkExy         ,
                   size, 99                   , &b_trkExy          );
  this->setBranch( "trkEz"                    , &  trkEz          ,
                   size, 99                   , &b_trkEz           );
  this->setBranch( "trkVtxPx"                 , &  trkVtxPx        ,
                   size, 99                   , &b_trkVtxPx         );
  this->setBranch( "trkVtxPy"                 , &  trkVtxPy        ,
                   size, 99                   , &b_trkVtxPy         );
  this->setBranch( "trkVtxPz"                 , &  trkVtxPz        ,
                   size, 99                   , &b_trkVtxPz         );
  this->setBranch( "trkQOverPError"                 , &  trkQOverPError        ,
                   size, 99                   , &b_trkQOverPError         );
  this->setBranch( "trkLambdaError"                 , &  trkLambdaError        ,
                   size, 99                   , &b_trkLambdaError         );
  this->setBranch( "trkPhiError"                 , &  trkPhiError        ,
                   size, 99                   , &b_trkPhiError         );
  this->setBranch( "trkDszError"                 , &  trkDszError        ,
                   size, 99                   , &b_trkDszError         );
  this->setBranch( "trkEtaError"                 , &  trkEtaError        ,
                   size, 99                   , &b_trkEtaError         );
  this->setBranch( "trkDxyDszCov"                 , &  trkDxyDszCov        ,
                   size, 99                   , &b_trkDxyDszCov         );
  this->setBranch( "trkLambdaDszCov"                 , &  trkLambdaDszCov        ,
                   size, 99                   , &b_trkLambdaDszCov         );
  this->setBranch( "trkPhiDxyCov"                 , &  trkPhiDxyCov        ,
                   size, 99                   , &b_trkPhiDxyCov         );
  return;
}


template <class T>
void PDNtupleBranch<T>::setBranches_trkper() {
  // track perigee parameters
  int size = 32768;
  this->setBranch( "nTrkPer"                  , &  nTrkPer     ,
                   "nTrkPer/I"                , &b_nTrkPer      );
  this->setBranch( "tppTrk"                   , &  tppTrk      ,
                   size, 99                   , &b_tppTrk       );
  this->setBranch( "tppQop"                   , &  tppQop      ,
                   size, 99                   , &b_tppQop       );
  this->setBranch( "tppLam"                   , &  tppLam      ,
                   size, 99                   , &b_tppLam       );
  this->setBranch( "tppPhi"                   , &  tppPhi      ,
                   size, 99                   , &b_tppPhi       );
  this->setBranch( "tppDxy"                   , &  tppDxy      ,
                   size, 99                   , &b_tppDxy       );
  this->setBranch( "tppDsz"                   , &  tppDsz      ,
                   size, 99                   , &b_tppDsz       );
  this->setBranch( "tppSQopQop"               , &  tppSQopQop  ,
                   size, 99                   , &b_tppSQopQop   );
  this->setBranch( "tppSQopLam"               , &  tppSQopLam  ,
                   size, 99                   , &b_tppSQopLam   );
  this->setBranch( "tppSQopPhi"               , &  tppSQopPhi  ,
                   size, 99                   , &b_tppSQopPhi   );
  this->setBranch( "tppSQopDxy"               , &  tppSQopDxy  ,
                   size, 99                   , &b_tppSQopDxy   );
  this->setBranch( "tppSQopDsz"               , &  tppSQopDsz  ,
                   size, 99                   , &b_tppSQopDsz   );
  this->setBranch( "tppSLamLam"               , &  tppSLamLam  ,
                   size, 99                   , &b_tppSLamLam   );
  this->setBranch( "tppSLamPhi"               , &  tppSLamPhi  ,
                   size, 99                   , &b_tppSLamPhi   );
  this->setBranch( "tppSLamDxy"               , &  tppSLamDxy  ,
                   size, 99                   , &b_tppSLamDxy   );
  this->setBranch( "tppSLamDsz"               , &  tppSLamDsz  ,
                   size, 99                   , &b_tppSLamDsz   );
  this->setBranch( "tppSPhiPhi"               , &  tppSPhiPhi  ,
                   size, 99                   , &b_tppSPhiPhi   );
  this->setBranch( "tppSPhiDxy"               , &  tppSPhiDxy  ,
                   size, 99                   , &b_tppSPhiDxy   );
  this->setBranch( "tppSPhiDsz"               , &  tppSPhiDsz  ,
                   size, 99                   , &b_tppSPhiDsz   );
  this->setBranch( "tppSDxyDxy"               , &  tppSDxyDxy  ,
                   size, 99                   , &b_tppSDxyDxy   );
  this->setBranch( "tppSDxyDsz"               , &  tppSDxyDsz  ,
                   size, 99                   , &b_tppSDxyDsz   );
  this->setBranch( "tppSDszDsz"               , &  tppSDszDsz  ,
                   size, 99                   , &b_tppSDszDsz   );
  return;
}


template <class T>
void PDNtupleBranch<T>::setBranches_pVertices() {
  // primary vertices
  int size = 8192;
  this->setBranch( "nPVTotal"                 , &  nPVTotal       ,
                   "nPVTotal/I"               , &b_nPVTotal        );
  this->setBranch( "nPVertices"               , &  nPVertices     ,
                   "nPVertices/I"             , &b_nPVertices      );
  this->setBranch( "pvtX"                     , &  pvtX           ,
                   size, 99                   , &b_pvtX            );
  this->setBranch( "pvtY"                     , &  pvtY           ,
                   size, 99                   , &b_pvtY            );
  this->setBranch( "pvtZ"                     , &  pvtZ           ,
                   size, 99                   , &b_pvtZ            );
  this->setBranch( "pvtSxx"                   , &  pvtSxx         ,
                   size, 99                   , &b_pvtSxx          );
  this->setBranch( "pvtSyy"                   , &  pvtSyy         ,
                   size, 99                   , &b_pvtSyy          );
  this->setBranch( "pvtSzz"                   , &  pvtSzz         ,
                   size, 99                   , &b_pvtSzz          );
  this->setBranch( "pvtSxy"                   , &  pvtSxy         ,
                   size, 99                   , &b_pvtSxy          );
  this->setBranch( "pvtSxz"                   , &  pvtSxz         ,
                   size, 99                   , &b_pvtSxz          );
  this->setBranch( "pvtSyz"                   , &  pvtSyz         ,
                   size, 99                   , &b_pvtSyz          );
//  this->setBranch( "pvtCovariance"            , &  pvtCovariance  ,
//                   size, 99                   , &b_pvtCovariance   );
  this->setBranch( "pvtNTracks"               , &  pvtNTracks     ,
                   size, 99                   , &b_pvtNTracks      );
//  this->setBranch( "pvtNormChi2"              , &  pvtNormChi2    ,
//                   size, 99                   , &b_pvtNormChi2     );
  this->setBranch( "pvtChi2"                  , &  pvtChi2        ,
                   size, 99                   , &b_pvtChi2         );
  this->setBranch( "pvtNDOF"                  , &  pvtNDOF        ,
                   size, 99                   , &b_pvtNDOF         );
  this->setBranch( "pvtBadQuality"            , &  pvtBadQuality  ,
                   size, 99                   , &b_pvtBadQuality   );
  return;
}


template <class T>
void PDNtupleBranch<T>::setBranches_sVertices() {
  // secondary vertices
  int size = 8192;
  this->setBranch( "nSVertices"               , &  nSVertices     ,
                   "nSVertices/I"             , &b_nSVertices      );
  this->setBranch( "svtX"                     , &  svtX           ,
                   size, 99                   , &b_svtX            );
  this->setBranch( "svtY"                     , &  svtY           ,
                   size, 99                   , &b_svtY            );
  this->setBranch( "svtZ"                     , &  svtZ           ,
                   size, 99                   , &b_svtZ            );
  this->setBranch( "svtSxx"                   , &  svtSxx         ,
                   size, 99                   , &b_svtSxx          );
  this->setBranch( "svtSyy"                   , &  svtSyy         ,
                   size, 99                   , &b_svtSyy          );
  this->setBranch( "svtSzz"                   , &  svtSzz         ,
                   size, 99                   , &b_svtSzz          );
  this->setBranch( "svtSxy"                   , &  svtSxy         ,
                   size, 99                   , &b_svtSxy          );
  this->setBranch( "svtSxz"                   , &  svtSxz         ,
                   size, 99                   , &b_svtSxz          );
  this->setBranch( "svtSyz"                   , &  svtSyz         ,
                   size, 99                   , &b_svtSyz          );
  this->setBranch( "svtPVtx"                  , &  svtPVtx        ,
                   size, 99                   , &b_svtPVtx         );
  this->setBranch( "svtDirX"                  , &  svtDirX        ,
                   size, 99                   , &b_svtDirX         );
  this->setBranch( "svtDirY"                  , &  svtDirY        ,
                   size, 99                   , &b_svtDirY         );
  this->setBranch( "svtDirZ"                  , &  svtDirZ        ,
                   size, 99                   , &b_svtDirZ         );
  this->setBranch( "svtType"                  , &  svtType        ,
                   size, 99                   , &b_svtType         );
  this->setBranch( "svtNTracks"               , &  svtNTracks     ,
                   size, 99                   , &b_svtNTracks      );
  this->setBranch( "svtChi2"                  , &  svtChi2        ,
                   size, 99                   , &b_svtChi2         );
  this->setBranch( "svtNDOF"                  , &  svtNDOF        ,
                   size, 99                   , &b_svtNDOF         );
  this->setBranch( "svtMass"                  , &  svtMass        ,
                   size, 99                   , &b_svtMass         );
  this->setBranch( "svtDist2D"                , &  svtDist2D      ,
                   size, 99                   , &b_svtDist2D       );
  this->setBranch( "svtSigma2D"               , &  svtSigma2D     ,
                   size, 99                   , &b_svtSigma2D      );
  this->setBranch( "svtDist3D"                , &  svtDist3D      ,
                   size, 99                   , &b_svtDist3D       );
  this->setBranch( "svtSigma3D"               , &  svtSigma3D     ,
                   size, 99                   , &b_svtSigma3D      );
  this->setBranch( "svtJet"                   , &  svtJet         ,
                   size, 99                   , &b_svtJet          );
  this->setBranch( "svtBadQuality"            , &  svtBadQuality  ,
                   size, 99                   , &b_svtBadQuality   );
  return;
}


template <class T>
void PDNtupleBranch<T>::setBranches_vsub() {
  // composite particle vertices
  int size = 8192;
  this->setBranch( "nCompVts"                 , &  nCompVts       ,
                   "nCompVts/I"               , &b_nCompVts        );
  this->setBranch( "subPart"                  , &  subPart        ,
                   size, 99                   , &b_subPart         );
  this->setBranch( "subSVtx"                  , &  subSVtx        ,
                   size, 99                   , &b_subSVtx         );
  this->setBranch( "subMass"                  , &  subMass        ,
                   size, 99                   , &b_subMass         );
  return;
}


template <class T>
void PDNtupleBranch<T>::setBranches_tkips() {
  // refitted vertices tracks
  int size = 8192;
  this->setBranch( "nTkIPs"                   , &  nTkIPs         ,
                   "nTkIPs/I"                 , &b_nTkIPs          );
  this->setBranch( "tipTrk"                   , &  tipTrk         ,
                   size, 99                   , &b_tipTrk          );
  this->setBranch( "tipSVtx"                  , &  tipSVtx        ,
                   size, 99                   , &b_tipSVtx         );
  this->setBranch( "tipMass"                  , &  tipMass        ,
                   size, 99                   , &b_tipMass         );
  this->setBranch( "tipDxy"                   , &  tipDxy         ,
                   size, 99                   , &b_tipDxy          );
  this->setBranch( "tipDz"                    , &  tipDz          ,
                   size, 99                   , &b_tipDz           );
  this->setBranch( "tipExy"                   , &  tipExy         ,
                   size, 99                   , &b_tipExy          );
  this->setBranch( "tipEz"                    , &  tipEz          ,
                   size, 99                   , &b_tipEz           );
  return;
}


template <class T>
void PDNtupleBranch<T>::setBranches_vtxps() {
  // momenta at vertices
  int size = 32768;
  this->setBranch( "nVtxPs"                   , &  nVtxPs         ,
                   "nVtxPs/I"                 , &b_nVtxPs          );
  this->setBranch( "tvpTip"                   , &  tvpTip         ,
                   size, 99                   , &b_tvpTip          );
  if ( use_vtxps_sphe ) {
  this->setBranch( "tvpPt"                    , &  tvpPt          ,
                   size, 99                   , &b_tvpPt           );
  this->setBranch( "tvpEta"                   , &  tvpEta         ,
                   size, 99                   , &b_tvpEta          );
  this->setBranch( "tvpPhi"                   , &  tvpPhi         ,
                   size, 99                   , &b_tvpPhi          );
  }
  if ( use_vtxps_cart ) {
  this->setBranch( "tvpPx"                    , &  tvpPx          ,
                   size, 99                   , &b_tvpPx           );
  this->setBranch( "tvpPy"                    , &  tvpPy          ,
                   size, 99                   , &b_tvpPy           );
  this->setBranch( "tvpPz"                    , &  tvpPz          ,
                   size, 99                   , &b_tvpPz           );
  }
  return;
}


template <class T>
void PDNtupleBranch<T>::setBranches_puwgt() {
  this->setBranch( "puWeight"                 , &  puWeight       ,
             bCat( "puWeight",    numberType ), &b_puWeight        );
  return;
}


template <class T>
void PDNtupleBranch<T>::setBranches_gen() {
  // gen particles
  int size = 65536;
  this->setBranch( "nGenP"                    , &  nGenP          ,
                   "nGenP/I"                  , &b_nGenP           );
  this->setBranch( "genId"                    , &  genId          ,
                   size, 99                   , &b_genId           );
  this->setBranch( "genStatus"                , &  genStatus      ,
                   size, 99                   , &b_genStatus       );
  this->setBranch( "genMother"                , &  genMother      ,
                   size, 99                   , &b_genMother       );
  this->setBranch( "genPartner"               , &  genPartner     ,
                   size, 99                   , &b_genPartner      );
  if ( use_gen_sphe ) {
  this->setBranch( "genPt"                    , &  genPt          ,
                   size, 99                   , &b_genPt           );
  this->setBranch( "genEta"                   , &  genEta         ,
                   size, 99                   , &b_genEta          );
  this->setBranch( "genPhi"                   , &  genPhi         ,
                   size, 99                   , &b_genPhi          );
  }
  if ( use_gen_cart ) {
  this->setBranch( "genPx"                    , &  genPx          ,
                   size, 99                   , &b_genPx           );
  this->setBranch( "genPy"                    , &  genPy          ,
                   size, 99                   , &b_genPy           );
  this->setBranch( "genPz"                    , &  genPz          ,
                   size, 99                   , &b_genPz           );
  }
  this->setBranch( "genE"                     , &  genE           ,
                   size, 99                   , &b_genE            );
  this->setBranch( "genCharge"                , &  genCharge      ,
                   size, 99                   , &b_genCharge       );
  this->setBranch( "genMass"                  , &  genMass        ,
                   size, 99                   , &b_genMass         );
//  this->setBranch( "genJet"                   , &  genJet         ,
//                   size, 99                   , &b_genJet          );
  this->setBranch( "genVx"                    , &  genVx          ,
                   size, 99                   , &b_genVx           );
  this->setBranch( "genVy"                    , &  genVy          ,
                   size, 99                   , &b_genVy           );
  this->setBranch( "genVz"                    , &  genVz          ,
                   size, 99                   , &b_genVz           );
  return;
}


template <class T>
void PDNtupleBranch<T>::setBranches_gpj() {
  // gen jets
  int size = 65536;
  this->setBranch( "nGenJets"                 , &  nGenJets       ,
                   "nGenJets/I"               , &b_nGenJets        );
  if ( use_gpj_sphe ) {
  this->setBranch( "gpjPt"                    , &  gpjPt          ,
                   size, 99                   , &b_gpjPt           );
  this->setBranch( "gpjEta"                   , &  gpjEta         ,
                   size, 99                   , &b_gpjEta          );
  this->setBranch( "gpjPhi"                   , &  gpjPhi         ,
                   size, 99                   , &b_gpjPhi          );
  }
  if ( use_gpj_cart ) {
  this->setBranch( "gpjPx"                    , &  gpjPx          ,
                   size, 99                   , &b_gpjPx           );
  this->setBranch( "gpjPy"                    , &  gpjPy          ,
                   size, 99                   , &b_gpjPy           );
  this->setBranch( "gpjPz"                    , &  gpjPz          ,
                   size, 99                   , &b_gpjPz           );
  }
  this->setBranch( "gpjE"                     , &  gpjE           ,
                   size, 99                   , &b_gpjE            );
  this->setBranch( "gpjNDau"                  , &  gpjNDau        ,
                   size, 99                   , &b_gpjNDau         );
  this->setBranch( "gpjReco"                  , &  gpjReco        ,
                   size, 99                   , &b_gpjReco         );
  return;
}


template <class T>
const char* PDNtupleBranch<T>::bCat( const std::string& name,
                                     const std::string& type ) {
  static std::string s;
  s = name;
  s += "/";
  s += type;
  return s.c_str();
}

