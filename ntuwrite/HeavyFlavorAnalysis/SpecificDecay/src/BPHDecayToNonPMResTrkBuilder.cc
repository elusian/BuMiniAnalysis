/*
 *  See header file for a description of this class.
 *
 *  \author Paolo Ronchese INFN Padova
 *
 */

//-----------------------
// This Class' Header --
//-----------------------
#include "HeavyFlavorAnalysis/SpecificDecay/interface/BPHDecayToNonPMResTrkBuilder.h"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "HeavyFlavorAnalysis/RecoDecay/interface/BPHRecoBuilder.h"
#include "HeavyFlavorAnalysis/RecoDecay/interface/BPHPlusMinusCandidate.h"
#include "HeavyFlavorAnalysis/RecoDecay/interface/BPHRecoCandidate.h"
#include "HeavyFlavorAnalysis/SpecificDecay/interface/BPHParticleNeutralVeto.h"

//---------------
// C++ Headers --
//---------------
using namespace std;

//-------------------
// Initializations --
//-------------------


//----------------
// Constructors --
//----------------
BPHDecayToNonPMResTrkBuilder::BPHDecayToNonPMResTrkBuilder( const edm::EventSetup& es,
    const std::string& resName, double resMass, double resWidth,
    const std::vector<BPHRecoConstCandPtr>& resCollection,
    const std::string& trkName, double trkMass, double trkSigma,
    const BPHRecoBuilder::BPHGenericCollection*  trkCollection ):
 BPHDecayNonPMConstrainedBuilder( es, resName, resMass, resWidth, resCollection ),
 tName( trkName ),
 tMass( trkMass ),
 tSigma( trkSigma ),
 tCollection(  trkCollection ),
 tknVeto( new BPHParticleNeutralVeto ),
 ptSel  ( new BPHParticlePtSelect (   0.0 ) ),
 etaSel ( new BPHParticleEtaSelect( 100.0 ) ) {
}

//--------------
// Destructor --
//--------------
BPHDecayToNonPMResTrkBuilder::~BPHDecayToNonPMResTrkBuilder() {
  delete tknVeto;
  delete   ptSel;
  delete  etaSel;
}

//--------------
// Operations --
//--------------
vector<BPHRecoConstCandPtr> BPHDecayToNonPMResTrkBuilder::build() {

  if ( updated ) return recList;

  recList.clear();

  BPHRecoBuilder brb( *evSetup );
  brb.setMinPDiffererence( minPDiff );
  brb.add( rName, *rCollection );
  brb.add( tName,  tCollection, tMass, tSigma );
  if ( resoSel->getMassMax() > 0.0 )
  brb.filter( rName, *resoSel );
  brb.filter( tName, *tknVeto );
  if (  ptSel->getPtMin   () >= 0.0 )
  brb.filter( tName, *  ptSel );
  if (  etaSel->getEtaMax () >= 0.0 )
  brb.filter( tName, * etaSel );

  if ( massSel->getMassMax() >= 0.0 )
  brb.filter( *massSel );
  if ( chi2Sel->getProbMin() >= 0.0 )
  brb.filter( *chi2Sel );
  if ( mFitSel->getMassMax() >= 0.0 )
  brb.filter( *mFitSel );

  recList = BPHRecoCandidate::build( brb );
  updated = true;
  return recList;

}

/// set cuts
void BPHDecayToNonPMResTrkBuilder::setTrkPtMin( double pt ) {
  updated = false;
  ptSel->setPtMin( pt );
  return;
}


void BPHDecayToNonPMResTrkBuilder::setTrkEtaMax( double eta ) {
  updated = false;
  etaSel->setEtaMax( eta );
  return;
}

