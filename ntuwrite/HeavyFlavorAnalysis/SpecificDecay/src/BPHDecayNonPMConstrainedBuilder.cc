/*
 *  See header file for a description of this class.
 *
 *  \author Paolo Ronchese INFN Padova
 *
 */

//-----------------------
// This Class' Header --
//-----------------------
#include "HeavyFlavorAnalysis/SpecificDecay/interface/BPHDecayNonPMConstrainedBuilder.h"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "HeavyFlavorAnalysis/RecoDecay/interface/BPHRecoBuilder.h"
#include "HeavyFlavorAnalysis/RecoDecay/interface/BPHPlusMinusCandidate.h"
#include "HeavyFlavorAnalysis/RecoDecay/interface/BPHRecoCandidate.h"

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
BPHDecayNonPMConstrainedBuilder::BPHDecayNonPMConstrainedBuilder( const edm::EventSetup& es,
    const std::string& resName, double resMass, double resWidth,
    const std::vector<BPHRecoConstCandPtr>& resCollection ):
 BPHDecayGenericBuilder( es, new BPHMassFitSelect( resName, resMass, resWidth,
                                                   -2.0e+06, -1.0e+06 ) ),
 rName( resName ),
 rMass( resMass ),
 rWidth( resWidth ),
 rCollection( &resCollection ),
 resoSel( new BPHMassSelect( -2.0e+06, -1.0e+06 ) ),
 massConstr( true ) {
}

//--------------
// Destructor --
//--------------
BPHDecayNonPMConstrainedBuilder::~BPHDecayNonPMConstrainedBuilder() {
  delete resoSel;
}

//--------------
// Operations --
//--------------
/// set cuts
void BPHDecayNonPMConstrainedBuilder::setResMassMin( double m ) {
  updated = false;
  resoSel->setMassMin( m );
  return;
}


void BPHDecayNonPMConstrainedBuilder::setResMassMax( double m ) {
  updated = false;
  resoSel->setMassMax( m );
  return;
}


void BPHDecayNonPMConstrainedBuilder::setResMassRange( double mMin, double mMax ) {
  updated = false;
  resoSel->setMassMin( mMin );
  resoSel->setMassMax( mMax );
  return;
}


void BPHDecayNonPMConstrainedBuilder::setConstr( bool flag ) {
  updated = false;
  if ( flag == massConstr ) return;
  double mMin = mFitSel->getMassMin();
  double mMax = mFitSel->getMassMax();
  delete mFitSel;
  massConstr = flag;
  if ( massConstr ) mFitSel = new BPHMassFitSelect    ( rName, rMass, rWidth,
                                                        mMin, mMax );
  else              mFitSel = new BPHMassFitSelect    ( mMin, mMax );
  return;
}

