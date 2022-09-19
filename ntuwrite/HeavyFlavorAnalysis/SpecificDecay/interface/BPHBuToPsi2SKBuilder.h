#ifndef HeavyFlavorAnalysis_SpecificDecay_BPHBuToPsi2SKBuilder_h
#define HeavyFlavorAnalysis_SpecificDecay_BPHBuToPsi2SKBuilder_h
/** \class BPHBuToPsi2SKBuilder
 *
 *  Description: 
 *     Class to build B+- to Psi2S K+- candidates
 *
 *  \author Paolo Ronchese INFN Padova
 *
 */

//----------------------
// Base Class Headers --
//----------------------
#include "HeavyFlavorAnalysis/SpecificDecay/interface/BPHDecayToNonPMResTrkBuilder.h"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
#include "HeavyFlavorAnalysis/RecoDecay/interface/BPHRecoBuilder.h"
#include "HeavyFlavorAnalysis/RecoDecay/interface/BPHRecoCandidate.h"
#include "HeavyFlavorAnalysis/SpecificDecay/interface/BPHParticleMasses.h"

#include "FWCore/Framework/interface/Event.h"

//---------------
// C++ Headers --
//---------------
#include <string>
#include <vector>

//              ---------------------
//              -- Class Interface --
//              ---------------------

class BPHBuToPsi2SKBuilder: public BPHDecayToNonPMResTrkBuilder {

 public:

  /** Constructor
   */
  BPHBuToPsi2SKBuilder( const edm::EventSetup& es,
      const std::vector<BPHRecoConstCandPtr>& psi2SCollection,
      const BPHRecoBuilder::BPHGenericCollection*  kaonCollection ):
   BPHDecayToNonPMResTrkBuilder( es,
                            "Psi2S", 
                            BPHParticleMasses::psi2Mass,
                            BPHParticleMasses::psi2MWidth, psi2SCollection,
                            "Kaon",
                            BPHParticleMasses::kaonMass,
                            BPHParticleMasses::kaonMSigma, kaonCollection ) {
    setResMassRange( 3.30, 4.00 );
    setTrkPtMin    (  0.7 );
    setTrkEtaMax   ( 10.0 );
    setMassRange   ( 3.50, 8.00 );
    setProbMin     ( 0.02 );
    setMassFitRange( 5.00, 6.00 );
    setConstr( true );
  }

  // deleted copy constructor and assignment operator
  BPHBuToPsi2SKBuilder           ( const BPHBuToPsi2SKBuilder& x ) = delete;
  BPHBuToPsi2SKBuilder& operator=( const BPHBuToPsi2SKBuilder& x ) = delete;

  /** Destructor
   */
  ~BPHBuToPsi2SKBuilder() override {}

  /** Operations
   */
  /// set cuts
  void setKPtMin     ( double pt  ) { setTrkPtMin  (  pt ); }
  void setKEtaMax    ( double eta ) { setTrkEtaMax ( eta ); }
  void setPsi2SMassMin( double m   ) { setResMassMin(   m ); }
  void setPsi2SMassMax( double m   ) { setResMassMax(   m ); }

  /// get current cuts
  double getKPtMin     () const { return getTrkPtMin  (); }
  double getKEtaMax    () const { return getTrkEtaMax (); }
  double getPsi2SMassMin() const { return getResMassMin(); }
  double getPsi2SMassMax() const { return getResMassMax(); }

};


#endif

