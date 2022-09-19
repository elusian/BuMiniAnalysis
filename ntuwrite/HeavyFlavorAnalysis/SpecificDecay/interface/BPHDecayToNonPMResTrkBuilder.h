#ifndef HeavyFlavorAnalysis_SpecificDecay_BPHDecayToNonPMResTrkBuilder_h
#define HeavyFlavorAnalysis_SpecificDecay_BPHDecayToNonPMResTrkBuilder_h
/** \class BPHDecayToNonPMResTrkBuilder
 *
 *  Description: 
 *     Class to build a particle decaying to a resonance, decaying itself
 *     to an opposite charged particles pair, and an additional track
 *
 *  \author Paolo Ronchese INFN Padova
 *
 */

//----------------------
// Base Class Headers --
//----------------------
#include "HeavyFlavorAnalysis/SpecificDecay/interface/BPHDecayNonPMConstrainedBuilder.h"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
#include "HeavyFlavorAnalysis/SpecificDecay/interface/BPHParticlePtSelect.h"
#include "HeavyFlavorAnalysis/SpecificDecay/interface/BPHParticleEtaSelect.h"

#include "HeavyFlavorAnalysis/RecoDecay/interface/BPHRecoBuilder.h"
#include "HeavyFlavorAnalysis/RecoDecay/interface/BPHRecoCandidate.h"
#include "HeavyFlavorAnalysis/RecoDecay/interface/BPHPlusMinusCandidate.h"

#include "FWCore/Framework/interface/Event.h"

class BPHParticleNeutralVeto;

//---------------
// C++ Headers --
//---------------
#include <string>
#include <vector>

//              ---------------------
//              -- Class Interface --
//              ---------------------

class BPHDecayToNonPMResTrkBuilder: public BPHDecayNonPMConstrainedBuilder {

 public:

  /** Constructor
   */
  BPHDecayToNonPMResTrkBuilder( const edm::EventSetup& es,
      const std::string& resName, double resMass, double resWidth,
      const std::vector<BPHRecoConstCandPtr>& resCollection,
      const std::string& trkName, double trkMass, double trkSigma,
      const BPHRecoBuilder::BPHGenericCollection*  trkCollection );

  // deleted copy constructor and assignment operator
  BPHDecayToNonPMResTrkBuilder           ( const BPHDecayToNonPMResTrkBuilder& x )
                                      = delete;
  BPHDecayToNonPMResTrkBuilder& operator=( const BPHDecayToNonPMResTrkBuilder& x )
                                      = delete;

  /** Destructor
   */
  ~BPHDecayToNonPMResTrkBuilder() override;

  /** Operations
   */
  /// build candidates
  std::vector<BPHRecoConstCandPtr> build();

  /// set cuts
  void setTrkPtMin  ( double pt  );
  void setTrkEtaMax ( double eta );

  /// get current cuts
  double getTrkPtMin  () const { return  ptSel->getPtMin (); }
  double getTrkEtaMax () const { return etaSel->getEtaMax(); }

 private:

  std::string tName;
  double tMass;
  double tSigma;

  const BPHRecoBuilder::BPHGenericCollection*  tCollection;

  BPHParticleNeutralVeto* tknVeto;
  BPHParticlePtSelect   *   ptSel;
  BPHParticleEtaSelect  *  etaSel;

  std::vector<BPHRecoConstCandPtr> recList;

};


#endif

