#ifndef HeavyFlavorAnalysis_SpecificDecay_BPHPsi2SToJPsiPiPiBuilder_h
#define HeavyFlavorAnalysis_SpecificDecay_BPHPsi2SToJPsiPiPiBuilder_h
/** \class BPHPsi2SToJPsiPiPiBuilder
 *
 *  Description: 
 *     Class to build Psi2S to JPsi pi+ pi- candidates
 *
 *  \author Paolo Ronchese INFN Padova
 *
 */

//----------------------
// Base Class Headers --
//----------------------


//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
#include "HeavyFlavorAnalysis/RecoDecay/interface/BPHRecoBuilder.h"
#include "HeavyFlavorAnalysis/RecoDecay/interface/BPHRecoCandidate.h"
#include "HeavyFlavorAnalysis/RecoDecay/interface/BPHPlusMinusCandidate.h"

#include "FWCore/Framework/interface/Event.h"

class BPHParticleChargeSelect;
class BPHParticlePtSelect;
class BPHParticleEtaSelect;
class BPHMassSelect;
class BPHChi2Select;
class BPHMassFitSelect;

//---------------
// C++ Headers --
//---------------
#include <string>
#include <vector>

//              ---------------------
//              -- Class Interface --
//              ---------------------

class BPHPsi2SToJPsiPiPiBuilder {

 public:

  /** Constructor
   */
  BPHPsi2SToJPsiPiPiBuilder( const edm::EventSetup& es,
      const std::vector<BPHPlusMinusConstCandPtr>& jpsiCollection,
      const BPHRecoBuilder::BPHGenericCollection*   posCollection,
      const BPHRecoBuilder::BPHGenericCollection*   negCollection );

  // deleted copy constructor and assignment operator
  BPHPsi2SToJPsiPiPiBuilder           ( const BPHPsi2SToJPsiPiPiBuilder& x ) = delete;
  BPHPsi2SToJPsiPiPiBuilder& operator=( const BPHPsi2SToJPsiPiPiBuilder& x ) = delete;

  /** Destructor
   */
  virtual ~BPHPsi2SToJPsiPiPiBuilder();

  /** Operations
   */
  /// build Psi2S candidates
  std::vector<BPHRecoConstCandPtr> build();

  /// set cuts
  void setPiPtMin    ( double pt  );
  void setPiEtaMax   ( double eta );
  void setJPsiMassMin( double m   );
  void setJPsiMassMax( double m   );
  void setMassMin    ( double m   );
  void setMassMax    ( double m   );
  void setProbMin    ( double p   );
  void setMassFitMin ( double m   );
  void setMassFitMax ( double m   );
  void setConstr     ( bool flag  );

  /// get current cuts
  double getPiPtMin    () const;
  double getPiEtaMax   () const;
  double getJPsiMassMin() const;
  double getJPsiMassMax() const;
  double getMassMin    () const;
  double getMassMax    () const;
  double getProbMin    () const;
  double getMassFitMin () const;
  double getMassFitMax () const;
  bool   getConstr     () const;

 private:

  std::string    jPsiName;
  std::string pionPosName;
  std::string pionNegName;

  const edm::EventSetup* evSetup;
  const std::vector<BPHPlusMinusConstCandPtr>* jCollection;
  const BPHRecoBuilder::BPHGenericCollection*  pCollection;
  const BPHRecoBuilder::BPHGenericCollection*  nCollection;

  BPHMassSelect          * jpsiSel;
  double ptMin;
  double etaMax;

  BPHMassSelect          * massSel;
  BPHChi2Select          * chi2Sel;
  BPHMassFitSelect       * mFitSel;

  bool massConstr;
  float minPDiff;
  bool updated;

  std::vector<BPHRecoConstCandPtr> psi2SList;

};


#endif

