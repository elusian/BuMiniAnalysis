#include "PDAnalysis/EDM/interface/PDEDMToNtuple.h"
#include "NtuAnalysis/Write/interface/NtuEDMWriteSteering.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

class PDEDMNtuplizer: public NtuEDMWriteSteering<PDEDMToNtuple> {

 public:

  explicit PDEDMNtuplizer( const edm::ParameterSet& ps );

};


PDEDMNtuplizer::PDEDMNtuplizer( const edm::ParameterSet& ps ):
  NtuEDMWriteSteering<PDEDMToNtuple>( ps ) {
}


//define this as a plug-in
DEFINE_FWK_MODULE( PDEDMNtuplizer );
