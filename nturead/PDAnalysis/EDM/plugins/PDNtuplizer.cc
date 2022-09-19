#include "PDAnalysis/EDM/interface/PDEDMToNtuple.h"
#include "NtuAnalysis/Write/interface/NtuWriteSteering.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

class PDNtuplizer: public NtuWriteSteering<PDEDMToNtuple> {

 public:

  explicit PDNtuplizer( const edm::ParameterSet& ps );

};


PDNtuplizer::PDNtuplizer( const edm::ParameterSet& ps ):
  NtuWriteSteering<PDEDMToNtuple>( ps ) {
}


//define this as a plug-in
DEFINE_FWK_MODULE( PDNtuplizer );
