#include "TMPAnalysis/EDM/interface/TMPEDMToNtuple.h"
#include "NtuAnalysis/Write/interface/NtuWriteSteering.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

class TMPNtuplizer: public NtuWriteSteering<TMPEDMToNtuple> {

 public:

  explicit TMPNtuplizer( const edm::ParameterSet& ps );

};


TMPNtuplizer::TMPNtuplizer( const edm::ParameterSet& ps ):
  NtuWriteSteering<TMPEDMToNtuple>( ps ) {
}


//define this as a plug-in
DEFINE_FWK_MODULE( TMPNtuplizer );
