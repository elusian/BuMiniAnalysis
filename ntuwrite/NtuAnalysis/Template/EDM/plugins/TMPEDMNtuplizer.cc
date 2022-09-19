#include "TMPAnalysis/EDM/interface/TMPEDMToNtuple.h"
#include "NtuAnalysis/Write/interface/NtuEDMWriteSteering.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

class TMPEDMNtuplizer: public NtuEDMWriteSteering<TMPEDMToNtuple> {

 public:

  explicit TMPEDMNtuplizer( const edm::ParameterSet& ps );

};


TMPEDMNtuplizer::TMPEDMNtuplizer( const edm::ParameterSet& ps ):
  NtuEDMWriteSteering<TMPEDMToNtuple>( ps ) {
}


//define this as a plug-in
DEFINE_FWK_MODULE( TMPEDMNtuplizer );
