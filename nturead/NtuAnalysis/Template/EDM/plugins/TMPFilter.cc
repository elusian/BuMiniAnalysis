#include "TMPAnalysis/EDM/interface/TMPEDMToNtuple.h"
#include "NtuAnalysis/Write/interface/NtuFilter.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

class TMPFilter: public NtuFilter<TMPEDMToNtuple> {

 public:

  explicit TMPFilter( const edm::ParameterSet& ps );

};


TMPFilter::TMPFilter( const edm::ParameterSet& ps ):
  NtuFilter<TMPEDMToNtuple>( ps ) {
}


//define this as a plug-in
DEFINE_FWK_MODULE( TMPFilter );
