#include "PDAnalysis/EDM/interface/PDEDMToNtuple.h"
#include "NtuAnalysis/Write/interface/NtuFilter.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

class PDFilter: public NtuFilter<PDEDMToNtuple> {

 public:

  explicit PDFilter( const edm::ParameterSet& ps );

};


PDFilter::PDFilter( const edm::ParameterSet& ps ):
  NtuFilter<PDEDMToNtuple>( ps ) {
}


//define this as a plug-in
DEFINE_FWK_MODULE( PDFilter );
