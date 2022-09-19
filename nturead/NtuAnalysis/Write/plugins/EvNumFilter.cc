#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "NtuAnalysis/Write/plugins/EvNumFilter.h"

// system include files
#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>

EvNumFilter::EvNumFilter( const edm::ParameterSet& ps ):
 eventList( "evtlist" ),
 listType ( "" ) {
  if ( ps.exists( "eventList" ) )
                   eventList = ps.getParameter<std::string>( "eventList" );
  if ( ps.exists( "listType"  ) )
                   listType  = ps.getParameter<std::string>( "listType"  );
  skipList = ( read( eventList ) && ( listType == "skip" ) );
}


EvNumFilter::~EvNumFilter() {
}


void EvNumFilter::beginJob() {
  return;
}


void EvNumFilter::endJob()  {
  return;
}


bool EvNumFilter::filter( edm::Event& ev,
                    const edm::EventSetup& es ) {
  unsigned int runNumber   = ev.id().run();
  unsigned int eventNumber = ev.id().event();
  return ( skipList != find( runNumber, eventNumber ) );
}

//define this as a plug-in
DEFINE_FWK_MODULE( EvNumFilter );
