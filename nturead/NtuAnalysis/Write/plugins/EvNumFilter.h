
#ifndef EvNumFilter_h
#define EvNumFilter_h

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "NtuAnalysis/Common/interface/NtuEventSelect.h"
#include "NtuTool/EDM/interface/EDFilterWrapper.h"

#include <string>
#include <iostream>
#include <fstream>

//
// class declaration
//

class EvNumFilter: public EDFilterWrapper,
                   public NtuEventSelect {

 public:

  explicit EvNumFilter( const edm::ParameterSet& ps );
  ~EvNumFilter();

  virtual void beginJob();
  virtual void endJob();

 private:

  // dummy copy constructor and assignment
  EvNumFilter           ( const EvNumFilter& c );
  EvNumFilter& operator=( const EvNumFilter& c );

  std::string eventList;
  std::string listType;
  bool skipList;

  virtual bool filter(         edm::Event& ev, const edm::EventSetup& es );

};

#endif
