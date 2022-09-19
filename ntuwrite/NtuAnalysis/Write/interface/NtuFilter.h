
#ifndef NtuFilter_h
#define NtuFilter_h

#include "NtuAnalysis/Common/interface/NtuAnalyzerUtil.h"

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "NtuAnalysis/Write/interface/NtuEDConsumer.h"
#include "NtuTool/EDM/interface/EDFilterWrapper.h"

#include <string>
#include <iostream>
#include <fstream>

//
// class declaration
//

template <class T>
class NtuFilter: public NtuEDConsumer<EDFilterWrapper>,
                 public T,
                 public virtual NtuAnalyzerUtil {

 public:

  explicit NtuFilter( const edm::ParameterSet& ps );
  ~NtuFilter();

  virtual void beginJob();
  virtual void endJob();

 private:

  // dummy copy constructor and assignment
  NtuFilter           ( const NtuFilter& c );
  NtuFilter& operator=( const NtuFilter& c );

  int analyzedFile;
  std::string histName;

  virtual void beginRun( const edm::Run& run , const edm::EventSetup& es );
  virtual void endRun(   const edm::Run& run , const edm::EventSetup& es );
  virtual bool filter(         edm::Event& ev, const edm::EventSetup& es );

};

#include "NtuAnalysis/Write/interface/NtuFilter.hpp"

#endif
