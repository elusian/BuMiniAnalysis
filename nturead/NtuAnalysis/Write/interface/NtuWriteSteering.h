
#ifndef NtuWriteSteering_h
#define NtuWriteSteering_h

#include "NtuAnalysis/Write/interface/NtuWriteInterface.h"
#include "NtuAnalysis/Common/interface/NtuAnalyzerUtil.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "NtuAnalysis/Write/interface/NtuEDConsumer.h"
#include "NtuTool/EDM/interface/EDAnalyzerWrapper.h"

#include <string>
#include <iostream>
#include <fstream>

//
// class declaration
//

template <class T>
class NtuWriteSteering: public NtuEDConsumer<EDAnalyzerWrapper>,
                        public NtuWriteInterface<T>,
                        public virtual NtuAnalyzerUtil {

 public:

  explicit NtuWriteSteering( const edm::ParameterSet& ps );
  ~NtuWriteSteering();

  virtual void beginJob();
  virtual void endJob();

 private:

  // dummy copy constructor and assignment
  NtuWriteSteering           ( const NtuWriteSteering& c );
  NtuWriteSteering& operator=( const NtuWriteSteering& c );

  int analyzedFile;
  std::string histName;

  virtual void beginRun( const edm::Run  & run, const edm::EventSetup& es );
  virtual void endRun(   const edm::Run  & run, const edm::EventSetup& es );
  virtual void analyze(  const edm::Event&  ev, const edm::EventSetup& es );

};

#include "NtuAnalysis/Write/interface/NtuWriteSteering.hpp"

#endif
