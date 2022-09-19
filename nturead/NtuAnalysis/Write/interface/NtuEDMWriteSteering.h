
#ifndef NtuEDMWriteSteering_h
#define NtuEDMWriteSteering_h

#include "NtuAnalysis/Common/interface/NtuAnalyzerUtil.h"

#include "NtuTool/EDM/interface/EDMTreeWriter.h"

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "NtuAnalysis/Write/interface/NtuEDConsumer.h"
#include "NtuTool/EDM/interface/EDProducerWrapper.h"

#include <string>
#include <iostream>
#include <fstream>

//
// class declaration
//

template <class T>
class NtuEDMWriteSteering: public EDMTreeWriter,
                           public NtuEDConsumer<EDProducerWrapper>,
                           public T,
                           public virtual NtuAnalyzerUtil {

 public:

  explicit NtuEDMWriteSteering( const edm::ParameterSet& ps );
  ~NtuEDMWriteSteering();

  virtual void beginJob();
  virtual void endJob();

 private:

  // dummy copy constructor and assignment
  NtuEDMWriteSteering           ( const NtuEDMWriteSteering& c );
  NtuEDMWriteSteering& operator=( const NtuEDMWriteSteering& c );

  int analyzedFile;
  std::string histName;

  virtual void beginRun(   edm::Run  & run, const edm::EventSetup& es );
  virtual void endRun(     edm::Run  & run, const edm::EventSetup& es );
  virtual bool fill( const edm::Event&  ev, const edm::EventSetup& es );

};

#include "NtuAnalysis/Write/interface/NtuEDMWriteSteering.hpp"

#endif
