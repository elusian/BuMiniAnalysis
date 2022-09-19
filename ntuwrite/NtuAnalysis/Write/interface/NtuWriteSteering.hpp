
#include "TFile.h"

// system include files
#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>

template <class T>
NtuWriteSteering<T>::NtuWriteSteering( const edm::ParameterSet& ps ):
  NtuWriteInterface<T>( ps, this ) {
  histName = ps.getUntrackedParameter<std::string>( "histName" );
}


template <class T>
NtuWriteSteering<T>::~NtuWriteSteering() {
}


template <class T>
void NtuWriteSteering<T>::beginJob() {
  NtuWriteInterface<T>::beginJob();
  analyzedFile = 0;
  return;
}


template <class T>
void NtuWriteSteering<T>::endJob()  {
  NtuWriteInterface<T>::endJob();
  TreeWrapper::save( histName );
  return;
}


template <class T>
void NtuWriteSteering<T>::beginRun( const edm::Run& run,
                                    const edm::EventSetup& es ) {
  currentRun     = &run;
  currentEvSetup = &es;
  NtuWriteInterface<T>::beginRun();
  return;
}


template <class T>
void NtuWriteSteering<T>::endRun( const edm::Run& run,
                                  const edm::EventSetup& es ) {
  NtuWriteInterface<T>::endRun();
  return;
}


template <class T>
void NtuWriteSteering<T>::analyze( const edm::Event& ev,
                                   const edm::EventSetup& es ) {
  currentEvent   = &ev;
  currentEvSetup = &es;
  int ientry = 0;
  NtuWriteInterface<T>::analyzeEDM( ev, ientry, analyzedFile++ );
  return;
}

