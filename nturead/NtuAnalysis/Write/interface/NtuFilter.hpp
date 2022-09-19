
#include "TFile.h"

// system include files
#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>

template <class T>
NtuFilter<T>::NtuFilter( const edm::ParameterSet& ps ):
  T( ps, this ) {
  histName = ps.getUntrackedParameter<std::string>( "histName" );
}


template <class T>
NtuFilter<T>::~NtuFilter() {
}


template <class T>
void NtuFilter<T>::beginJob() {
  T::beginJob();
  analyzedFile = 0;
  return;
}


template <class T>
void NtuFilter<T>::endJob()  {
  T::endJob();
  TreeWrapper::save( histName );
  return;
}


template <class T>
void NtuFilter<T>::beginRun( const edm::Run& run,
                             const edm::EventSetup& es ) {
  currentRun     = &run;
  currentEvSetup = &es;
  T::beginRun();
  return;
}


template <class T>
void NtuFilter<T>::endRun( const edm::Run& run,
                           const edm::EventSetup& es ) {
  T::endRun();
  return;
}


template <class T>
bool NtuFilter<T>::filter( edm::Event& ev,
                     const edm::EventSetup& es ) {
  currentEvent   = &ev;
  currentEvSetup = &es;
  int ientry = 0;
  return T::analyzeEDM( ev, ientry, analyzedFile++ );
}

