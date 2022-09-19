
// system include files
#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>

template <class T>
NtuEDMWriteSteering<T>::NtuEDMWriteSteering( const edm::ParameterSet& ps ):
  T( ps, this ) {
  histName = ps.getUntrackedParameter<std::string>( "histName" );
  initWrite();
}


template <class T>
NtuEDMWriteSteering<T>::~NtuEDMWriteSteering() {
}


template <class T>
void NtuEDMWriteSteering<T>::beginJob() {
  T::beginJob();
  analyzedFile = 0;
  return;
}


template <class T>
void NtuEDMWriteSteering<T>::endJob()  {
  T::endJob();
  TreeWrapper::save( histName );
  return;
}


template <class T>
void NtuEDMWriteSteering<T>::beginRun( edm::Run& run,
                                 const edm::EventSetup& es ) {
  currentRun     = &run;
  currentEvSetup = &es;
  T::beginRun();
  return;
}


template <class T>
void NtuEDMWriteSteering<T>::endRun( edm::Run& run,
                               const edm::EventSetup& es ) {
  T::endRun();
  return;
}


template <class T>
bool NtuEDMWriteSteering<T>::fill( const edm::Event& ev,
                                   const edm::EventSetup& es ) {
  currentEvent   = &ev;
  currentEvSetup = &es;
  int ientry = 0;
  return T::analyzeEDM( ev, ientry, analyzedFile++ );
}

