#ifndef NtuEDMReader_H
#define NtuEDMReader_H

#include "NtuTool/Read/interface/EDMTreeReader.h"
#include <iostream>
#include <string>

template<class T>
class NtuEDMReader: public T, public EDMTreeReader {

 public:

  NtuEDMReader( const std::string& process, const std::string& producer ) {
    std::cout << "NtuEDMReader" << std::endl;
    setProcessName( process );
    setProducerName( producer );
  }

  virtual void beginJob() {
    T::setupNtuple();
    T::beginJob();
    return;
  }

  virtual bool getEntry( int ientry ) {
    return T::getEntry( ientry );
  };

 private:

  // dummy copy constructor and assignment
  NtuEDMReader           ( const NtuEDMReader& );
  NtuEDMReader& operator=( const NtuEDMReader& );

};

#endif
