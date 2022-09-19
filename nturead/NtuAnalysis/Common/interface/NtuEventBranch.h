#ifndef NtuEventBranch_H
#define NtuEventBranch_H

#include "NtuAnalysis/Common/interface/NtuEventHeader.h"
#include <string>
#include <vector>

template <class T>
class NtuEventBranch: public virtual NtuEventHeader, public virtual T {

 public:

  virtual ~NtuEventBranch();

 protected:

  NtuEventBranch() {
    this->setBranch( "runNumber"    , &  runNumber      ,
                     "runNumber/i"  , &b_runNumber       );
    this->setBranch( "lumiSection"  , &  lumiSection    ,
                     "lumiSection/i", &b_lumiSection     );
    this->setBranch( "eventNumber"  , &  eventNumber    ,
                     "eventNumber/i", &b_eventNumber     );
  }

 private:

  // dummy copy constructor and assignment
  NtuEventBranch           ( const NtuEventBranch& td );
  NtuEventBranch& operator=( const NtuEventBranch& td );

};

template <class T>
NtuEventBranch<T>::~NtuEventBranch() {
}

#endif

