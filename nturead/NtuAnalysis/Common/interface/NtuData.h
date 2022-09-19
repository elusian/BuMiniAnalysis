#ifndef NtuData_H
#define NtuData_H

#include "NtuAnalysis/Common/interface/NtuUtil.h"
#include "NtuAnalysis/Common/interface/NtuEventHeader.h"

#include <string>
#include <vector>

template<class T>
class NtuData: public NtuUtil<T>,
               public virtual NtuEventHeader {

 public:

  typedef T number;

 protected:

  NtuData();
  const char* typedName( const std::string& name ) {
    static std::string s;
    s = name + "/" + numberType;
    return s.c_str();
  }
  std::string numberType;

 private:

  // dummy copy constructor and assignment
  NtuData           ( const NtuData& );
  NtuData& operator=( const NtuData& );

};

#endif

