#ifndef PDNtuParameterSetWrapper_H
#define PDNtuParameterSetWrapper_H
/** \class PDNtuParameterSetWrapper
 *
 *  Description: 
 *
 *
 *  $Date: 2019-05-28 14:55:56 $
 *  $Revision: 1.1 $
 *  \author Paolo Ronchese INFN Padova
 *
 */

//----------------------
// Base Class Headers --
//----------------------


//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

//---------------
// C++ Headers --
//---------------


//              ---------------------
//              -- Class Interface --
//              ---------------------

class NtuParameterSetWrapper {

 public:

  /** Operations
   */
  static std::shared_ptr<edm::ParameterSet> readPSetsFrom(
                                            const std::string& fileOrString ) {
    return edm::readPSetsFrom( fileOrString );
  }

};

//#include "PDAnalysis/EDM/interface/PDNtuParameterSetWrapper.hpp"

#endif // PDNtuParameterSetWrapper_H

