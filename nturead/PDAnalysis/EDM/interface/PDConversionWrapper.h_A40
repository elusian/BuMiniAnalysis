#ifndef PDConversionWrapper_H
#define PDConversionWrapper_H
/** \class PDConversionWrapper
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
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

//---------------
// C++ Headers --
//---------------
#include <iostream>

//              ---------------------
//              -- Class Interface --
//              ---------------------

class PDConversionWrapper {

 public:

  /** Operations
   */
  static const reco::Conversion* matchedConversion(
               const pat::Electron& electron,
               const edm::Handle<reco::ConversionCollection> &convCollection,
               const math::XYZPoint& beamspot,
               bool allowCkfMatch = true,
               float lxyMin = 2.0,
               float probMin = 1.0e-06,
               unsigned int nHitsBeforeVtxMax = 0 ) {
    return ConversionTools::matchedConversion( electron, *convCollection,
                                               beamspot,
                                               allowCkfMatch, lxyMin, probMin,
                                               nHitsBeforeVtxMax );
  }
  
};

//#include "PDAnalysis/EDM/interface/PDConversionWrapper.hpp"

#endif // PDConversionWrapper_H

