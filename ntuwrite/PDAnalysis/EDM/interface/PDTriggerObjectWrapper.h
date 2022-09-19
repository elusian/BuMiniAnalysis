#ifndef PDTriggerObjectWrapper_H
#define PDTriggerObjectWrapper_H
/** \class PDTriggerObjectWrapper
 *
 *  Description: 
 *
 *
 *  $Date: 2015-12-30 11:15:56 $
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


//---------------
// C++ Headers --
//---------------
#include <iostream>

//              ---------------------
//              -- Class Interface --
//              ---------------------

class PDTriggerObjectWrapper {

 public:

  /** Operations
   */
  /// get number of hits
  template<class O, class E, class R>
  static void unpackFilterLabels( O& to, const E& eb, const R& tr ) {
    to.unpackFilterLabels( eb, tr );
    return;
  }

};

#endif // PDTriggerObjectWrapper_H

