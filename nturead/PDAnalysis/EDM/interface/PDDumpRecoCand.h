#ifndef PDDumpRecoCand_H
#define PDDumpRecoCand_H
/** \class PDDumpRecoCand
 *
 *  Description: 
 *
 *
 *  $Date: 2015-07-24 11:29:20 $
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
#include "FWCore/Framework/interface/Event.h"

class BPHRecoCandidate;
class BPHPlusMinusCandidate;

//---------------
// C++ Headers --
//---------------
#include <vector>

//              ---------------------
//              -- Class Interface --
//              ---------------------

class PDDumpRecoCand {

 public:

  /** Constructor
   */
  PDDumpRecoCand();

  /** Destructor
   */
  virtual ~PDDumpRecoCand();

  /** Operations
   */
  static void dumpRecoCand( const std::string& name,
                            const BPHRecoCandidate* cand );

 private:

  // private copy and assigment constructors
  PDDumpRecoCand           ( const PDDumpRecoCand& x );
  PDDumpRecoCand& operator=( const PDDumpRecoCand& x );

  static float energy ( float px, float py, float pz, float mm );
  static float invMass( float px, float py, float pz, float ee );


};


#endif // PDDumpRecoCand_H

