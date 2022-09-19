#define UTIL_USE FULL
#include "NtuAnalysis/Read/interface/NtuReader.h"
#include "NtuAnalysis/Read/interface/NtuEDMReader.h"

/** \class NtuFlexibleAnalyzer
 *
 *  Description: 
 *    Class to steer the tree analysis:
 *          - open all files in sequence
 *          - loop over the tree for each file
 *
 *  $Date: 2011/07/27 10:32:24 $
 *  $Revision: 1.3 $
 *  \author Paolo Ronchese INFN Padova
 *
 */

//----------------------
// Base Class Headers --
//----------------------
#include "NtuTool/Read/interface/TreeStandardAnalyzer.h"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
#include "NtuTool/Read/interface/TreeReader.h"
#include "TChain.h"

//---------------
// C++ Headers --
//---------------
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
class TreeReader;

//              ---------------------
//              -- Class Interface --
//              ---------------------

template<class T>
class NtuFlexibleAnalyzer: public TreeStandardAnalyzer {

 public:

  /** Constructor
   */
  NtuFlexibleAnalyzer( const std::string& process,
                       const std::string& producer ):
   ntuProdProcess ( process  ),
   ntuProdProducer( producer ) {
  }

  /** Destructor
   */
  virtual ~NtuFlexibleAnalyzer() {}

  /** Operations
   */
  /// run the application
  virtual int run( int argc, char* argv[] ) {

    std::cout << "NtuFlexibleAnalyzer::run" << std::endl;

    char** argp = argv;
    char** argl = argp + argc;
    std::string key = "";
    std::string val = "";
    enum ntuType { ntu, edm };
    ntuType type = ntu;
    while ( argp < argl ) {
      std::string args( *argp++ );
      if ( args == "-v" ) {
        key = *argp++;
        val = *argp++;
        if ( key == "process"  ) ntuProdProcess  = val;
        if ( key == "producer" ) ntuProdProducer = val;
        if ( key == "ntuType"  ) {
          if      ( val == "ntu" ) type = ntu;
          else if ( val == "edm" ) type = edm;
          else                     std::cout << "invalid ntuple type: " << val
                                             << " choose \"ntu\" or \"edm\""
                                             << std::endl;
        }
      }
    }
    
    if ( type == ntu ) new NtuReader<T>;
    if ( type == edm ) new NtuEDMReader<T>( ntuProdProcess, ntuProdProducer );

    return TreeStandardAnalyzer::run( argc, argv );

  }

 private:

  NtuFlexibleAnalyzer( const NtuFlexibleAnalyzer& t );
  NtuFlexibleAnalyzer& operator=( const NtuFlexibleAnalyzer& t );

  std::string ntuProdProcess;
  std::string ntuProdProducer;

};

