#include <iostream>
#include <sstream>
#include <string>
#include <math.h>

#include "TMPAnalyzerUtil.h"

TMPAnalyzerUtil::TMPAnalyzerUtil() {

  // user parameters are set as names associated to a string, 
  // default values can be set in the analyzer class contructor
  setUserParameter( "use_muons" , "false" );
  setUserParameter( "use_jets"  , "false" );

}


TMPAnalyzerUtil::~TMPAnalyzerUtil() {
}


void TMPAnalyzerUtil::setupNtuple() {

  getUserParameter( "use_muons", use_muons );
  getUserParameter( "use_jets" , use_jets  );
  initTree();

  return;

}


void TMPAnalyzerUtil::beginJob() {

  NtuInterface< TMPNtupleBranch<WrapperBase> >::beginJob();

  // user parameters are retrieved as strings by using their names;
  // numeric parameters ( int, float or whatever ) can be directly set
  // by passing the corresponding variable,
  // e.g. getUserParameter( "name", x )

  return;

}


// pre-selection, with minimal data process 
// before full ntuple entry reading
bool TMPAnalyzerUtil::preSelect( int ientry ) {
  // get number of muons
  getEntry( b_nMuons, ientry );
  process( b_nMuons, ientry );
  std::cout << nMuons  << " muons" << std::endl;
  // skip events with no muons
  if ( !nMuons ) return false;
  currentTree->GetEntry( ientry );
  return true;
}
