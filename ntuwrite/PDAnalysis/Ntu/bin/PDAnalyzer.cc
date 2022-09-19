#include <iostream>
#include <sstream>
#include <string>
#include <math.h>

#include "PDAnalyzer.h"

#include "TDirectory.h"
#include "TBranch.h"
#include "TTree.h"
#include "TCanvas.h"
#include "Math/LorentzVector.h"

using namespace std;

PDAnalyzer::PDAnalyzer() {

  std::cout << "new PDAnalyzer" << std::endl;

  // user parameters are set as names associated to a string, 
  // default values can be set in the analyzer class contructor

  setUserParameter( "verbose", "f" );

}


PDAnalyzer::~PDAnalyzer() {
}



void PDAnalyzer::beginJob() {

  PDAnalyzerUtil::beginJob();

  // user parameters are retrieved as strings by using their names;
  // numeric parameters ( int, float or whatever ) can be directly set
  // by passing the corresponding variable,
  // e.g. getUserParameter( "name", x )

  getUserParameter( "verbose", verbose );

  return;

}


void PDAnalyzer::book() {

  // putting "autoSavedObject" in front of the histo creation 
  // it's automatically marked for saving on file; the option 
  // is uneffective when not using the full utility

  autoSavedObject =
  hptmumax        = new TH1D( "hptmumax"    , "ptmumax"    ,  50, 0.0, 100.0 );
  autoSavedObject =
  hptmuall        = new TH1D( "hptmuall"    , "ptmuall"    ,  50, 0.0, 100.0 );

  return;

}


void PDAnalyzer::reset() {
// automatic reset
  autoReset();
  return;
}


bool PDAnalyzer::analyze( int entry, int event_file, int event_tot ) {

  if ( verbose ) {
    cout << " +++++++++++++++++++++++++++ " << endl;
    cout << "entry: "
         << entry << " " << event_file << " " << event_tot << endl;
    cout << "run: " <<   runNumber << " , "
         << "evt: " << eventNumber << endl;
  }
  else {
//    if ( !( event_file % 10000 ) || !( event_tot % 10000 ) )
    if ( !( event_tot % 10000 ) && event_tot )
      cout << event_file << " " << event_tot << endl;
  }

/*
  int iMuon;
  float ptmu;
  float ptmumax = -1.0;
  for ( iMuon = 0; iMuon < nMuons; ++iMuon ) {
    ptmu = muoPt->at( iMuon );
    hptmuall->Fill( ptmu );
    if( ptmu > ptmumax ) {
      ptmumax = ptmu;
    }
  }
  hptmumax->Fill( ptmumax );
*/

  return true;

}


void PDAnalyzer::endJob() {
  return;
}

