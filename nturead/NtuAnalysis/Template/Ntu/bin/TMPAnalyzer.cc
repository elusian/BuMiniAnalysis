#include <iostream>
#include <sstream>
#include <string>
#include <math.h>

#include "TMPAnalyzer.h"

#include "TDirectory.h"
#include "TBranch.h"
#include "TTree.h"
#include "TCanvas.h"
#include "Math/LorentzVector.h"

using namespace std;

TMPAnalyzer::TMPAnalyzer() {

  cout << "new TMPAnalyzer" << endl;

  // user parameters are set as names associated to a string, 
  // default values can be set in the analyzer class contructor

  setUserParameter( "verbose", "f" );
  setUserParameter( "ptCut", "40.0" );

}


TMPAnalyzer::~TMPAnalyzer() {
}



void TMPAnalyzer::beginJob() {

  TMPAnalyzerUtil::beginJob();

  // user parameters are retrieved as strings by using their names;
  // numeric parameters ( int, float or whatever ) can be directly set
  // by passing the corresponding variable,
  // e.g. getUserParameter( "name", x )

  getUserParameter( "verbose", verbose );
  getUserParameter( "ptCut", ptCut );

  return;

}


void TMPAnalyzer::book() {

  // putting "autoSavedObject" in front of the histo creation 
  // it's automatically marked for saving on file; the option 
  // is uneffective when not using the full utility

  autoSavedObject =
  hptmumax        = new TH1D( "hptmumax"    , "ptmumax"    ,  50, 0.0, 100.0 );
  autoSavedObject =
  hptmu2nd        = new TH1D( "hptmu2nd"    , "ptmu 2nd"   ,  50, 0.0, 100.0 );
  autoSavedObject =
  hptmu           = new TH1D( "hptmu"       , "ptmu"       ,  50, 0.0, 100.0 );

  return;

}


void TMPAnalyzer::reset() {
#  if UTIL_USE == FULL
  // automatic reset of ntu content at zero
  autoReset();
#elif UTIL_USE == BARE
  // if ntu content is to be reset at zero that should be done manually
  // when not using the full utility
#endif
  return;
}


bool TMPAnalyzer::analyze( int entry, int event_file, int event_tot ) {

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

  // flag to be set "true" or "false" for events accepted or rejected
  bool flag;

  flag = true;

  int iMuon;
  float ptmu;
  float ptmumax = -1.0;
  float ptmu2nd = -1.0;
  for ( iMuon = 0; iMuon < nMuons; ++iMuon ) {
    ptmu = muoPt->at( iMuon );
    hptmu->Fill( ptmu );
    if( ptmu > ptmu2nd ) {
      if( ptmu > ptmumax ) {
        ptmu2nd = ptmumax;
        ptmumax = ptmu;
      }
      else {
        ptmu2nd = ptmu;
      }
    }
  }
  hptmumax->Fill( ptmumax );
  hptmu2nd->Fill( ptmu2nd );	

  flag = ptmumax > ptCut;

  return flag;

}


void TMPAnalyzer::endJob() {
  return;
}


void TMPAnalyzer::save() {
#  if UTIL_USE == FULL
  // explicit saving not necessary for "autoSavedObjects"
  autoSave();
#elif UTIL_USE == BARE
  // explicit save histos when not using the full utility
  hptmumax->Write();
  hptmu2nd->Write();
  hptmu   ->Write();
#endif
  return;
}


void TMPAnalyzer::plot() {
/*
  TCanvas* can = new TCanvas( "muoPt", "muoPt", 800, 600 );
  can->cd();
  can->Divide( 1, 2 );
  can->cd( 1 );
  hptmumax->Draw();
  hptmu2nd->Draw();
*/
  return;
}

