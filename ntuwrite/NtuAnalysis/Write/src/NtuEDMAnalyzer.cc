#include "NtuAnalysis/Write/interface/NtuEDMAnalyzer.h"
#include "TFile.h"
#include "TROOT.h"

#include <iostream>
#include <sstream>
#include <map>

using namespace std;

NtuEDMAnalyzer::NtuEDMAnalyzer( const edm::ParameterSet& ps ) {
  parameterSet = &ps;
  if ( ps.exists( "eventList" ) )
  setUserParameter( "eventList", ps.getParameter<string>( "eventList" ) );
  if ( ps.exists( "listType" ) )
  setUserParameter( "listType" , ps.getParameter<string>( "listType"  ) );
}


NtuEDMAnalyzer::~NtuEDMAnalyzer() {
}


void NtuEDMAnalyzer::beginJob() {
  TDirectory* current = gDirectory;
  gROOT->cd();
  book();
  current->cd();
  return;
}


void NtuEDMAnalyzer::beginRun() {
  return;
}


bool NtuEDMAnalyzer::analyzeEDM( const edm::EventBase& ev,
                                 int ientry, int event_file ) {

  reset();

  // store event id
  runNumber   = ev.id().run();
  lumiSection = ev.id().luminosityBlock();
  eventNumber = ev.id().event();

  // event select
  if ( skipList == find( runNumber, eventNumber ) ) return false;
  currentEvBase = &ev;
  read( ev );
  bool select = analyze( ientry, event_file, analyzedEvts++ );
  if ( select ) acceptedEvts++;

  return select;

}


void NtuEDMAnalyzer::read( const edm::EventBase& ev ) {
  return;
}


void NtuEDMAnalyzer::endRun() {
  return;
}


void NtuEDMAnalyzer::endJob() {
  cout << endl << analyzedEvents() << " events read, "
               << acceptedEvents() << " accepted" << endl;
  return;
}

