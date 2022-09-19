#include "TMPAnalysis/EDM/interface/TMPEDMToNtuple.h"
#include "TFile.h"
#include "TROOT.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include <iostream>
#include <sstream>
#include <map>

using pat::Muon;
using pat::Jet;

using namespace std;

void TMPEDMToNtuple::build( const edm::ParameterSet& ps ) {

  cout << "TMPEDMToNtuple::TMPEDMToNtuple" << endl;

  GET_UNTRACKED( verbose, "false" );

  // get label of muons collection and
  // switch on muon blocklet if not null
//  getParameter( "labelMuons", labelMuons, "" );
  GET_PARAMETER( labelMuons, "" );
  use_muons = ( labelMuons != "" );
  setUserParameter( "use_muons", use_muons ? "t" : "f" );

  // get label of jets collection and
  // switch on muon blocklet if not null
  GET_PARAMETER( labelJets , "" );
  use_jets = ( labelJets != "" );
  setUserParameter( "use_jets", use_jets ? "t" : "f" );

  GET_PARAMETER( ptCut, "10.0" );

  setupNtuple();

}


TMPEDMToNtuple::~TMPEDMToNtuple() {
}


void TMPEDMToNtuple::beginJob() {

  TMPAnalyzer::beginJob();
  NtuEDMAnalyzer::beginJob();

  return;

}


void TMPEDMToNtuple::beginRun() {
  return;
}


void TMPEDMToNtuple::read( const edm::EventBase& ev ) {

  nMuons      = 0;
  nJets       = 0;
  if ( use_muons ) fillMuons();
  if ( use_jets  ) fillJets ();

  return;

}


void TMPEDMToNtuple::endRun() {
  return;
}


void TMPEDMToNtuple::endJob() {
  TMPAnalyzer::endJob();
  NtuEDMAnalyzer::endJob();
  return;
}


void TMPEDMToNtuple::fillMuons() {

  if ( labelMuons == "RANDOM" ) {
    float pt;
    float pm = 1.0e+30;
    while ( ( pt = -50.0 * log ( random() * 1.0 / RAND_MAX ) ) < pm ) {
      if ( ( pt < 5.0 ) && !nMuons ) break;
      ++nMuons;
      muoPt->push_back( pm = pt );
    }
    return;
  }

  // get muons through an interface to access data by label or by token
  // according to CMSSW version
  gt_muons.get( currentEvBase, muons );
  bool vMuons = muons.isValid();

  // store muons info

  int iObj;
  int nObj = ( vMuons ? muons->size() : 0 );
  muoPt->resize( nObj );
  if ( !vMuons ) {
    cout << "invalid muons: " << getUserParameter( "labelMuons" ) << endl;
    return;
  }

  vector<const Muon*> muonPtr;
  muonPtr.resize( nObj );
  for ( iObj = 0; iObj < nObj; ++iObj ) muonPtr[iObj] = &( muons->at( iObj ) );

  CompareByPt<Muon> muoComp;
  sort( muonPtr.begin(), muonPtr.end(), muoComp );

  nMuons = nObj;
  muoPt->resize( nMuons );
  for ( iObj = 0; iObj < nObj; ++iObj ) {
    const Muon& muon = *muonPtr[iObj];
    muoPt->at( iObj ) = muon.pt();
  }

  return;

}


void TMPEDMToNtuple::fillJets() {

  // get muons through an interface to access data by label or by token
  // according to CMSSW version
  gt_jets.get( currentEvBase, jets );
  bool vJets = jets.isValid();

  // store jets info

  int iObj;
  int nObj = ( vJets ? jets->size() : 0 );
  jetPt  ->resize( 0 );
  if ( !vJets ) {
    cout << "invalid jets" << getUserParameter( "labelJets" ) << endl;
    return;
  }

  vector<const Jet*> jetPtr;
  jetPtr.resize( nObj );
  for ( iObj = 0; iObj < nObj; ++iObj ) jetPtr[iObj] = &( jets->at( iObj ) );

  CompareByPt<Jet> jetComp;
  sort( jetPtr.begin(), jetPtr.end(), jetComp );

  nJets = nObj;
  jetPt->resize( nJets );
  for ( iObj = 0; iObj < nObj; ++iObj ) {
    const Jet& jet = *jetPtr[iObj];
    jetPt->at( iObj ) = jet.pt();
  }

  return;

}


