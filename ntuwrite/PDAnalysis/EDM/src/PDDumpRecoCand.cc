/*
 *  See header file for a description of this class.
 *
 *  $Date: 2015-07-24 11:29:20 $
 *  $Revision: 1.1 $
 *  \author Paolo Ronchese INFN Padova
 *
 */

//-----------------------
// This Class' Header --
//-----------------------
#include "PDAnalysis/EDM/interface/PDDumpRecoCand.h"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "HeavyFlavorAnalysis/RecoDecay/interface/BPHRecoBuilder.h"
#include "HeavyFlavorAnalysis/RecoDecay/interface/BPHRecoSelect.h"
#include "HeavyFlavorAnalysis/RecoDecay/interface/BPHPlusMinusCandidate.h"
#include "HeavyFlavorAnalysis/RecoDecay/interface/BPHMomentumSelect.h"
#include "HeavyFlavorAnalysis/RecoDecay/interface/BPHVertexSelect.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include <TMath.h>


//---------------
// C++ Headers --
//---------------
#include <iostream>

using namespace std;

//-------------------
// Initializations --
//-------------------


//----------------
// Constructors --
//----------------
PDDumpRecoCand::PDDumpRecoCand() {
}

//--------------
// Destructor --
//--------------
PDDumpRecoCand::~PDDumpRecoCand() {
}

//--------------
// Operations --
//--------------
void PDDumpRecoCand::dumpRecoCand( const std::string& name,
                                   const BPHRecoCandidate* cand ) {

  static string cType = " cowboy";
  static string sType = " sailor";
  static string dType = "";
  string* type;
  const BPHPlusMinusCandidate* pmCand =
        dynamic_cast<const BPHPlusMinusCandidate*>( cand );
  if ( pmCand != 0 ) {
    if ( pmCand->isCowboy() ) type = &cType;
    else                      type = &sType;
  }
  else                        type = &dType;

  bool constrMass = ( cand->constrMass() > 0.0 );

  cout << "****** " << name << "   cand mass: "
       << cand->composite().mass() << " momentum "
       << cand->composite().px() << " "
       << cand->composite().py() << " "
       << cand->composite().pz() << *type << endl;

  const reco::Vertex& vx = cand->vertex();
  cout << "****** " << name << " vertex mass: "
       << vx.p4().M() << " momentum "
       << vx.p4().px() << " "
       << vx.p4().py() << " "
       << vx.p4().pz() << endl;


  cout << "****** " << name << " constr mass: "
       << cand->p4().mass() << " momentum "
       << cand->p4().px() << " "
       << cand->p4().py() << " "
       << cand->p4().pz() << endl;

  const reco::Vertex::Point& vp = vx.position();
  double chi2 = vx.chi2();
  int    ndof = lround( vx.ndof() );
  double prob = TMath::Prob( chi2, ndof );
  string tdca = "";
  if ( pmCand != 0 ) {
    stringstream sstr;
    sstr << " - " << pmCand->cAppInRPhi().distance();
    tdca = sstr.str();
  }
  cout << "****** " << name << " vertex: "
       << vx.isFake() << " " << vx.isValid() << " - "
       << chi2 << " " << ndof << " " << prob << " - "
       << vp.X() << " " << vp.Y() << " " << vp.Z() << tdca << endl;

  const vector<string>& dl = cand->daugNames();
  int i;
  int n = dl.size();
  float dMass;
  float totPx = 0.0;
  float totPy = 0.0;
  float totPz = 0.0;
  float totEn = 0.0;
  float tvtPx = 0.0;
  float tvtPy = 0.0;
  float tvtPz = 0.0;
  float tvtEn = 0.0;
  for ( i = 0; i < n; ++i ) {
    const string& name = dl[i];
    const reco::Candidate* dp = cand->getDaug( name );
    reco::TransientTrack* tt = cand->getTransientTrack( dp );
    GlobalPoint gp( vp.X(), vp.Y(), vp.Z() ); 
    TrajectoryStateClosestToPoint tscp =
                                  tt->trajectoryStateClosestToPoint( gp );
    GlobalVector gv = tscp.momentum();
//    TrajectoryStateOnSurface tsos = tt->stateOnSurface( gp );
//    GlobalVector gv = tsos.globalMomentum();
    dMass = dp->mass();
    cout << "daughter " << i
         << " " << name << " " << dMass
         << " momentum: "
         << dp->px() << " "
         << dp->py() << " "
         << dp->pz() << " - at vertex: "
         << gv.x() << " "
         << gv.y() << " "
         << gv.z() << endl;
    totPx += dp->px();
    totPy += dp->py();
    totPz += dp->pz();
    totEn += energy( dp->px(), dp->py(), dp->pz(), dMass );
    tvtPx += gv.x();
    tvtPy += gv.y();
    tvtPz += gv.z();
    tvtEn += energy( gv.x(), gv.y(), gv.z(), dMass );
  }
  const vector<string>& dc = cand->compNames();
  int j;
  int m = dc.size();
  for ( j = 0; j < m; ++j ) {
    const string& name = dc[j];
    const BPHRecoCandidate* dp = cand->getComp( name ).get();
    cout << "composite daughter " << j
         << " " << name
         << " momentum: "
         << dp->composite().px() << " "
         << dp->composite().py() << " "
         << dp->composite().pz() << endl;
  }
  if ( !m ) {
    float totMass = invMass( totPx, totPy, totPz, totEn );
    float tvtMass = invMass( tvtPx, tvtPy, tvtPz, tvtEn );
    cout << "total momentum: "
         << totPx << " "
         << totPy << " "
         << totPz << " - at vertex: "
         << tvtPx << " "
         << tvtPy << " "
         << tvtPz << endl;
    cout << "invariant mass: "
         << totMass << " - at vertex: "
         << tvtMass << endl;
  }

  if ( constrMass ) {
  const RefCountedKinematicTree& kt = cand->kinematicTree();
  const RefCountedKinematicVertex& kd = kt->currentDecayVertex();
  GlobalPoint gp = kd->position(); 
  cout << "   kin fit vertex: "
       << gp.x() << " "
       << gp.y() << " "
       << gp.z() << endl;
  const KinematicState& ks = kt->currentParticle()->currentState();
  GlobalVector gv = ks.globalMomentum();
  cout << "   kin fit momentum: "
       << ks.mass() << " - "
       << gv.x() << " "
       << gv.y() << " "
       << gv.z() << " - deltaM: "
       << ks.mass() - cand->constrMass() << endl;
  vector<RefCountedKinematicParticle> dk = kt->finalStateParticles();
  int j;
  int m = dk.size();
  for ( j = 0; j < m; ++j ) {
    const reco::TransientTrack& tt = dk[j]->refittedTransientTrack();
    TrajectoryStateOnSurface tsos = tt.stateOnSurface( gp );
    GlobalVector dm = tsos.globalMomentum();
    cout << "daughter " << j << " refitted: "
         << dm.x() << " "
         << dm.y() << " "
         << dm.z() << endl;
  }
  }

  return;

}


float PDDumpRecoCand::energy ( float px, float py, float pz, float mm ) {
  return sqrt( ( px * px ) + ( py * py ) + ( pz * pz ) + ( mm * mm ) );
}


float PDDumpRecoCand::invMass( float px, float py, float pz, float ee ) {
  float mq = ( px * px ) + ( py * py ) + ( pz * pz ) - ( ee * ee );
  if ( mq > 0.0 ) return 0.0;
  return sqrt( -mq );
}

