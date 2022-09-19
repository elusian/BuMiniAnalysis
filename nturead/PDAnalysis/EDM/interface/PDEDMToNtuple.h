#ifndef PDEDMToNtuple_h
#define PDEDMToNtuple_h

#include "NtuAnalysis/Write/interface/NtuEDMAnalyzer.h"
#include "PDAnalysis/Ntu/bin/PDAnalyzer.h"
#include "PDAnalysis/Ntu/interface/PDEnumString.h"
#include "NtuAnalysis/Write/interface/NtuEDConsumer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"

#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "HeavyFlavorAnalysis/RecoDecay/interface/BPHVertexCompositePtrCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "HeavyFlavorAnalysis/RecoDecay/interface/BPHTrackReference.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "HeavyFlavorAnalysis/RecoDecay/interface/BPHTrackReference.h"
#include "PDAnalysis/EDM/interface/PDVertexWrapper.h"

#include <string>
#include <vector>
#include <set>
#include <iostream>

class BPHRecoCandidate;
class VertexDistanceXY;

class TFile;

namespace edm {
  class LumiReWeighting;
}

class HLTConfigProvider;
class PileupSummaryInfo;

class PDEDMToNtuple: public PDAnalyzer,
                      public NtuEDMAnalyzer {

 public:

  template<class T>
  PDEDMToNtuple( const edm::ParameterSet& ps, NtuEDConsumer<T>* c ):
   NtuEDMAnalyzer( ps ) {
    build( ps );
    ObjectConsumer<T>& oc = getConsumer( c );
    oc.template consume< edm::TriggerResults                  >(
                           gt_trigResults  ,
                         labelTrigResults   );
    oc.template consume< pat::TriggerEvent                    >(
                           gt_trigEvent    ,
                         labelTrigEvent     );
    oc.template consume< std::vector<pat::TriggerObjectStandAlone> >(
                           gt_trigObjects,
                         labelTrigObjects );
    oc.template consume< reco::BeamSpot                       >(
                           gt_beamSpot     ,
                         labelBeamSpot      );
    oc.template consume< std::vector<pat::MET               > >(
                           gt_mets         ,
                         labelMets          );
    oc.template consume< std::vector<pat::Muon              > >(
                           gt_muons        ,
                         labelMuons         );
    oc.template consume< std::vector<pat::Electron          > >(
                           gt_electrons    ,
                         labelElectrons     );
    oc.template consume< reco::ConversionCollection           >(
                           gt_conversions  ,
                         labelConversions   );
    oc.template consume< std::vector<pat::Tau               > >(
                           gt_taus         ,
                         labelTaus          );
    oc.template consume< std::vector<pat::Jet               > >(
                           gt_jets         ,
                         labelJets          );
    oc.template consume< std::vector<reco::PFCandidate      > >(
                           gt_pfCandidates ,
                         labelPFCandidates  );
    oc.template consume< std::vector<pCandidate             > >(
                           gt_pcCandidates ,
                         labelPCCandidates  );
    oc.template consume< std::vector<pat::GenericParticle   > >(
                           gt_gpCandidates ,
                         labelGPCandidates  );
    oc.template consume< std::vector<pat::CompositeCandidate> >(
                           gt_mmCandidates ,
                         labelMMCandidates  );
    oc.template consume< std::vector<pat::CompositeCandidate> >(
                           gt_mmCandidates ,
                         labelMMCandidates  );
    oc.template consume< std::vector<reco::Track            > >(
                           gt_generalTracks,
                         labelGeneralTracks );
    oc.template consume< std::vector<reco::Vertex           > >(
                           gt_pVertices    ,
                         labelPVertices     );
    oc.template consume< std::vector<svtInfo                > >(
                           gt_sVertices    ,
                         labelSVertices     );
    oc.template consume< std::vector<PileupSummaryInfo      > >(
                           gt_PUInfo       ,
                         labelPUInfo        );
    oc.template consume< std::vector<reco::GenParticle      > >(
                           gt_genParticles ,
                         labelGen           );
    oc.template consume< std::vector<reco::GenJet           > >(
                           gt_genJets      ,
                         labelGPJ           );
    std::map<int,std::string>::const_iterator ei_iter = labelMapEvtInfo.begin();
    std::map<int,std::string>::const_iterator ei_iend = labelMapEvtInfo.end();
    while ( ei_iter != ei_iend ) {
      const std::map<int,std::string>::value_type& entry = *ei_iter++;
      EvtInfoHT* ht = evtIht[entry.first] = new EvtInfoHT;
      oc.template consume< double >(
                           ht->gt_evtInfo, entry.second );
    }
    std::map<int,std::string>::const_iterator v0_iter = labelMapV0.begin();
    std::map<int,std::string>::const_iterator v0_iend = labelMapV0.end();
    while ( v0_iter != v0_iend ) {
      const std::map<int,std::string>::value_type& entry = *v0_iter++;
      V0HT* ht = v0ht[entry.first] = new V0HT;
      oc.template consume< std::vector<reco::VertexCompositeCandidate> >(
                           ht->gt_v0, entry.second );
    }
    std::map<int,std::string>::const_iterator vs_iter = labelMapVS.begin();
    std::map<int,std::string>::const_iterator vs_iend = labelMapVS.end();
    while ( vs_iter != vs_iend ) {
      const std::map<int,std::string>::value_type& entry = *vs_iter++;
      VSHT* ht = vsht[entry.first] = new VSHT;
      oc.template consume< std::vector<reco::VertexCompositePtrCandidate> >(
                           ht->gt_v0, entry.second );
    }
  }

  virtual ~PDEDMToNtuple();

  virtual void beginJob();
  virtual void beginRun();
  virtual void read( const edm::EventBase& ev );
  virtual void endRun();
  virtual void endJob();

 private:

  template<class PATObject> class compareByPt {
   public:
    bool operator()( const PATObject* l, const PATObject* r ) {
      return l->pt() > r->pt();
    }
  };

  bool verbose;
  bool evStamp;

  std::string labelTrigResults;
  std::string labelTrigEvent;
  std::string labelTrigObjects;
  std::string labelBeamSpot;
  std::string labelMets;
  std::string labelMuons;
  std::string labelElectrons;
  std::string labelConversions;
  std::string labelTaus;
  std::string labelJets;
  std::string labelPFCandidates;
  std::string labelPCCandidates;
  std::string labelGPCandidates;
  std::string labelMMCandidates;
  std::string labelGeneralTracks;
  std::string labelPVertices;
  std::string labelSVertices;
  std::string labelSVTagInfo;
  std::string labelPUInfo;
  std::string labelGen;
  std::string labelGPJ;
//  std::vector<edm::ParameterSet> eventInfo;

  std::string labelCSV;
  std::string labelTCHE;

//  std::vector< std::string> labelTags;
//  std::vector<edm::ParameterSet> labelTags;
  std::map<int,std::string> labelMapTags;
  std::map<int,std::string> labelMapEvtInfo;
  std::map<int,std::string> labelMapV0;
  std::map<int,std::string> labelMapVS;

  bool selectAssociatedPF;
  bool selectAssociatedTk;
  bool recoverMuonTracks;
  bool recoverCandTracks;
  bool writeAllPrimaryVertices;

  bool read_hlts;
  bool read_hlto;
  bool read_hltm;
  bool read_bspot;
  bool read_met;
  bool read_muons;
  bool read_electrons;
  bool read_taus;
  bool read_jets;
  bool read_pflow;
  bool read_mm;
  bool read_gentrk;
  bool read_tracks;
  bool read_pvts;
  bool read_svts;
  bool read_svtj;
  bool reco_svts;
  bool read_puwgt;
  bool read_gen;
  bool read_gpj;

  bool readTE;
  bool readTO;
  bool readPF;
  bool readPC;
  bool readGP;

  bool fill_svts;
  bool fill_svtt;

  typedef BPHTrackReference::candidate pCandidate;

  typedef PDVertexWrapBase::svtInfo   svtInfo;
  typedef PDVertexWrapBase::sVtx      sVtx;
  typedef PDVertexWrapBase::svtTkList svtTkList;

  typedef PDVertexWrapBase::VtxComponent VtxComponent;

  struct VtxInfo {
    int jetId;
    PDVertexContainer vtx;
    const GlobalVector* dir;
    Measurement1D  d2d;
    Measurement1D  d3d;
    std::vector<VtxComponent> compList;
  };

  std::map< const pat::Jet*,std::vector<VtxInfo> > svtMap;

  edm::Handle< edm::TriggerResults                  >    trigResults;
  NtuEDToken < edm::TriggerResults                  > gt_trigResults;
  edm::Handle< pat::TriggerEvent                    >    trigEvent;
  NtuEDToken < pat::TriggerEvent                    > gt_trigEvent;
  edm::Handle< std::vector<pat::TriggerObjectStandAlone> >    trigObjects;
  NtuEDToken < std::vector<pat::TriggerObjectStandAlone> > gt_trigObjects;
  edm::Handle< reco::BeamSpot                       >    beamSpot;
  NtuEDToken < reco::BeamSpot                       > gt_beamSpot;
  edm::Handle< std::vector<pat::MET               > >    mets;
  NtuEDToken < std::vector<pat::MET               > > gt_mets;
  edm::Handle< std::vector<pat::Muon              > >    muons;
  NtuEDToken < std::vector<pat::Muon              > > gt_muons;
  edm::Handle< std::vector<pat::Electron          > >    electrons;
  NtuEDToken < std::vector<pat::Electron          > > gt_electrons;
  edm::Handle< reco::ConversionCollection           >    conversions;
  NtuEDToken < reco::ConversionCollection           > gt_conversions;
  edm::Handle< std::vector<pat::Tau               > >    taus;
  NtuEDToken < std::vector<pat::Tau               > > gt_taus;
  edm::Handle< std::vector<pat::Jet               > >    jets;
  NtuEDToken < std::vector<pat::Jet               > > gt_jets;
  edm::Handle< std::vector<reco::PFCandidate      > >    pfCandidates;
  NtuEDToken < std::vector<reco::PFCandidate      > > gt_pfCandidates;
  edm::Handle< std::vector<pCandidate             > >    pcCandidates;
  NtuEDToken < std::vector<pCandidate             > > gt_pcCandidates;
  edm::Handle< std::vector<pat::GenericParticle   > >    gpCandidates;
  NtuEDToken < std::vector<pat::GenericParticle   > > gt_gpCandidates;
  edm::Handle< std::vector<pat::CompositeCandidate> >    mmCandidates;
  NtuEDToken < std::vector<pat::CompositeCandidate> > gt_mmCandidates;
  edm::Handle< std::vector<reco::Track            > >    generalTracks;
  NtuEDToken < std::vector<reco::Track            > > gt_generalTracks;
  edm::Handle< std::vector<reco::Vertex           > >    pVertices;
  NtuEDToken < std::vector<reco::Vertex           > > gt_pVertices;
  edm::Handle< std::vector<svtInfo                > >    sVertices;
  NtuEDToken < std::vector<svtInfo                > > gt_sVertices;
  edm::Handle< std::vector<PileupSummaryInfo      > >    PUInfo;
  NtuEDToken < std::vector<PileupSummaryInfo      > > gt_PUInfo;
  edm::Handle< std::vector<reco::GenParticle      > >    genParticles;
  NtuEDToken < std::vector<reco::GenParticle      > > gt_genParticles;
  edm::Handle< std::vector<reco::GenJet           > >    genJets;
  NtuEDToken < std::vector<reco::GenJet           > > gt_genJets;

  struct EvtInfoHT {
    edm::Handle<double>    evtInfo;
    NtuEDToken <double> gt_evtInfo;
  };
  std::map<int,EvtInfoHT*> evtIht;
  struct V0HT {
    edm::Handle< std::vector<reco::VertexCompositeCandidate> >    v0;
    NtuEDToken < std::vector<reco::VertexCompositeCandidate> > gt_v0;
  };
  std::map<int,V0HT*> v0ht;
  struct VSHT {
    edm::Handle< std::vector<reco::VertexCompositePtrCandidate> >    v0;
    NtuEDToken < std::vector<reco::VertexCompositePtrCandidate> > gt_v0;
  };
  std::map<int,VSHT*> vsht;

  const edm::TriggerNames* triggerNames;
  std::vector<std::string> savedTriggerPaths;
  std::vector<std::string> savedTriggerObjects;
  std::vector<std::string> savedJetInfo;
  std::vector<std::string> savedEleInfo;
  HLTConfigProvider* hltConfigProvider;

  std::vector<double> puWgt_mc;
  std::vector<double> puWgt_data;
  edm::LumiReWeighting* LumiWeights;

  std::vector<bool> vRecoFlag;

  bool fixedTrigPathList;
  int lastTrigMapSize;

  float jetPtMin;
  float jetEtaMax;
  float trkDzMax;
  float trkPtMin;
  float trkEtaMax;

  float dRmatchHLT;
  float dPmatchHLT;

  float pvRefitPtMin;

  std::map<int,int> hltPathMap;
  std::map<int,int> hltObjMap;
  std::     map<const pat ::Muon       *,int                   > muoMap;
  std::     map<const pat ::Electron   *,int                   > eleMap;
  std::     map<const pat ::Tau        *,int                   > tauMap;
  std::     map<const pat ::Jet        *,int                   > jetMap;
  std::     map<const reco::Candidate  *,int                   > pcjMap;
  std::     map<const reco::Candidate  *,int                   > pfcMap;
//  std::     map<const reco::Track      *,int                   > tkmMap;
//  std::     map<const reco::Track      *,int                   > tkeMap;
//  std::     map<const reco::Track      *,int                   > tkpMap;
  std::     map<const reco::Track      *,int                   > ptjMap;
  std::     map<const reco::Track      *,int                   > tkvMap;
  std::     map<const reco::Track      *,int                   > trkMap;
  std::     map<const reco::Candidate  *,int                   > genMap;
  std::     map<const reco::GenJet     *,int                   > gpjMap;
  std::     set<const reco::Track      *                       > tkrSet;
  std::     set<const reco::Track      *                       > allPTk;
  std::multimap<const pat ::Muon       *,const reco::Track    *> mtkMap;
  std::multimap<const pat ::Muon       *,const reco::Candidate*> mrcMap;
  std::     map<const pat ::Electron   *,const reco::Track    *> etkMap;
  std::multimap<const reco::Candidate  *,const reco::Track    *> ptkMap;
  std::     map<const reco::Track      *,const pat ::Muon     *> tmuMap;
  std::     map<const reco::Track      *,const pat ::Electron *> telMap;
  std::     map<const reco::Track      *,const reco::Candidate*> tpfMap;
//  std::     map<const reco::Candidate  *,const reco::Vertex   *> mcvMap;
  std::multimap<const reco::Candidate  *,const pat ::Muon     *> cmuMap;
  std::     map<const reco::Candidate  *,const reco::Vertex   *> pcvMap;
  std::     map<const reco::Track      *,const reco::Vertex   *> tvpMap;
  std::multimap<const reco::Vertex     *,const reco::Track    *> vtkMap;
  std::  vector<const pat ::Muon       *                       > muonList;
  std::  vector<const reco::Candidate  *                       > muRCList;
  std::  vector<const reco::Candidate  *                       > candList;
  std::  vector<const reco::Track      *                       > trkpList;
//  std::  vector<PDEnumString::trackType                        > candType;
  std::     map<const BPHRecoCandidate *,int                   > rsvMap;
  std::  vector<const reco::Track      *                       > pTrack;
  std::  vector<const reco::Track      *                       > pTList;
  std::     map<const reco::Track      *,int                   > pTType;
  std::  vector<const reco::Track      *                       > mTrack;
  std::  vector<const reco::Track      *                       > mTList;
  std::     map<const reco::Track      *,int                   > mTType;
//  std::  vector<const reco::Track      *                       > pfcTrkP;
  std::  vector<const reco::Vertex     *                       > vtxList;
  std::     map<const reco::Vertex     *,int                   > vtxMap;
  std::     map<const reco::Vertex*,
                             std::vector<const reco::Track *> > vvtMap;

  std::set<int> linkedTrk;

  std::set<int> savedTrigObjType;
  std::set<int> savedTrk;
  std::set<int> savedPFC;
  std::set<int> savedVtx;
  std::set<int> savedJet;

  typedef
  std::multimap<const reco::Candidate  *,const reco::Track    *>
     ::const_iterator ptk_iter;

  void build( const edm::ParameterSet& ps );

  bool dumpNtuple;
  void fillHLTStatus   ();
  void fillHLTObjects  ();
  void fillBeamSpot    ();
  void fillMet         ();
  void fillOniaMuons   ();
  void fillMuons       ();
  void fillElectrons   ();
  void fillTaus        ();
  void fillJets        ();
  void fillVtxTrkMap   ();
  void fillPFCandidates();
  void fillTracks      ();
  void fillPVertices   ();
  void fillSVertices   ();
  void fillSVertReco   ();
  void fillPUWeight    ();
  void fillEvtInfo     ();
  void fillGenParticles();
  void fillGenJets     ();
  void fillUserInfo( int obj, int index,
                     int info , number value );
  void fillSVTInfoMap();
  void fillSVTIJetMap();

  void fillSaved();
  void cleanTracks();
  void cleanPFCand();

  const reco::Vertex* nearestPV( const reco::Track* );

  static void vertexDir( const reco::Vertex& pvtx,
                         const reco::Vertex& svtx,
                         const pat::CompositeCandidate& cand,
                         GlobalVector & dir,
                         Measurement1D& d2d,
                         Measurement1D& d3d );

  const reco::Vertex* vertexDir( const reco::Vertex& svtx,
                                 const pat::CompositeCandidate& cand,
                                 GlobalVector & dir,
                                 Measurement1D& d2d,
                                 Measurement1D& d3d );

  int addSecondaryVertex( const PDVertexWrapBase& vtx,
                          const reco::Vertex* pvt,
                          double vMass,
                          const std::map<const reco::         Track*,
                                               reco::TransientTrack*>& ttm,
                          const std::map<const reco::         Track*,
                                               float                >& tmm,
                          const  GlobalVector& dir,
                          const Measurement1D& d2d,
                          const Measurement1D& d3d,
                          const std::string& type, int jetId );
  int addSecondaryVertex( const PDVertexWrapBase& vtx,
                          int pvtIndex,
                          double vMass,
                          const std::map<const reco::         Track*,
                                               reco::TransientTrack*>& ttm,
                          const std::map<const reco::         Track*,
                                               float                >& tmm,
//                          const std::map<const reco::         Track*,
//                                         const reco::Candidate     *>& rcm,
//                          const std::map<const reco::         Track*,
//                                         const reco::Candidate     *>& rom,
                          const  GlobalVector& dir,
                          const Measurement1D& d2d,
                          const Measurement1D& d3d,
                          const std::string& type, int jetId );
  int addTrackIP(   int trkIndex, const reco::Track & trk, float mass,
                    int vtxIndex, const PDVertexWrapBase& vtx );
  int addTrackVtxP( int tipIndex, const reco::TransientTrack* trk,
                                  const PDVertexWrapBase& vtx );

  void refitPV( const reco::Vertex* pvt,
                const BPHRecoCandidate* cand, int svtIndex );
  void refitPV( const std::vector<reco::TransientTrack>& trTracks,
                const reco::BeamSpot* bs,
                const std::map<const reco::Track*,reco::TransientTrack*> ttm,
                const std::map<const reco::Track*,float                > tmm,
                int pvtIndex, int svtIndex, const std::string& type );

  void linkMTracks();
  void linkETracks();
  void linkPTracks();

  void fillPerigeePar();

  int nearestHLT( PDEnumString::trigObject type,
                  double pt, double eta, double phi );

  bool wildMatch( const std::string& name, std::string model );

  void setTrkType( const reco::Track* trk, int type,
                   std::vector<const reco::Track *   >& cTList,
                   std::map   <const reco::Track*,int>& cTType );
  int lastTrack( int i );

//  void fillTrigEvent();
//  void fillTrigObjects();
//  int fillTrigObject( int type,
//                      const pat::TriggerObject          & obj );
  int fillTrigObject( const pat::TriggerObjectRef         ref );
//  int fillTrigObject( const pat::TriggerObject          & obj );
  int fillTrigObject(       pat::TriggerObjectStandAlone  obj );

  template<class Container,class Token>
  bool fillPFList( Container& container, const Token& token,
                   PDEnumString::trackType type ) {
    token.get( currentEvBase, container );
    if ( !container.isValid() ) {
      switch ( type ) {
      case PDEnumString::pflink:
        std::cout << "invalid particle flow" << std::endl;
        break;
      case PDEnumString::packed:
        std::cout << "invalid packed candidates" << std::endl;
        break;
      case PDEnumString::gplink:
        std::cout << "invalid generic candidates" << std::endl;
        break;
      default:
        break;
      }
      return false;
    }
    int pos = 0;
    int iPF = candList.size();
    int lPF = iPF + container->size();
    candList.resize( lPF );
    while ( iPF < lPF ) candList[iPF++] = &container->at( pos++ );
/*
    candType.resize( lPF );
    while ( iPF < lPF ) {
      candList[iPF] = &container->at( pos++ );
      candType[iPF] = type;
      ++iPF;
    }
*/
    return true;
  }

};

#endif
