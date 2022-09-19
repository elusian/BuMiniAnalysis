#ifndef TMPEDMToNtuple_h
#define TMPEDMToNtuple_h

#include "NtuAnalysis/Write/interface/NtuEDMAnalyzer.h"
#include "NtuAnalysis/Write/interface/NtuEDConsumer.h"
#include "TMPAnalysis/Ntu/bin/TMPAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <string>
#include <vector>
class TFile;

namespace pat {
  class Muon;
  class Jet;
}

class TMPEDMToNtuple: public TMPAnalyzer,
                      public NtuEDMAnalyzer {

 public:

  template<class T>
  TMPEDMToNtuple( const edm::ParameterSet& ps, NtuEDConsumer<T>* c  ):
   NtuEDMAnalyzer( ps ) {
    // parse ParameterSet
    build( ps );
    // interface to allow uniform access to data in different CMSSW versions
    ObjectConsumer<T>& oc = getConsumer( c );
    oc.template consume< std::vector<pat::Muon> >( gt_muons,
                                                 labelMuons );
    oc.template consume< std::vector<pat::Jet > >( gt_jets ,
                                                 labelJets  );
  }
  virtual ~TMPEDMToNtuple();

  virtual void beginJob();
  virtual void beginRun();
  virtual void read( const edm::EventBase& ev );
  virtual void endRun();
  virtual void endJob();

 private:

  template<class PATObject> class CompareByPt {
   public:
    bool operator()( const PATObject* l, const PATObject* r ) {
      return l->pt() > r->pt();
    }
  };

  std::string labelMuons;
  std::string labelJets;

  edm::Handle< std::vector<pat::Muon> >    muons;
  edm::Handle< std::vector<pat::Jet > >    jets;

  // interfaces to allow uniform access to data in different CMSSW versions
  NtuEDToken < std::vector<pat::Muon> > gt_muons;
  NtuEDToken < std::vector<pat::Jet > > gt_jets;

  void build( const edm::ParameterSet& ps );

  void fillMuons();
  void fillJets();

};

#endif
