#ifndef NtuEDMAnalyzer_h
#define NtuEDMAnalyzer_h

#define UTIL_USE FULL

#include "NtuAnalysis/Common/interface/NtuEventHeader.h"
#include "NtuAnalysis/Common/interface/NtuAnalyzerUtil.h"
#include "NtuAnalysis/Write/interface/NtuEDConsumer.h"
#include "NtuTool/Common/interface/TreeWrapper.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"

#include <string>
#include <vector>
#include <map>
class TFile;

#define GET_PARAMETER( X, D ) getParameter( #X, X, D, true )
#define GET_UNTRACKED( X, D ) getParameter( #X, X, D, false )

class NtuEDMAnalyzer: public virtual NtuEventHeader,
                      public virtual NtuAnalyzerUtil,
                      public virtual TreeWrapper {

 protected:

  NtuEDMAnalyzer( const edm::ParameterSet& ps );
  virtual ~NtuEDMAnalyzer();

  template<class T>
  void getParameter( const std::string& key, T& val,
                     const std::string& def, bool tracked ) {
    setUserParameter( key, def );
    T tmp;
    getUserParameter( key, tmp );
    val = ( parameterSet->exists( key ) ?
            ( tracked ? 
              parameterSet->getParameter<T>( key ) :
              parameterSet->getUntrackedParameter<T>( key ) ) : tmp );
    setUserParameter( key, val );
    return;
  }

  virtual void beginJob();
  virtual void beginRun();
  virtual bool analyzeEDM( const edm::EventBase& ev,
                           int ientry, int event_file );
  virtual void read( const edm::EventBase& ev );
  virtual void endRun();
  virtual void endJob();

  template<class Consumer, class Obj>
  static void consume( NtuEDConsumer<Consumer>* c,
                       typename NtuEDToken<Obj>::type& token,
                       const std::string& label ) {
    if ( label == "" ) return;
    edm::InputTag tag( label );
    if ( c != 0 ) c->template consume<Obj>( token, tag );
    return;
  }
  template<class Consumer, class Obj>
  static void consume( NtuEDConsumer<Consumer>* c,
                       typename NtuEDToken<Obj>::type& token,
                       const edm::InputTag tag ) {
    if ( c != 0 ) c->template consume<Obj>( token, tag );
    return;
  }

  template<class T>
  class ObjectConsumer {
   public:
    ObjectConsumer( NtuEDMAnalyzer* a, NtuEDConsumer<T>* c ): 
     analyzer( a ),
     consumer( c ) {}
    template<class Obj>
    void consume( NtuEDToken<Obj>& nt,
                  const std::string& label ) {
      analyzer->template consume<T,Obj>( consumer, nt.token, label );
      return;
    }
    template<class Obj>
    void consume( NtuEDToken<Obj>& nt,
                  const edm::InputTag tag ) {
      analyzer->template consume<T,Obj>( consumer, nt.token, tag );
      return;
    }
   private:
    NtuEDMAnalyzer* analyzer;
    NtuEDConsumer<T>* consumer;
  };
  template<class T>
  ObjectConsumer<T>& getConsumer( NtuEDConsumer<T>* c ) {
    static ObjectConsumer<T>* ptr = new ObjectConsumer<T>( this, c );
    return *ptr;
  }

 private:

  // dummy copy constructor and assignment
  NtuEDMAnalyzer           ( const NtuEDMAnalyzer& c );
  NtuEDMAnalyzer& operator=( const NtuEDMAnalyzer& c );

};

#endif
