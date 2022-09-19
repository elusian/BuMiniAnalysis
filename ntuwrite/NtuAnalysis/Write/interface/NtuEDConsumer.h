#ifndef NtuEDConsumer_H
#define NtuEDConsumer_H
/** \class NtuEDConsumer
 *
 *  Description: 
 *    common interface to forward calls to EDConsumerBase
 *
 *  $Date: 2016-04-15 17:47:56 $
 *  $Revision: 1.1 $
 *  \author Paolo Ronchese INFN Padova
 *
 */

//----------------------
// Base Class Headers --
//----------------------


//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"

//---------------
// C++ Headers --
//---------------
#include <string>

//              ---------------------
//              -- Class Interface --
//              ---------------------

template<class Obj>
class NtuEDToken {
 public:
  typedef typename edm::EDGetTokenT<Obj> type;
  static bool get( const edm::EventBase* ev,
                   typename NtuEDToken<Obj>::type& t,
                   edm::Handle<Obj>& obj ) {
    const edm::Event* ep = dynamic_cast<const edm::Event*>( ev );
    if ( ep == 0 ) return false;
    return ep->getByToken( t, obj );
  }
  bool get( const edm::EventBase* ev,
            edm::Handle<Obj>& obj ) const {
    const edm::Event* ep = dynamic_cast<const edm::Event*>( ev );
    if ( ep == 0 ) return false;
    return ep->getByToken( token, obj );
  }
  type token;
};


template<class T>
class NtuEDConsumer: public virtual T {

 public:

  /** Constructor
   */
  NtuEDConsumer() {}

  /** Destructor
   */
  virtual ~NtuEDConsumer() {}

  /** Operations
   */
  /// forward call to EDConsumerBase
  template<class Obj>
  void consume( typename NtuEDToken<Obj>::type& token,
                const edm::InputTag& tag ) {
    token = this->template consumes<Obj>( tag );
    return;
  }

 private:

  // private copy and assigment constructors
  NtuEDConsumer           ( const NtuEDConsumer& x );
  NtuEDConsumer& operator=( const NtuEDConsumer& x );

};

//#include "NtuAnalysis/Write/interface/NtuEDConsumer.hpp"

#endif // NtuEDConsumer_H

