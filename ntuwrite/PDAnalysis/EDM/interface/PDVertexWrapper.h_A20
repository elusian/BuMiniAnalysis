#ifndef PDVertexWrapper_H
#define PDVertexWrapper_H
/** \class PDVertexWrapper
 *
 *  Description: 
 *
 *
 *  $Date: 2015-12-30 11:15:56 $
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


//---------------
// C++ Headers --
//---------------
#include <iostream>

//              ---------------------
//              -- Class Interface --
//              ---------------------

class PDVertexWrapBase {

 public:

  PDVertexWrapBase() {}
  virtual ~PDVertexWrapBase() {}

  /** typedefs
   */
  typedef reco::TemplatedSecondaryVertexTagInfo<reco::IPTagInfo<std::vector<edm::Ptr<reco::Candidate> >,reco::JetTagInfo>,reco::VertexCompositePtrCandidate> svtInfo;
  typedef reco::VertexCompositePtrCandidate sVtx;
  typedef std::vector< edm::Ptr<reco::Candidate> > svtTkList;

  struct VtxComponent {
    const reco::Candidate* cand;
    const reco::Track* track;
  };

  /** Operations
   */
  virtual bool isFake () const = 0;
  virtual bool isValid() const = 0;
  virtual const reco::Candidate::Point& position() const = 0;
  virtual double covariance( int i, int j ) const = 0;
  virtual const math::XYZTLorentzVector p4() const = 0;
  virtual unsigned int nTracks( float minWeight = 0.5 ) const = 0;
  virtual double chi2() const = 0;
  virtual double ndof() const = 0;
  virtual double normalizedChi2() const = 0;
  virtual PDVertexWrapBase* clone() const = 0;

  static void fillComponent( const svtInfo& secVtxTagInfo, int iVtx,
                             std::vector<VtxComponent>& compList ) {
    const svtTkList& tks = secVtxTagInfo.vertexTracks( iVtx );
    int itk = 0;
    int ntk = tks.size();
    compList.resize( ntk );
    for ( itk = 0; itk < ntk; ++itk ) {
      VtxComponent& vc = compList[itk];
      const reco::Candidate* pcp = tks[itk].get();
      const reco::Track* tkp = BPHTrackReference::getTrack( *pcp, "cbfp" );
      vc.cand = pcp;
      vc.track = tkp;
    }
    return;
  }
  static const reco::Vertex* vertexRef( const reco::Candidate* cand ) {
    const pat::PackedCandidate* pkc =
          dynamic_cast<const pat::PackedCandidate*>( cand );
    if ( pkc == 0 ) return 0;
    const reco::VertexRef vtr = pkc->vertexRef();
    try {
      return ( vtr.isNonnull() && vtr.isAvailable() ? vtr.get() : 0 );
    }
    catch ( edm::Exception e ) {
      return 0;
    }
  }

 private:

  PDVertexWrapBase           ( const PDVertexWrapBase& x );
  PDVertexWrapBase& operator=( const PDVertexWrapBase& x );

};


template<class T>
class PDVertexWrapper: public PDVertexWrapBase {

 public:

  /** Constructor
   */
  PDVertexWrapper( const T& v ): vtx( &v ) {}

  /** Destructor
   */
  virtual ~PDVertexWrapper() {}

  /** Operations
   */
  /// get vertex informations
  bool isFake () const { return isFake ( vtx ); }
  bool isValid() const { return isValid( vtx ); }
  const reco::Candidate::Point& position() const { return vtx->position(); }
  double covariance( int i, int j ) const { return covariance( vtx, i, j ); }
  const math::XYZTLorentzVector p4() const { return vtx->p4(); }
  unsigned int nTracks( float minWeight = 0.5 ) const { return nTracks( vtx, minWeight ); }
  double chi2() const { return chi2( vtx ); }
  double ndof() const { return ndof( vtx ); }
  double normalizedChi2() const { return normalizedChi2( vtx ); }
  PDVertexWrapBase* clone() const { return new PDVertexWrapper<T>( *vtx ); }

 private:

  // private copy and assigment constructors
  PDVertexWrapper           ( const PDVertexWrapper& x );
  PDVertexWrapper& operator=( const PDVertexWrapper& x );

  const T* vtx;

  bool isFake(  const reco::Vertex* v ) const { return v->isFake(); };
  bool isValid( const reco::Vertex* v ) const { return v->isValid(); };
  double covariance( const reco::Vertex* v, int i, int j ) const { return v->covariance( i, j ); }
  unsigned int nTracks( const reco::Vertex* v, float minWeight ) const { return vtx->nTracks( minWeight ); }
  double chi2( const reco::Vertex* v ) const { return v->chi2(); }
  double ndof( const reco::Vertex* v ) const { return v->ndof(); }
  double normalizedChi2( const reco::Vertex* v ) const { return v->normalizedChi2(); }

  bool isFake( const reco::VertexCompositePtrCandidate* v ) const { return false; };
  bool isValid( const reco::VertexCompositePtrCandidate* v ) const { return true; };
  double covariance( const reco::VertexCompositePtrCandidate* v, int i, int j ) const { return v->vertexCovariance( i, j ); }
  unsigned int nTracks( const reco::VertexCompositePtrCandidate* v, float minWeight ) const { return vtx->numberOfDaughters(); }
  double chi2( const reco::VertexCompositePtrCandidate* v ) const { return v->vertexChi2(); }
  double ndof( const reco::VertexCompositePtrCandidate* v ) const { return v->vertexNdof(); }
  double normalizedChi2( const reco::VertexCompositePtrCandidate* v ) const { return v->vertexNormalizedChi2(); }

};

class PDVertexContainer {
 public:
  PDVertexContainer(): vwb( 0 ) {} 
  PDVertexContainer           ( const PDVertexContainer& x ) { vwb = x.vwb->clone(); }
  PDVertexContainer& operator=( const PDVertexContainer& x ) { vwb = x.vwb->clone();
                                                               return *this; }
  ~PDVertexContainer() { delete vwb; }
  template<class T> void set( const T& v ) {
    delete vwb;
    vwb = new PDVertexWrapper<T>( v );
  }
  PDVertexWrapBase* get() const { return vwb; }
 private:
  PDVertexWrapBase* vwb;
};

#endif // PDVertexWrapper_H

