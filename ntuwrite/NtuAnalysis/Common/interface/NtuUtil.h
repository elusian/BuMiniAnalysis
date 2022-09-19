#ifndef NtuUtil_H
#define NtuUtil_H

#include "map"
#include "vector"
#include "math.h"

template<class T>
class NtuUtil {

 public:

  NtuUtil();
  virtual ~NtuUtil();

  static T energy( T px, T py, T pz, T mm );
  static T invMass( T px , T  py , T  pz,  T ee );
  static T invMass( T px1, T  py1, T  pz1, T ee1,
                    T px2, T  py2, T  pz2, T ee2 );
  static T invMass( T px1, T  py1, T  pz1, T ee1,
                    T px2, T  py2, T  pz2, T ee2,
                    T px3, T  py3, T  pz3, T ee3 );
  static T invMass( T px1, T  py1, T  pz1, T ee1,
                    T px2, T  py2, T  pz2, T ee2,
                    T px3, T  py3, T  pz3, T ee3,
                    T px4, T  py4, T  pz4, T ee4 );
  static T invMSph( T pt , T eta , T phi,  T ee );
  static T invMSph( T pt1, T eta1, T phi1, T ee1,
                    T pt2, T eta2, T phi2, T ee2 );
  static T invMSph( T pt1, T eta1, T phi1, T ee1,
                    T pt2, T eta2, T phi2, T ee2,
                    T pt3, T eta3, T phi3, T ee3 );
  static T invMSph( T pt1, T eta1, T phi1, T ee1,
                    T pt2, T eta2, T phi2, T ee2,
                    T pt3, T eta3, T phi3, T ee3,
                    T pt4, T eta4, T phi4, T ee4 );
  static T invMass_M( T px1, T  py1, T  pz1, T mm1,
                      T px2, T  py2, T  pz2, T mm2 );
  static T invMass_M( T px1, T  py1, T  pz1, T mm1,
                      T px2, T  py2, T  pz2, T mm2,
                      T px3, T  py3, T  pz3, T mm3 );
  static T invMass_M( T px1, T  py1, T  pz1, T mm1,
                      T px2, T  py2, T  pz2, T mm2,
                      T px3, T  py3, T  pz3, T mm3,
                      T px4, T  py4, T  pz4, T mm4 );
  static T invMSph_M( T pt1, T eta1, T phi1, T mm1,
                      T pt2, T eta2, T phi2, T mm2 );
  static T invMSph_M( T pt1, T eta1, T phi1, T mm1,
                      T pt2, T eta2, T phi2, T mm2,
                      T pt3, T eta3, T phi3, T mm3 );
  static T invMSph_M( T pt1, T eta1, T phi1, T mm1,
                      T pt2, T eta2, T phi2, T mm2,
                      T pt3, T eta3, T phi3, T mm3,
                      T pt4, T eta4, T phi4, T mm4 );
  static T dSqua( T a1, T b1,
                  T a2, T b2 );
  static T dSqua( T a1, T b1, T c1,
                  T a2, T b2, T c2 );
  static T delta( T a1, T b1,
                  T a2, T b2 );
  static T delta( T a1, T b1, T c1,
                  T a2, T b2, T c2 );
  static T dRSqua( T eta1, T phi1,
                   T eta2, T phi2 );
  static T deltaR( T eta1, T phi1,
                   T eta2, T phi2 );

  static T modSqua( T x, T y, T z );
  static T modSqua( T t, T e );
  static T modCart( T x, T y, T z )
    { return sqrt( modSqua( x, y, z ) ); }
  static T modSphe( T t, T e )
    { return sqrt( modSqua( t, e ) ); }

  static T angleCart( T x1, T y1, T z1,
                      T x2, T y2, T z2 );
  static T angleSphe( T t1, T e1, T p1,
                      T t2, T e2, T p2 );
  static T angle( T x1, T y1, T z1,
                  T x2, T y2, T z2 ) 
    { return angleCart( x1, y1, z1, x2, y2, z2 ); }
  static T angle( T e1, T p1,
                  T e2, T p2 )
    { return angleSphe( 1.0, e1, p1, 1.0, e2, p2 ); }

  static T ptCart( T px, T py, T pz,
                   T qx, T qy, T qz );
  static T ptSphe( T pt, T pe, T pp,
                         T qe, T qp );
  static T ptRel ( T px, T py, T pz,
                   T qx, T qy, T qz )
    { return ptCart( px, py, pz, qx, qy, qz ); }
  static T ptRel ( T pt, T pe, T pp,
                         T qe, T qp )
    { return ptSphe( pt, pe, pp, qe, qp ); }
  static T plCart( T px, T py, T pz,
                   T qx, T qy, T qz );
  static T plSphe( T pt, T pe, T pp,
                         T qe, T qp );
  static T plRel ( T px, T py, T pz,
                   T qx, T qy, T qz )
    { return plCart( px, py, pz, qx, qy, qz ); }
  static T plRel ( T pt, T pe, T pp,
                         T qe, T qp )
    { return plSphe( pt, pe, pp, qe, qp ); }

  static T transf( T  px, T  py, T  pz, T  ep,
                   T  qx, T  qy, T  qz, T  eq,
                   T& tx, T& ty, T& tz, T& et );

  static T ptCMin( T px, T py, T pz,
                   T qx, T qy, T qz );
  static T ptSMin( T pt, T pe, T pp,
                   T qt, T qe, T qp );
  static T ptCMax( T px, T py, T pz,
                   T qx, T qy, T qz );
  static T ptSMax( T pt, T pe, T pp,
                   T qt, T qe, T qp );

  static void convCartSphe( T  x, T  y, T  z,
                            T& t, T& e, T& p );
  static void convSpheCart( T  t, T  e, T  p,
                            T& x, T& y, T& z );
  static void convCartSphe( const std::vector<T>* x,
                            const std::vector<T>* y,
                            const std::vector<T>* z,
                                  std::vector<T>* t,
                                  std::vector<T>* e,
                                  std::vector<T>* p, int i = -1 );
  static void convSpheCart( const std::vector<T>* t,
                            const std::vector<T>* e,
                            const std::vector<T>* p,
                                  std::vector<T>* x,
                                  std::vector<T>* y,
                                  std::vector<T>* z, int i = -1 );

  struct ObjectAssociation {
    int iObj;
    float dist;
  };
  class ObjectSelection {
   public:
    virtual ~ObjectSelection() {}
    virtual bool operator()( int obj ) { return true; }
  };
//  static ObjectSelection dummySel;
  class ObjectDistance {
   public:
    virtual ~ObjectDistance() {}
    virtual T operator()( T lObjPt , T lObjEta, T lObjPhi,
                          T lObjPx , T lObjPy , T lObjPz ,
                          T lObjE  , int    lObjCha,
                          T rObjPt , T rObjEta, T rObjPhi,
                          T rObjPx , T rObjPy , T rObjPz ,
                          T rObjE  , int    rObjCha )
                          = 0;
    virtual T dMax() { return 1.0e+37; }
  };
  static void associateObjects( const std::vector<T>* lObjPt,
                                const std::vector<T>* lObjEta,
                                const std::vector<T>* lObjPhi,
                                const std::vector<T>* lObjPx,
                                const std::vector<T>* lObjPy,
                                const std::vector<T>* lObjPz,
                                const std::vector<T>* lObjE,
                                const std::vector<int   >* lObjCharge,
                                const std::vector<T>* rObjPt,
                                const std::vector<T>* rObjEta,
                                const std::vector<T>* rObjPhi,
                                const std::vector<T>* rObjPx,
                                const std::vector<T>* rObjPy,
                                const std::vector<T>* rObjPz,
                                const std::vector<T>* rObjE,
                                const std::vector<int   >* rObjCharge,
                                ObjectDistance* objDist,
                                ObjectSelection* lObjSel,
                                ObjectSelection* rObjSel,
                                bool singleAssociation,
                                std::vector<ObjectAssociation>* assoc );

 private:

  // dummy copy constructor and assignment
  NtuUtil           ( const NtuUtil& );
  NtuUtil& operator=( const NtuUtil& );

};


template<class Key,class Val>
class ConstMap {
 public:
  ConstMap<Key,Val>() { cmp = 0; }
  ConstMap<Key,Val>( const Val& v ) { cmp = 0; val = v; }
  ConstMap<Key,Val>( const std::map<Key,Val>& map ) { cmp = &map; }
  ConstMap<Key,Val>( const std::map<Key,Val>& map, const Val& v ) {
    cmp = &map; val = v; }
  const Val& operator[]( Key key ) {
    if ( cmp == 0 ) return val;
    typename std::map<Key,Val>::const_iterator iter = cmp->find( key );
    typename std::map<Key,Val>::const_iterator iend = cmp->end();
    if ( iter != iend ) return iter->second;
    return val;
  }
  ConstMap& operator=( ConstMap& map ) {
    cmp = map.cmp;
    return *this;
  }
  ConstMap& operator=( const std::map<Key,Val>& map ) {
    cmp = &map;
    return *this;
  }

 private:

  const std::map<Key,Val>* cmp;
  Val val;

};

#include "NtuAnalysis/Common/interface/NtuUtil.hpp"

#endif

