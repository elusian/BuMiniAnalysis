#include <iostream>
#include <sstream>
#include <string>
#include <set>
#include <math.h>

#include "TDirectory.h"
#include "TBranch.h"
#include "TTree.h"
#include "TCanvas.h"
//#include "Math/LorentzVector.h"

//NtuUtil::ObjectSelection NtuUtil::dummySel;

template<class T>
NtuUtil<T>::NtuUtil() {
}


template<class T>
NtuUtil<T>::~NtuUtil() {
}


template<class T>
T NtuUtil<T>::energy( T px, T py, T pz, T mm ) {
  return sqrt( ( px * px ) + ( py * py ) + ( pz * pz ) + ( mm * mm ) );
}


template<class T>
T NtuUtil<T>::invMass( T px, T py, T pz, T ee ) {
  T mq = ( px * px ) + ( py * py ) + ( pz * pz ) - ( ee * ee );
  if ( mq > 0.0 ) return 0.0;
  return sqrt( -mq );
}


template<class T>
T NtuUtil<T>::invMass( T px1, T py1, T pz1, T ee1,
                       T px2, T py2, T pz2, T ee2 ) {
  T pxt = px1 + px2;
  T pyt = py1 + py2;
  T pzt = pz1 + pz2;
  T eet = ee1 + ee2;
  return sqrt( ( eet * eet ) -
             ( ( pxt * pxt ) + ( pyt * pyt ) + ( pzt * pzt ) ) );
}


template<class T>
T NtuUtil<T>::invMass( T px1, T py1, T pz1, T ee1,
                       T px2, T py2, T pz2, T ee2,
                       T px3, T py3, T pz3, T ee3 ) {
  T pxt = px1 + px2 + px3;
  T pyt = py1 + py2 + py3;
  T pzt = pz1 + pz2 + pz3;
  T eet = ee1 + ee2 + ee3;
  return sqrt( ( eet * eet ) -
             ( ( pxt * pxt ) + ( pyt * pyt ) + ( pzt * pzt ) ) );
}


template<class T>
T NtuUtil<T>::invMass( T px1, T py1,T pz1, T ee1,
                       T px2, T py2, T pz2, T ee2,
                       T px3, T py3, T pz3, T ee3,
                       T px4, T py4, T pz4, T ee4 ) {
  T pxt = px1 + px2 + px3 + px4;
  T pyt = py1 + py2 + py3 + py4;
  T pzt = pz1 + pz2 + pz3 + pz4;
  T eet = ee1 + ee2 + ee3 + ee4;
  return sqrt( ( eet * eet ) -
             ( ( pxt * pxt ) + ( pyt * pyt ) + ( pzt * pzt ) ) );
}


template<class T>
T NtuUtil<T>::invMSph( T  pt, T eta, T phi, T  ee ) {
  T px;
  T py;
  T pz;
  convSpheCart( pt, eta, phi, px, py, pz );
  T mq = ( px * px ) + ( py * py ) + ( pz * pz ) - ( ee * ee );
  if ( mq > 0.0 ) return 0.0;
  return sqrt( -mq );
}


template<class T>
T NtuUtil<T>::invMSph( T  pt1, T eta1, T phi1, T  ee1,
                       T  pt2, T eta2, T phi2, T  ee2 ) {
  T px1;
  T py1;
  T pz1;
  T px2;
  T py2;
  T pz2;
  convSpheCart( pt1, eta1, phi1, px1, py1, pz1 );
  convSpheCart( pt2, eta2, phi2, px2, py2, pz2 );
  return invMass( px1, py1, pz1, ee1,
                  px2, py2, pz2, ee2 );
}


template<class T>
T NtuUtil<T>::invMSph( T  pt1, T eta1, T phi1, T  ee1,
                       T  pt2, T eta2, T phi2, T  ee2,
                       T  pt3, T eta3, T phi3, T  ee3 ) {
  T px1;
  T py1;
  T pz1;
  T px2;
  T py2;
  T pz2;
  T px3;
  T py3;
  T pz3;
  convSpheCart( pt1, eta1, phi1, px1, py1, pz1 );
  convSpheCart( pt2, eta2, phi2, px2, py2, pz2 );
  convSpheCart( pt3, eta3, phi3, px3, py3, pz3 );
  return invMass( px1, py1, pz1, ee1,
                  px2, py2, pz2, ee2,
                  px3, py3, pz3, ee3 );
}


template<class T>
T NtuUtil<T>::invMSph( T  pt1, T eta1, T phi1, T  ee1,
                       T  pt2, T eta2, T phi2, T  ee2,
                       T  pt3, T eta3, T phi3, T  ee3,
                       T  pt4, T eta4, T phi4, T  ee4 ) {
  T px1;
  T py1;
  T pz1;
  T px2;
  T py2;
  T pz2;
  T px3;
  T py3;
  T pz3;
  T px4;
  T py4;
  T pz4;
  convSpheCart( pt1, eta1, phi1, px1, py1, pz1 );
  convSpheCart( pt2, eta2, phi2, px2, py2, pz2 );
  convSpheCart( pt3, eta3, phi3, px3, py3, pz3 );
  convSpheCart( pt4, eta4, phi4, px4, py4, pz4 );
  return invMass( px1, py1, pz1, ee1,
                  px2, py2, pz2, ee2,
                  px3, py3, pz3, ee3,
                  px4, py4, pz4, ee4 );
}


template<class T>
T NtuUtil<T>::invMass_M( T px1, T py1, T pz1, T mm1,
                         T px2, T py2, T pz2, T mm2 ) {
  T pxt = px1 + px2;
  T pyt = py1 + py2;
  T pzt = pz1 + pz2;
  T eet = energy( px1, py1, pz1, mm1 ) +
          energy( px2, py2, pz2, mm2 );
  return sqrt( ( eet * eet ) -
             ( ( pxt * pxt ) + ( pyt * pyt ) + ( pzt * pzt ) ) );
}


template<class T>
T NtuUtil<T>::invMass_M( T px1, T py1, T pz1, T mm1,
                         T px2, T py2, T pz2, T mm2,
                         T px3, T py3, T pz3, T mm3 ) {
  T pxt = px1 + px2 + px3;
  T pyt = py1 + py2 + py3;
  T pzt = pz1 + pz2 + pz3;
  T eet = energy( px1, py1, pz1, mm1 ) +
          energy( px2, py2, pz2, mm2 ) +
          energy( px3, py3, pz3, mm3 );
  return sqrt( ( eet * eet ) -
             ( ( pxt * pxt ) + ( pyt * pyt ) + ( pzt * pzt ) ) );
}


template<class T>
T NtuUtil<T>::invMass_M( T px1, T py1, T pz1, T mm1,
                         T px2, T py2, T pz2, T mm2,
                         T px3, T py3, T pz3, T mm3,
                         T px4, T py4, T pz4, T mm4 ) {
  T pxt = px1 + px2 + px3 + px4;
  T pyt = py1 + py2 + py3 + py4;
  T pzt = pz1 + pz2 + pz3 + pz4;
  T eet = energy( px1, py1, pz1, mm1 ) +
          energy( px2, py2, pz2, mm2 ) +
          energy( px3, py3, pz3, mm3 ) +
          energy( px4, py4, pz4, mm4 );
  return sqrt( ( eet * eet ) -
             ( ( pxt * pxt ) + ( pyt * pyt ) + ( pzt * pzt ) ) );
}


template<class T>
T NtuUtil<T>::invMSph_M( T  pt1, T eta1, T phi1, T  mm1,
                         T  pt2, T eta2, T phi2, T  mm2 ) {
  T px1;
  T py1;
  T pz1;
  T px2;
  T py2;
  T pz2;
  convSpheCart( pt1, eta1, phi1, px1, py1, pz1 );
  convSpheCart( pt2, eta2, phi2, px2, py2, pz2 );
  return invMass_M( px1, py1, pz1, mm1,
                    px2, py2, pz2, mm2 );
}


template<class T>
T NtuUtil<T>::invMSph_M( T  pt1, T eta1, T phi1, T  mm1,
                         T  pt2, T eta2, T phi2, T  mm2,
                         T  pt3, T eta3, T phi3, T  mm3 ) {
  T px1;
  T py1;
  T pz1;
  T px2;
  T py2;
  T pz2;
  T px3;
  T py3;
  T pz3;
  convSpheCart( pt1, eta1, phi1, px1, py1, pz1 );
  convSpheCart( pt2, eta2, phi2, px2, py2, pz2 );
  convSpheCart( pt3, eta3, phi3, px3, py3, pz3 );
  return invMass_M( px1, py1, pz1, mm1,
                    px2, py2, pz2, mm2,
                    px3, py3, pz3, mm3 );
}


template<class T>
T NtuUtil<T>::invMSph_M( T  pt1, T eta1, T phi1, T  mm1,
                         T  pt2, T eta2, T phi2, T  mm2,
                         T  pt3, T eta3, T phi3, T  mm3,
                         T  pt4, T eta4, T phi4, T  mm4 ) {
  T px1;
  T py1;
  T pz1;
  T px2;
  T py2;
  T pz2;
  T px3;
  T py3;
  T pz3;
  T px4;
  T py4;
  T pz4;
  convSpheCart( pt1, eta1, phi1, px1, py1, pz1 );
  convSpheCart( pt2, eta2, phi2, px2, py2, pz2 );
  convSpheCart( pt3, eta3, phi3, px3, py3, pz3 );
  convSpheCart( pt4, eta4, phi4, px4, py4, pz4 );
  return invMass_M( px1, py1, pz1, mm1,
                    px2, py2, pz2, mm2,
                    px3, py3, pz3, mm3,
                    px4, py4, pz4, mm4 );
}


template<class T>
T NtuUtil<T>::dSqua( T a1, T b1, T a2, T b2 ) {
  T aD = a1 - a2;
  T bD = b1 - b2;
  return ( ( aD * aD ) + ( bD * bD ) );
}


template<class T>
T NtuUtil<T>::dSqua( T a1, T b1, T c1,
                     T a2, T b2, T c2 ) {
  T aD = a1 - a2;
  T bD = b1 - b2;
  T cD = c1 - c2;
  return ( ( aD * aD ) + ( bD * bD ) + ( cD * cD ) );
}


template<class T>
T NtuUtil<T>::delta( T a1, T b1,
                     T a2, T b2 ) {
  T aD = a1 - a2;
  T bD = b1 - b2;
  return sqrt( ( aD * aD ) + ( bD * bD ) );
}


template<class T>
T NtuUtil<T>::delta( T a1, T b1, T c1,
                     T a2, T b2, T c2 ) {
  T aD = a1 - a2;
  T bD = b1 - b2;
  T cD = c1 - c2;
  return sqrt( ( aD * aD ) + ( bD * bD ) + ( cD * cD ) );
}

template<class T>
T NtuUtil<T>::dRSqua( T eta1, T phi1,
                      T eta2, T phi2 ) {
  T etaD = eta1 - eta2;
  T phiD = phi1 - phi2;
  phiD = fabs( M_PI - fabs( M_PI - fabs( phiD ) ) );
  return ( etaD * etaD ) + ( phiD * phiD );
}

template<class T>
T NtuUtil<T>::deltaR( T eta1, T phi1,
                      T eta2, T phi2 ) {
  T etaD = eta1 - eta2;
  T phiD = phi1 - phi2;
  phiD = fabs( M_PI - fabs( M_PI - fabs( phiD ) ) );
  return sqrt( ( etaD * etaD ) + ( phiD * phiD ) );
}

template<class T>
T NtuUtil<T>::modSqua( T x, T y, T z ) {
  return ( x * x ) + ( y * y ) + ( z * z );
}

template<class T>
T NtuUtil<T>::modSqua( T t, T e ) {
  T z = t / tan( 2 * atan( exp( -e ) ) );
  return ( t * t ) + ( z * z );
}

template<class T>
T NtuUtil<T>::angleCart( T x1, T y1, T z1,
                         T x2, T y2, T z2 ) {
  T s = ( x1 * x2 ) + ( y1 * y2 ) + ( z1 * z2 );
  T p = modCart( x1, y1, z1 ) * modCart( x2, y2, z2 );
  if ( fabs( p ) <= fabs( s ) ) return ( s * p > 0 ? 0 : M_PI );
  return acos( s / p );
}

template<class T>
T NtuUtil<T>::angleSphe( T t1, T e1, T p1,
                         T t2, T e2, T p2 ) {
  T x1;
  T y1;
  T z1;
  T x2;
  T y2;
  T z2;
  convSpheCart( t1, e1, p1, x1, y1, z1 );
  convSpheCart( t2, e2, p2, x2, y2, z2 );
  return angleCart( x1, y1, z1, x2, y2, z2 );
}

template<class T>
T NtuUtil<T>::ptCart( T px, T py, T pz,
                      T qx, T qy, T qz ) {
  T s = ( px * qx ) + ( py * qy ) + ( pz * qz );
  T p = modSqua( px, py, pz );
  T q = modSqua( qx, qy, qz );
  T z = ( s * s ) / ( p * q );
  T f = ( z < 1.0 ? sqrt( 1.0 - z ) : 0.0 );
  return sqrt( p ) * f;
}

template<class T>
T NtuUtil<T>::ptSphe( T pt, T pe, T pp,
                      T qe, T qp ) {
  T px;
  T py;
  T pz;
  T qx;
  T qy;
  T qz;
  convSpheCart(  pt, pe, pp, px, py, pz );
  convSpheCart( 1.0, qe, qp, qx, qy, qz );
  return ptCart( px, py, pz, qx, qy, qz );
}

template<class T>
T NtuUtil<T>::plCart( T px, T py, T pz,
                      T qx, T qy, T qz ) {
  return ( ( px * qx ) + ( py * qy ) + ( pz * qz ) ) / modCart( qx, qy, qz );
}

template<class T>
T NtuUtil<T>::plSphe( T pt, T pe, T pp,
                      T qe, T qp ) {
  T px;
  T py;
  T pz;
  T qx;
  T qy;
  T qz;
  convSpheCart(  pt, pe, pp, px, py, pz );
  convSpheCart( 1.0, qe, qp, qx, qy, qz );
  return plCart( px, py, pz, qx, qy, qz );
}

template<class T>
T NtuUtil<T>::transf( T  px, T  py, T  pz, T  ep,
                      T  qx, T  qy, T  qz, T  eq,
                      T& tx, T& ty, T& tz, T& et ) {
  T s = ( px * qx ) + ( py * qy ) + ( pz * qz );
  T p = modSqua( px, py, pz );
  T q = modSqua( qx, qy, qz );
  T z = ( s * s ) / ( p * q );
  T f = ( z < 1.0 ? sqrt( 1.0 - z ) : 0.0 );
  p = sqrt( p );
  q = sqrt( q );
  T b = q / eq;
  T g = 1.0 / sqrt( 1.0 - ( b * b ) );
  T pt = p * f;
  T pl = s / q;
  T tl = g * ( pl - ( b * ep ) );
  T at = atan( pt / tl );
  et = g * ( ep - ( b * pl ) );
  if ( at < 0.0 ) at += M_PI;
  tx = px + ( ( tl - pl ) * qx / q );
  ty = py + ( ( tl - pl ) * qy / q );
  tz = pz + ( ( tl - pl ) * qz / q );
  return at;
}

template<class T>
T NtuUtil<T>::ptCMin( T px, T py, T pz,
                      T qx, T qy, T qz ) {
  T s = ( px * qx ) + ( py * qy ) + ( pz * qz );
  T p = modSqua( px, py, pz );
  T q = modSqua( qx, qy, qz );
  T t = p * q;
  T f = ( fabs ( t ) > fabs( s ) ? sqrt( 1.0 - ( s * s / t ) ) : 1.0 );
  return ( p < q ? sqrt( p ) * f : sqrt( q ) * f );
}

template<class T>
T NtuUtil<T>::ptSMin( T pt, T pe, T pp,
                      T qt, T qe, T qp ) {
  T px;
  T py;
  T pz;
  T qx;
  T qy;
  T qz;
  convSpheCart(  pt, pe, pp, px, py, pz );
  convSpheCart(  qt, qe, qp, qx, qy, qz );
  return ptCMin( px, py, pz, qx, qy, qz );
}

template<class T>
T NtuUtil<T>::ptCMax( T px, T py, T pz,
                      T qx, T qy, T qz ) {
  T s = ( px * qx ) + ( py * qy ) + ( pz * qz );
  T p = modSqua( px, py, pz );
  T q = modSqua( qx, qy, qz );
  T t = p * q;
  T f = ( fabs ( t ) > fabs( s ) ? sqrt( 1.0 - ( s * s / t ) ) : 1.0 );
  return ( p > q ? sqrt( p ) * f : sqrt( q ) * f );
}

template<class T>
T NtuUtil<T>::ptSMax( T pt, T pe, T pp,
                      T qt, T qe, T qp ) {
  T px;
  T py;
  T pz;
  T qx;
  T qy;
  T qz;
  convSpheCart(  pt, pe, pp, px, py, pz );
  convSpheCart(  qt, qe, qp, qx, qy, qz );
  return ptCMax( px, py, pz, qx, qy, qz );
}

template<class T>
void NtuUtil<T>::convCartSphe( T  x, T  y, T  z,
                               T& t, T& e, T& p ) {
  t = sqrt( ( x * x ) + ( y * y ) );
  T
  r = sqrt( ( x * x ) + ( y * y ) + ( z * z ) );
  e =  log( ( r + z ) / ( r - z ) ) / 2.0;
  p = atan2( y, x );
  return;
}


template<class T>
void NtuUtil<T>::convSpheCart( T  t, T  e, T  p,
                               T& x, T& y, T& z ) {
  x = t * cos( p );
  y = t * sin( p );
  z = t / tan( 2 * atan( exp( -e ) ) );
  return;
}


template<class T>
void NtuUtil<T>::convCartSphe( const std::vector<T>* x,
                               const std::vector<T>* y,
                               const std::vector<T>* z,
                                     std::vector<T>* t,
                                     std::vector<T>* e,
                                     std::vector<T>* p,
                               int i ) {
  if ( i < 0 ) {
    int j;
    int n = x->size();
    t->resize( n );
    e->resize( n );
    p->resize( n );
    for ( j = 0; j < n; ++j )
    convCartSphe( x->at( j ), y->at( j ), z->at( j ),
                  t->at( j ), e->at( j ), p->at( j ) );
  }
  else {
    convCartSphe( x->at( i ), y->at( i ), z->at( i ),
                  t->at( i ), e->at( i ), p->at( i ) );
  }
  return;
}


template<class T>
void NtuUtil<T>::convSpheCart( const std::vector<T>* t,
                               const std::vector<T>* e,
                               const std::vector<T>* p,
                                     std::vector<T>* x,
                                     std::vector<T>* y,
                                     std::vector<T>* z,
                               int i ) {
  if ( i < 0 ) {
    int j;
    int n = t->size();
    x->resize( n );
    y->resize( n );
    z->resize( n );
    for ( j = 0; j < n; ++j )
    convSpheCart( t->at( j ), e->at( j ), p->at( j ),
                  x->at( j ), y->at( j ), z->at( j ) );
  }
  else {
    convSpheCart( t->at( i ), e->at( i ), p->at( i ),
                  x->at( i ), y->at( i ), z->at( i ) );
  }
  return;
}


template<class T>
void NtuUtil<T>::associateObjects( const std::vector<T>* lObjPt,
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
                                   std::vector<ObjectAssociation>* assoc ) {
  int nl = lObjPt->size();
  int nr = rObjPt->size();
  int il;
  int ir;
  T lPt;
  T lEta;
  T lPhi;
  T lPx;
  T lPy;
  T lPz;
  T lE;
  int    lCharge;
  T rPt;
  T rEta;
  T rPhi;
  T rPx;
  T rPy;
  T rPz;
  T rE;
  int    rCharge;
  T dMin;
  T dCur;
  int iAss;
  ObjectAssociation noAss;
  noAss.iObj = -1;
  noAss.dist = objDist->dMax();
  assoc->clear();
  assoc->resize( nl, noAss );
  std::set<int> used;
  for ( il = 0; il < nl; ++il ) {
    if ( ( lObjSel != 0 ) && !(*lObjSel)( il ) ) continue;
    iAss = -1;
    dMin = objDist->dMax();
    lPt     = ( lObjPt    != 0 ? lObjPt    ->at( il ) : 0.0 );
    lEta    = ( lObjEta   != 0 ? lObjEta   ->at( il ) : 0.0 );
    lPhi    = ( lObjPhi   != 0 ? lObjPhi   ->at( il ) : 0.0 );
    lPx     = ( lObjPx    != 0 ? lObjPx    ->at( il ) : 0.0 );
    lPy     = ( lObjPy    != 0 ? lObjPy    ->at( il ) : 0.0 );
    lPz     = ( lObjPz    != 0 ? lObjPz    ->at( il ) : 0.0 );
    lE      = ( lObjE     != 0 ? lObjE     ->at( il ) : 0.0 );
    lCharge = ( lObjCharge!= 0 ? lObjCharge->at( il ) : 0   );
    for ( ir = 0; ir < nr; ++ir ) {
      if ( ( rObjSel != 0 ) && !(*rObjSel)( ir ) ) continue;
      rPt     = ( rObjPt    != 0 ? rObjPt    ->at( ir ) : 0.0 );
      rEta    = ( rObjEta   != 0 ? rObjEta   ->at( ir ) : 0.0 );
      rPhi    = ( rObjPhi   != 0 ? rObjPhi   ->at( ir ) : 0.0 );
      rPx     = ( rObjPx    != 0 ? rObjPx    ->at( ir ) : 0.0 );
      rPy     = ( rObjPy    != 0 ? rObjPy    ->at( ir ) : 0.0 );
      rPz     = ( rObjPz    != 0 ? rObjPz    ->at( ir ) : 0.0 );
      rE      = ( rObjE     != 0 ? rObjE     ->at( ir ) : 0.0 );
      rCharge = ( rObjCharge!= 0 ? rObjCharge->at( ir ) : 0   );
      if ( singleAssociation && ( used.find( ir ) != used.end() ) ) continue;
      dCur = (*objDist)(    lPt    ,    lEta       ,    lPhi    ,
                            lPx    ,    lPy        ,    lPz     ,
                            lE     ,    lCharge    ,
                            rPt    ,    rEta       ,    rPhi    ,
                            rPx    ,    rPy        ,    rPz     ,
                            rE     ,    rCharge );
      if ( dCur < dMin ) {
        iAss = ir;
        dMin = dCur;
      }
    }
    if ( iAss < 0 ) continue;
    used.insert( iAss );
    assoc->at( il ).dist = dMin;
    assoc->at( il ).iObj = iAss;
  }
  return;
}

