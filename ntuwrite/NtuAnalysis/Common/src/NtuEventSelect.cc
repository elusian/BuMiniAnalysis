#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>

#include "NtuAnalysis/Common/interface/NtuEventSelect.h"

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

NtuEventSelect::ev_compare NtuEventSelect::ev_compare::instance;

NtuEventSelect::NtuEventSelect() {
}


NtuEventSelect::~NtuEventSelect() {
}


bool NtuEventSelect::read( const string& name ) {
  if ( name == "" ) return false;
  ifstream file( name.c_str() );
  stringstream sstr;
  char* line = new char[1000];
  vector<ev_id> evList;
  evList.reserve( 100000 );
  bool flag = false;
  while ( file.getline( line, 1000 ) ) {
    flag = true;
    char* lptr = line;
    while ( *lptr == ' ' ) ++lptr;
    if ( ( *lptr < '0' ) ||
         ( *lptr > '9' ) ) continue;
    sstr.clear();
    sstr.str( lptr );
    ev_id evId;
    sstr >> evId.run >> evId.event;
    evList.push_back( evId );
  }
  sort( evList.begin(), evList.end(), ev_compare::instance );
  eventSet.insert( evList.begin(), evList.end() );
  return flag;
}


bool NtuEventSelect::find( unsigned int run, unsigned int event ) {
  if ( eventSet.size() == 0 ) return true;
  ev_id evId;
  evId.run   = run;
  evId.event = event;
  if ( eventSet.find( evId ) != eventSet.end() ) return true;
  return false;
}



