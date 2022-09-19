#ifndef NtuEventSelect_H
#define NtuEventSelect_H

#include <set>
#include <string>

class NtuEventSelect {

 public:

  NtuEventSelect();
  virtual ~NtuEventSelect();

  virtual bool read( const std::string& name );
  virtual bool find( unsigned int run, unsigned int event );

 private:

  // dummy copy constructor and assignment
  NtuEventSelect           ( const NtuEventSelect& );
  NtuEventSelect& operator=( const NtuEventSelect& );

  struct ev_id {
    unsigned int run;
    unsigned int event;
  };
  class ev_compare {
   public:
    bool operator()( const ev_id& l, const ev_id& r ) const {
      if ( l.run   < r.run   ) return true;
      if ( l.run   > r.run   ) return false;
      if ( l.event < r.event ) return true;
      return false;
    }
    static ev_compare instance;
  };
  std::set<ev_id, ev_compare> eventSet;

};


#endif

