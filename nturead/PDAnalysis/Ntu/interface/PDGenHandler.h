#ifndef PDGenHandler_h
#define PDGenHandler_h

#include <vector>
#include <set>
#include <map>
#include <string>

#include "PDAnalysis/Ntu/interface/PDNtupleData.h"
#include "PDAnalysis/Ntu/interface/PDEnumString.h"

class PDGenHandler: public virtual PDNtupleData {

 public:

  PDGenHandler();
  ~PDGenHandler();

  int chkId( int id );

  const std::vector<ObjectAssociation>& genAssociation(
        ObjectDistance* dist, PDEnumString::recoObject obj,
        ObjectSelection* sel = 0,
        bool cache = false );

  static void print( const std::vector<int>* genId,
                     const std::vector<int>* genMother );

  void reshuffleFw( std::vector<int>* newMother,
		    std::vector<int>* newPartner );
  void reshuffleBw( std::vector<int>* newMother,
                    std::vector<int>* newPartner );

  const std::vector<int>& primaryHadrons();
  const std::vector<int>& finalPartons();
  const std::set<int>& primaryHadrSet();
  const std::set<int>& finalPartSet();

  const std::vector<int>& allMothers      ( int iGen );
  const std::vector<int>& allDaughters    ( int iGen );
  const std::vector<int>& finalDescendants( int iGen );

  int hasOscillated( int iGen );
  int willOscillate( int iGen );

  int ancestors( int iGen );
  class TypeSelect {
    friend PDGenHandler;
   public:
    TypeSelect() { sel = 0; m = '\0'; reset = true; }
    TypeSelect( int id, char mode = 't' ) { sel = id; m = mode; reset = true; }
    virtual ~TypeSelect() {}
    virtual bool operator()( int id ) const {
      switch ( m ) {
      case 't': return genId->at( id ) == sel;
      case 'p': return id == sel;
      }
      return false;
    }
    virtual int dmax() const { return -1; }
    static std::vector<int>* genId;
   private:
    int sel;
    char m;
    mutable bool reset;
  };
  std::multimap<int,int>& findAncestor( int iGen, const TypeSelect& sel );
  bool hasAncestors( int iGen, int ns, const TypeSelect** sel );
  std::multimap<int,int>& findDescendant( int iGen, const TypeSelect& sel );
  bool hasDescendant( int iGen, int ns, const TypeSelect** sel );

  int  sameMother  ( int iGen );
  int  sameDaughter( int iGen );
  int firstMother  ( int iGen );
  int  lastDaughter( int iGen );

/*
  static const TypeSelect& isLepton       ;
  static const TypeSelect& isMuon         ;
  static const TypeSelect& isElectron     ;
  static const TypeSelect& isW            ;
  static const TypeSelect& isDirectW      ;
  static const TypeSelect& isTau          ;
  static const TypeSelect& isDirectTau    ;
  static const TypeSelect& isBHadron      ;
  static const TypeSelect& isDirectBHadron;
  static const TypeSelect& isBQuark       ;
  static const TypeSelect& isCHadron      ;
  static const TypeSelect& isDirectCHadron;
  static const TypeSelect& isCQuark       ;
  static const TypeSelect& isLHadron      ;
  static const TypeSelect& isDirectLHadron;
  static const TypeSelect& isLQuark       ;
  static const TypeSelect& isQuarkNoBeauty;
  static const TypeSelect& isQuarkNoTop   ;
  static const TypeSelect& isQuark        ;
  static const TypeSelect& isTQuark       ;
  static const TypeSelect& isDirectTQuark ;
  static const TypeSelect& isNotTop       ;
  static const TypeSelect& isGluon        ;
*/
 private:

  void associateMuon(     ObjectDistance* dist );
  void associateElectron( ObjectDistance* dist );
  void associateTau(      ObjectDistance* dist );
  void associateJet(      ObjectDistance* dist,
                          ObjectSelection* sel );

  static void print( const std::string& head, int iPart, int level = 0 );

  void newEvent();
  void addDaughters( int iGen, std::set<int>& dSet ); 
  void fillPartonMothers( int iHadr, std::vector<int>& hadMothers );

  unsigned int lastRun;
  unsigned int lastLumi;
  unsigned int lastEvent;

  std::map< ObjectDistance*, 
            std::map< PDEnumString::recoObject,
                      std::vector<ObjectAssociation> > > associationCache;

  std::map< int, std::vector<int> >   mothersCache;
  std::map< int, std::vector<int> > daughtersCache;
  std::map< int, std::vector<int> >     finalCache;
  std::map< const TypeSelect*, std::map< int, std::multimap<int,int> > >
      ancestorsCache;
  std::map< const TypeSelect*, std::map< int, std::multimap<int,int> > >
    descendantsCache;
  std::vector<int> hadronsCache;
  std::vector<int> partonsCache;
  std::set<int> hadronsSet;
  std::set<int> partonsSet;

  const static std::vector<int>* ids;
  const static std::vector<int>* mothers;
  static std::map< int, std::vector<int> > daughters;

  PDGenHandler           ( const PDGenHandler& a );
  PDGenHandler& operator=( const PDGenHandler& a );

};


class IsLepton: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    if ( iGen < 0 ) return false;
    int id = genId->at( iGen );
    if ( id ==  11 ) return true;
    if ( id == -11 ) return true;
    if ( id ==  13 ) return true;
    if ( id == -13 ) return true;
    return false;
  }
};

class IsMuon: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    if ( iGen < 0 ) return false;
    int id = genId->at( iGen );
    if ( id ==  13 ) return true;
    if ( id == -13 ) return true;
    return false;
  }
};

class IsElectron: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    if ( iGen < 0 ) return false;
    int id = genId->at( iGen );
    if ( id ==  11 ) return true;
    if ( id == -11 ) return true;
    return false;
  }
};

class IsW: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    if ( iGen < 0 ) return false;
    int id = genId->at( iGen );
    if ( id ==  24 ) return true;
    if ( id == -24 ) return true;
    return false;
  }
};

class IsDirectW: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    if ( iGen < 0 ) return false;
    int id = genId->at( iGen );
    if ( id ==  24 ) return true;
    if ( id == -24 ) return true;
    return false;
  }
  virtual int dmax() const { return 1; }
};

class IsTau: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    if ( iGen < 0 ) return false;
    int id = genId->at( iGen );
    if ( id ==  15 ) return true;
    if ( id == -15 ) return true;
    return false;
  }
};

class IsDirectTau: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    if ( iGen < 0 ) return false;
    int id = genId->at( iGen );
    if ( id ==  15 ) return true;
    if ( id == -15 ) return true;
    return false;
  }
  virtual int dmax() const { return 1; }
};

class IsHadron: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    if ( iGen < 0 ) return false;
    int id = genId->at( iGen );
    if ( id < 0 ) id = -id;
    id %= 10000;
    id /= 100;
    if ( ( id >=1 ) && ( id <= 5 ) ) return true;
    id /= 10;
    if ( ( id >=1 ) && ( id <= 5 ) ) return true;
    return false;
  }
};

class IsBHadron: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    if ( iGen < 0 ) return false;
    int id = genId->at( iGen );
    if ( id < 0 ) id = -id;
    id %= 10000;
    if ( ( id / 100  ) == 5 ) return true;
    if (   id < 1000 )       return false;
    if ( ( id / 1000 ) == 5 ) return true;
    return false;
  }
};

class IsDirectBHadron: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    static IsBHadron ibh;
    return ibh( iGen );
  }
  virtual int dmax() const { return 1; }
};

class IsBQuark: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    if ( iGen < 0 ) return false;
    int id = genId->at( iGen );
    if ( id ==  5 ) return true;
    if ( id == -5 ) return true;
    return false;
  }
};

class IsCHadron: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    if ( iGen < 0 ) return false;
    int id = genId->at( iGen );
    if ( id < 0 ) id = -id;
    id %= 10000;
    if ( ( id / 100  ) == 4 ) return true;
    if (   id < 1000 )       return false;
    if ( ( id / 1000 ) == 4 ) return true;
    return false;
  }
};

class IsDirectCHadron: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    static IsCHadron ich;
    return ich( iGen );
  }
  virtual int dmax() const { return 1; }
};

class IsCharmonium: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    if ( iGen < 0 ) return false;
    int id = genId->at( iGen );
    if ( id < 0 ) id = -id;
    id /= 10;
    return ( ( id % 100 ) == 44 );
  }
};

class IsDirectCharmonium: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    static IsCharmonium icc;
    return icc( iGen );
  }
  virtual int dmax() const { return 1; }
};

class IsCHNonCharmonium: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    static IsCharmonium icc;
    static IsCHadron ich;
    return ich( iGen ) && !icc( iGen );
  }
};

class IsDirectCHNonCharmonium: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    static IsCHNonCharmonium inr;
    return inr( iGen );
  }
  virtual int dmax() const { return 1; }
};

class IsCQuark: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    if ( iGen < 0 ) return false;
    int id = genId->at( iGen );
    if ( id ==  4 ) return true;
    if ( id == -4 ) return true;
    return false;
  }
};

class IsLHadron: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    if ( iGen < 0 ) return false;
    int id = genId->at( iGen );
    if ( id < 0 ) id = -id;
    if ( id < 100 ) return false;
    id %= 10000;
    if ( ( id / 100  ) < 4 ) return true;
    if (   id < 1000 )       return false;
    if ( ( id / 1000 ) < 4 ) return true;
    return false;
  }
};

class IsDirectLHadron: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    static IsLHadron ilh;
    return ilh( iGen );
    return false;
  }
  virtual int dmax() const { return 1; }
};

class IsLQuark: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    if ( iGen < 0 ) return false;
    int id = genId->at( iGen );
    if ( id < -3 ) return false;
    if ( id >  3 ) return false;
    return true;
  }
};

class IsQuarkNoBeauty: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    if ( iGen < 0 ) return false;
    int id = genId->at( iGen );
    if ( id < -4 ) return false;
    if ( id >  4 ) return false;
    return true;
  }
};

class IsQuarkNoTop: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    if ( iGen < 0 ) return false;
    int id = genId->at( iGen );
    if ( id < -5 ) return false;
    if ( id >  5 ) return false;
    return true;
  }
};

class IsQuark: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    if ( iGen < 0 ) return false;
    int id = genId->at( iGen );
    if ( id < -6 ) return false;
    if ( id >  6 ) return false;
    return true;
  }
};

class IsTQuark: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    if ( iGen < 0 ) return false;
    int id = genId->at( iGen );
    if ( id ==  6 ) return true;
    if ( id == -6 ) return true;
    return false;
  }
};

class IsDirectTQuark: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    if ( iGen < 0 ) return false;
    int id = genId->at( iGen );
    if ( id ==  6 ) return true;
    if ( id == -6 ) return true;
    return false;
  }
  virtual int dmax() const { return 1; }
};

class IsNotTop: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    if ( iGen < 0 ) return false;
    int id = genId->at( iGen );
    if ( id ==  6 ) return false;
    if ( id == -6 ) return false;
    return true;
  }
};

class IsGluon: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    if ( iGen < 0 ) return false;
    int id = genId->at( iGen );
    if ( id != 21 ) return false;
    return true;
  }
};

class IsString: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    if ( iGen < 0 ) return false;
    int id = genId->at( iGen );
    if ( id != 92 ) return false;
    return true;
  }
};

class IsCluster: public PDGenHandler::TypeSelect {
 public:
  virtual bool operator()( int iGen ) const {
    if ( iGen < 0 ) return false;
    int id = genId->at( iGen );
    if ( id != 91 ) return false;
    return true;
  }
};

extern const IsLepton        isLepton       ;
extern const IsMuon          isMuon         ;
extern const IsElectron      isElectron     ;
extern const IsW             isW            ;
extern const IsDirectW       isDirectW      ;
extern const IsTau           isTau          ;
extern const IsDirectTau     isDirectTau    ;
extern const IsHadron        isHadron       ;
extern const IsBHadron       isBHadron      ;
extern const IsDirectBHadron isDirectBHadron;
extern const IsBQuark        isBQuark       ;
extern const IsCHadron       isCHadron      ;
extern const IsDirectCHadron isDirectCHadron;
extern const IsCharmonium    isCharmonium   ;
extern const IsDirectCharmonium isDirectCharmonium;
extern const IsCHNonCharmonium isCHNonCharmonium;
extern const IsDirectCHNonCharmonium isDirectCHNonCharmonium;
extern const IsCQuark        isCQuark       ;
extern const IsLHadron       isLHadron      ;
extern const IsDirectLHadron isDirectLHadron;
extern const IsLQuark        isLQuark       ;
extern const IsQuarkNoBeauty isQuarkNoBeauty;
extern const IsQuarkNoTop    isQuarkNoTop   ;
extern const IsQuark         isQuark        ;
extern const IsTQuark        isTQuark       ;
extern const IsDirectTQuark  isDirectTQuark ;
extern const IsNotTop        isNotTop       ;
extern const IsGluon         isGluon        ;
extern const IsString        isString       ;
extern const IsCluster       isCluster      ;

#endif

