#ifndef NtuLightReader_H
#define NtuLightReader_H

#include "TROOT.h"
#include "TTree.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

class NtuLightReader {

 public:

  NtuLightReader() {
    currentTree = 0;
    analyzedEvts = acceptedEvts = 0;
    dum = "";
  };
  virtual ~NtuLightReader() {};

  TTree* currentTree;
  std::string currentFile;
  std::string treeName;

  virtual void initTree() = 0;
  void openTree( TTree* tree ) {
    currentTree = tree;
    initTree();
    return;
  }

  virtual bool analyze( int entry, int event_file, int event_tot ) {
    return true;
  }

  int analyzedEvents() { return analyzedEvts; }
  int acceptedEvents() { return acceptedEvts; }
  virtual void plot() { return; }

  virtual void setConfiguration( const std::string& file ) {
    std::ifstream cfg( file.c_str() );
    int length;
    int lenMax = 1000;
    char* line = new char[lenMax];
    char* lptr;
    while ( cfg.getline( line, lenMax ) ) {
      lptr = line;
      while ( *lptr == ' ' ) ++lptr;
      if    ( *lptr == '#' ) continue;
      std::string key( lptr );
      length = key.find( " " );
      if ( length < 0 ) {
        std::cout << "invalid configuration input: " << line << std::endl;
        continue;
      }
      key = key.substr( 0, length );
      lptr += length;
      while ( *lptr == ' ' ) ++lptr;
      std::string val( lptr );
      length = val.find( " " );
      if ( length >= 0 )
	val = val.substr( 0, length );
      setUserParameter( key, val );
    }
    delete[] line;
    return;
  }

  void setUserParameter( const std::string& key,
                         const std::string& val ) {
    std::map<std::string,std::string>::iterator
                                       iter = userParameters.find( key );
    std::map<std::string,std::string>::iterator
                                       iend = userParameters.end();
    if ( iter != iend ) iter->second = val;
    else                userParameters.insert( make_pair( key, val ) );
    return;
  }

  const
  std::string& getUserParameter( const std::string& key ) {
    std::map<std::string,std::string>::iterator
                                       iter = userParameters.find( key );
    std::map<std::string,std::string>::iterator
                                       iend = userParameters.end();
    if ( iter != iend ) return iter->second;
    return dum;
  }


  template<class T>
  void getUserParameter( const std::string& key, T& val ) {
    std::stringstream sstr;
    sstr.str( getUserParameter( key ) );
    sstr >> val;
    return;
  }

  void getUserParameter( const std::string& key, bool& val ) {
    const char* flag = getUserParameter( key ).c_str();
    val = (   ( *flag == 't' ) || ( *flag == 'T' ) ||
            ( ( *flag >= '1' ) && ( *flag <= '9' ) ) );
    return;
  }

 protected:

  class AutoSavedObject {
   public:
    typedef std::vector<const TObject*> obj_list;
    typedef obj_list::const_iterator    obj_iter;
    typedef std::map<const TObject*, TDirectory*> dir_map;
    typedef dir_map::const_iterator               dir_iter;
    void insert( const TObject* obj, TDirectory* dir ) {
      dir_iter iter = directoryMap.find( obj );
      dir_iter iend = directoryMap.end();
      if ( iter != iend ) return;
      directoryMap.insert( std::pair<const TObject*,
                                     TDirectory*>( obj, dir ) );
      return;
    }
    obj_iter objBegin()                    { return   objectList.begin();     }
    obj_iter objEnd()                      { return   objectList.end();       }
    dir_iter dirBegin()                    { return directoryMap.begin();     }
    dir_iter dirFind( const TObject* obj ) { return directoryMap.find( obj ); }
    dir_iter dirEnd()                      { return directoryMap.end();       }
    AutoSavedObject& operator=( TObject* obj ) {
      objectList.push_back( obj );
      std::string type = obj->ClassName();
      if ( type.substr( 0, 10 ) == "TDirectory" ) 
           dynamic_cast<TDirectory*>( obj )->cd();
      else if ( obj == gROOT )
           gROOT->cd();
      return *this;
    }
    template<class T>
    AutoSavedObject& operator=( const std::vector<T*>& vObj ) {
      int i;
      int n = vObj.size();
      for ( i = 0; i < n; ++i ) *this = vObj[i];
      return *this;
    }
    template<class T>
    AutoSavedObject& operator=( const std::vector<T*>* vObj ) {
      *this = *vObj;
      return *this;
    }
   private:
    obj_list objectList;
    dir_map directoryMap;
  };
  AutoSavedObject autoSavedObject;
  virtual bool writable( const TObject* obj ) {
    return writable( obj->ClassName() );
  }
  virtual bool writable( const std::string& type ) {
    if ( type.substr( 0, 2 ) == "TH"       ) return true;
    if ( type                == "TProfile" ) return true;
    if ( type                == "TCanvas"  ) return true;
    return false;
  }
  void autoSave( TList* list = 0 ) {
    std::string type = gDirectory->ClassName();
    if ( ( type != "TFile" ) && 
       ( type != "TDirectoryFile" ) ) {
      std::cout << "File not open" << std::endl;
      return;
    }
    autoSavedObject.insert( gROOT, gDirectory );
    AutoSavedObject::obj_iter o_iter = autoSavedObject.objBegin();
    AutoSavedObject::obj_iter o_iend = autoSavedObject.objEnd();
    while ( o_iter != o_iend ) {
      const TObject* obj = *o_iter++;
      type = obj->ClassName();
      AutoSavedObject::dir_iter d_iter = autoSavedObject.dirFind( obj );
      AutoSavedObject::dir_iter d_iend = autoSavedObject.dirEnd();
      bool dirFound = ( d_iter != d_iend );
      bool dirType  = ( type.substr( 0, 10 ) == "TDirectory" );
      if ( dirType || dirFound ) {
        TDirectory* dir;
        if ( !dirFound ) {
          std::string name = obj->GetName();
          dir = gDirectory->mkdir( name.c_str() );
          autoSavedObject.insert( obj, dir );
        }
        else {
          dir = d_iter->second;
        }
        dir->cd();
      }
      else if ( writable( obj ) ) {
        if ( list == 0 ) obj->Write();
        else             list->FindObject( obj->GetName() )->Write();
      }
    }
    return;
  }

  const char* dummyString;
  std::string dum;
  int dummySize;
  int dummyLevel;

  int analyzedEvts;
  int acceptedEvts;

  void setBranch( const char* branchName, void* dataPtr,
                  const char* branchData ) {
    if ( currentTree == 0 ) return;
    currentTree->SetBranchAddress( branchName, dataPtr  );
    dummyString = branchData;
    return;
  }
  void setBranch( const char* branchName, void* dataPtr,
                  const char* branchData,
                  TBranch** branchPtr ) {
    if ( currentTree == 0 ) return;
    currentTree->SetBranchAddress( branchName, dataPtr, branchPtr );
    dummyString = branchData;
    return;
  }

  template<class T>
  void setBranch( const char* branchName, T** dataPtr,
                  int bufferSize = 32000,
                  int splitLevel = 99 ) {
    if ( currentTree == 0 ) return;
    currentTree->SetBranchAddress( branchName, dataPtr  );
    dummySize  = bufferSize;
    dummyLevel = splitLevel;
    return;
  }
  template<class T>
  void setBranch( const char* branchName, T** dataPtr,
                  int bufferSize,
                  int splitLevel,
                  TBranch** branchPtr ) {
    if ( currentTree == 0 ) return;
    currentTree->SetBranchAddress( branchName, dataPtr, branchPtr );
    dummySize  = bufferSize;
    dummyLevel = splitLevel;
    return;
  }

  bool edmNtuple;
  virtual void process( TBranch* b, int ientry ) {}

 private:

  // dummy copy constructor and assignment
  NtuLightReader           ( const NtuLightReader& );
  NtuLightReader& operator=( const NtuLightReader& );

  std::map<std::string,std::string> userParameters;

};

#endif

