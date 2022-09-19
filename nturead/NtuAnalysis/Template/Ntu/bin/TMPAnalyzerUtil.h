#ifndef TMPAnalyzerUtil_H
#define TMPAnalyzerUtil_H

#include "TMPAnalysis/Ntu/interface/TMPNtupleBranch.h"
#include "NtuAnalysis/Common/interface/NtuInterface.h"

class TMPAnalyzerUtil: public NtuInterface< TMPNtupleBranch<WrapperBase> > {

 public:

  TMPAnalyzerUtil();
  virtual ~TMPAnalyzerUtil();

 protected:

  void setupNtuple();
  virtual void beginJob();

  // pre-selection, with minimal data process 
  // before full ntuple entry reading
  virtual bool preSelect( int ientry );

 private:

  // dummy copy constructor and assignment
  TMPAnalyzerUtil           ( const TMPAnalyzerUtil& );
  TMPAnalyzerUtil& operator=( const TMPAnalyzerUtil& );

};


#endif

