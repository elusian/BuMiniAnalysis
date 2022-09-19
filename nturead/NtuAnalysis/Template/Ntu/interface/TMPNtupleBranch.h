#ifndef TMPNtupleBranch_H
#define TMPNtupleBranch_H

#include "NtuAnalysis/Common/interface/NtuEventBranch.h"
#include "TMPAnalysis/Ntu/interface/TMPNtupleData.h"
#include <string>
#include <vector>

class TBranch;

template <class T>
class TMPNtupleBranch: public virtual TMPNtupleData, public virtual T,
                      public NtuEventBranch<T> {

 public:

  TMPNtupleBranch();
  virtual ~TMPNtupleBranch();

  void initTree();

 protected:

  // List of branches

  // muons
  void setBranches_muons();
  TBranch* b_nMuons;
  TBranch* b_muoPt;

  // jets
  void setBranches_jets();
  TBranch* b_nJets;
  TBranch* b_jetPt;

 private:

  // dummy copy constructor and assignment
  TMPNtupleBranch           ( const TMPNtupleBranch& td );
  TMPNtupleBranch& operator=( const TMPNtupleBranch& td );

};

#include "TMPNtupleBranch.hpp"

#endif

