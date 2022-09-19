#include <iostream>

template <class T>
TMPNtupleBranch<T>::TMPNtupleBranch() {
  T::treeName = "TMPtree";
}


template <class T>
TMPNtupleBranch<T>::~TMPNtupleBranch() {
}


template <class T>
void TMPNtupleBranch<T>::initTree() {

  std::cout << "TMPNtupleBranch<T>::initTree()" << std::endl;

  // define branches
  // use local function "setBranch" (defined in the utility)
  // in place of ROOT "Branch" & "SetBranchAddress"

  if ( use_muons ) setBranches_muons();
  if ( use_jets  ) setBranches_jets ();

  return;

}


template <class T>
void TMPNtupleBranch<T>::setBranches_muons() {
  // muons
  int size = 8192;
  this->setBranch( "nMuons"         , &  nMuons         ,
                   "nMuons/I"       , &b_nMuons          );
  this->setBranch( "muoPt"          , &  muoPt          ,
                   size, 99         , &b_muoPt           );
  return;
}


template <class T>
void TMPNtupleBranch<T>::setBranches_jets() {
  // jets
  int size = 8192;
  this->setBranch( "nJets"          , &  nJets          ,
                   "nJets/I"        , &b_nJets           );
  this->setBranch( "jetPt"          , &  jetPt          ,
                   size, 99         , &b_jetPt           );
  return;
}

