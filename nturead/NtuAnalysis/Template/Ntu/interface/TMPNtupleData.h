#ifndef TMPNtupleData_H
#define TMPNtupleData_H

#include <string>
#include <vector>

#include "NtuAnalysis/Common/interface/NtuData.h"

class TMPNtupleData: public NtuData<float> {

 public:

  TMPNtupleData();
  virtual ~TMPNtupleData();

 protected:

  // Declaration of leaf types

  // muons
  bool use_muons;
  int nMuons;
  std::vector<number>* muoPt;

  // jets
  bool use_jets;
  int nJets;
  std::vector<number>* jetPt;

 private:

  // dummy copy constructor and assignment
  TMPNtupleData           ( const TMPNtupleData& td );
  TMPNtupleData& operator=( const TMPNtupleData& td );

};

#endif

