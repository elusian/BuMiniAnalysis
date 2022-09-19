#include <iostream>

using namespace std;

#include "TMPAnalysis/Ntu/interface/TMPNtupleData.h"

TMPNtupleData::TMPNtupleData() {

  // muons
  nMuons = 0;
  muoPt = new vector<number>;

  // jets
  nJets = 0;
  jetPt = new vector<number>;

}

TMPNtupleData::~TMPNtupleData() {
}


