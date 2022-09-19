#include <iostream>
#include <sstream>
#include <string>
#include <math.h>

#include "NtuAnalysis/Common/interface/NtuAnalyzerUtil.h"

#include "TDirectory.h"
#include "TBranch.h"
#include "TTree.h"

NtuAnalyzerUtil::NtuAnalyzerUtil() {
  currentRun     = 0;
  currentEvent   = 0;
  currentEvBase  = 0;
  currentEvSetup = 0;
}


NtuAnalyzerUtil::~NtuAnalyzerUtil() {
}


void NtuAnalyzerUtil::beginJob() {
  return;
}


void NtuAnalyzerUtil::getHeader( int ientry ) {
//  if ( currentEvBase != 0 ) return;
//  if ( currentEvent  != 0 ) return;
  if ( this->b_runNumber   != 0 ) this->b_runNumber  ->GetEntry( ientry );
  if ( this->b_lumiSection != 0 ) this->b_lumiSection->GetEntry( ientry );
  if ( this->b_eventNumber != 0 ) this->b_eventNumber->GetEntry( ientry );
  return;
}

