#include <iostream>
#include <sstream>
#include <string>
#include <math.h>

#include "PDAnalyzer.h"

#include "TDirectory.h"
#include "TBranch.h"
#include "TTree.h"
#include "TCanvas.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TMatrix.h"
#include "TVector3.h"
#include "TVector.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "PDSecondNtupleWriter.h"

using namespace std;

PDAnalyzer::PDAnalyzer() {

  std::cout << "new PDAnalyzer" << std::endl;

  // user parameters are set as names associated to a string, 
  // default values can be set in the analyzer class contructor

  setUserParameter( "verbose", "f" );

}


PDAnalyzer::~PDAnalyzer() {
}



void PDAnalyzer::beginJob() {

  PDAnalyzerUtil::beginJob();

  // user parameters are retrieved as strings by using their names;
  // numeric parameters ( int, float or whatever ) can be directly set
  // by passing the corresponding variable,
  // e.g. getUserParameter( "name", x )

  getUserParameter( "verbose", verbose );
  
  tWriter = new PDSecondNtupleWriter; // second ntuple
  tWriter->open( getUserParameter("outputFile"), "RECREATE" ); // second ntuple
  
  decayOut = new ofstream("decays");

  return;

}


void PDAnalyzer::book() {

  // putting "autoSavedObject" in front of the histo creation 
  // it's automatically marked for saving on file; the option 
  // is uneffective when not using the full utility

  autoSavedObject =
  hptmumax        = new TH1D( "hptmumax"    , "ptmumax"    ,  50, 0.0, 100.0 );
  autoSavedObject =
  hptmuall        = new TH1D( "hptmuall"    , "ptmuall"    ,  50, 0.0, 100.0 );
  
  

  return;

}


void PDAnalyzer::reset() {
// automatic reset
  autoReset();
  return;
}

const float MUMASS  = 0.105658;
const float KMASS   = 0.493677;
const float BSMASS  = 5.36689;
const float PIMASS  = 0.139570;
const float B0MASS  = 5.27963;

static size_t nGenBu = 0;
static size_t nBu = 0;
static size_t nMatchedBuRecoToSim = 0;
static size_t nMatchedBuSimToReco = 0;
static size_t nDupBu = 0;
static size_t nUnmatchedBuForTrkReasons = 0;



bool PDAnalyzer::analyze( int entry, int event_file, int event_tot ) {

  if ( verbose ) {
    clog << " +++++++++++++++++++++++++++ " << endl;
    clog << "entry: "
         << entry << " " << event_file << " " << event_tot << endl;
    clog << "run: " <<   runNumber << " , "
         << "evt: " << eventNumber << endl;
  }
  else {
//    if ( !( event_file % 10000 ) || !( event_tot % 10000 ) )
    if ( !( event_tot % 10000 ) && event_tot )
      clog << event_file << " " << event_tot << endl;
  }
  tWriter->Reset();
  
  //ofstream buOut(Form("bu%d.log", eventNumber));
  //ofstream jPsiOut(Form("jPsi%d.log", eventNumber));
  //ofstream psi2SOut(Form("psi2S%d.log", eventNumber));
  
  convSpheCart(muoPt, muoEta, muoPhi, muoPx, muoPy, muoPz);
  convSpheCart(trkPt, trkEta, trkPhi, trkPx, trkPy, trkPz);
  convSpheCart(pfcPt, pfcEta, pfcPhi, pfcPx, pfcPy, pfcPz);
  
  tWriter->evtNumber = eventNumber;
  
  map <int, int> vertexMap;
  size_t saved = 0;
  for (int iSvt = 0; iSvt < nSVertices; iSvt++) {
  
    auto type = svtType->at(iSvt);
    if (type != PDEnumString::svtBuPsi2SK and type != PDEnumString::svtPsi2S and type != PDEnumString::svtJPsi) continue;
    
    vertexMap[iSvt] = saved;
    saved++;
  }
  
  vector<int> genB;
  vector<int> genBDau;
  for( uint i = 0; i < genId->size(); ++i ){
    int id = abs(genId->at(i));
    if( id != 521 ) continue;
    
    auto vertexComp = getBuVertexComponents(i);
    if (vertexComp) {
      genB.push_back(i);
      for (int iDau = 0; iDau < 5; iDau++) {
        genBDau.push_back((*vertexComp)[iDau]);
      }
    }
    
  }
  
  vector<int> bus;
  
  for (int iSvt = 0; iSvt < nSVertices; iSvt++) {
  
    auto type = svtType->at(iSvt);
    if (type != PDEnumString::svtBuPsi2SK and type != PDEnumString::svtPsi2S and type != PDEnumString::svtJPsi) continue;
    
    ROOT::Math::XYZVector svtP;
    
    auto&& tkSvt = tracksFromSV(iSvt);
    
    for (auto&& iSvtTrk: tkSvt) {
      ROOT::Math::RhoEtaPhiVector trkP(trkPt->at(iSvtTrk), trkEta->at(iSvtTrk), trkPhi->at(iSvtTrk));
      
      svtP += trkP;
    }
    
    tWriter->svtPt->push_back(svtP.Rho());
    tWriter->svtEta->push_back(svtP.Eta());
    tWriter->svtPhi->push_back(svtP.Phi());
    tWriter->svtMass->push_back(svtMass->at(iSvt));
    tWriter->svtChi2->push_back(svtChi2->at(iSvt));
    tWriter->svtNdof->push_back(svtNDOF->at(iSvt));
    tWriter->svtDist2D->push_back(svtDist2D->at(iSvt));
    tWriter->svtSigma2D->push_back(svtSigma2D->at(iSvt));
    tWriter->svtDist3D->push_back(svtDist3D->at(iSvt));
    tWriter->svtSigma3D->push_back(svtSigma3D->at(iSvt));
    tWriter->svtType->push_back(svtType->at(iSvt));
    tWriter->svtIsGoodMatch->push_back(0);
    
    tWriter->svtGenId->push_back(0);
    tWriter->svtDecayId->push_back(-1);
    
    if (type == PDEnumString::svtBuPsi2SK) {
//      clog << subVtxFromSV(iSvt).size() << endl;
//      for (unsigned int i = 0; i < subVtxFromSV(iSvt).size(); i++) {
//        int iPsi2S = (subVtxFromSV(iSvt)).at(i);
//        clog << svtType->at(iPsi2S) << " ";
//      }
//      clog << endl;
       bus.push_back(iSvt);
    
      int iPsi2S = (subVtxFromSV(iSvt)).at(0);
      if(svtType->at(iPsi2S) != PDEnumString::svtPsi2S){
          clog<<"(subVtxFromSV(iSvt)).at(0) NOT A PSI2S but "<<svtType->at(iPsi2S)<<endl;
          return false;
      }
      tWriter->svtPsi2SIndex->push_back(vertexMap[iPsi2S]);
      
      int iJPsi  = (subVtxFromSV(iPsi2S)).at(0);
      if(svtType->at(iJPsi) != PDEnumString::svtJPsi){
          clog<<"(subVtxFromSV(iPsi2S)).at(0) NOT A JPSI but "<<svtType->at(iJPsi)<<endl;
          return false;
      }
      tWriter->svtJPsiIndex->push_back(vertexMap[iJPsi]);
      
      const vector<int>& mu_idx = tracksFromSV(iJPsi);
      vector<int> pi_idx  = tracksFromSV(iPsi2S);
      vector<int> k_idx = tracksFromSV(iSvt);
      
      pi_idx.erase(remove_if(pi_idx.begin(), pi_idx.end(), [&mu_idx](int trkId){return trkId == mu_idx[0] or trkId == mu_idx[1];}), pi_idx.end());
      if (pi_idx.size() != 2) {
        throw std::runtime_error(Form("pi collection size is not 2 but %zd", pi_idx.size()));
      }
      
      k_idx.erase(remove_if(k_idx.begin(), k_idx.end(), [&pi_idx, &mu_idx](int trkId){return trkId == mu_idx[0] or trkId == mu_idx[1] or trkId == pi_idx[0] or trkId == pi_idx[1];}), k_idx.end());
      if (k_idx.size() != 1) {
        throw std::runtime_error(Form("k collection size is not 1 but %zd", k_idx.size()));
      }
      
      tWriter->svtPi1Index->push_back(pi_idx[0]);
      tWriter->svtPi2Index->push_back(pi_idx[1]);
      tWriter->svtMu1Index->push_back(mu_idx[0]);
      tWriter->svtMu2Index->push_back(mu_idx[1]);
      tWriter->svtKIndex->push_back(k_idx[0]);
        
      auto iGen = GetClosestGen( svtP.Eta(), svtP.Phi(), svtP.Rho() );
        
      if (iGen >= 0) {
        tWriter->svtGenId->back() = abs(genId->at(iGen));
//        if (abs(genId->at(iGen)) != 521 and abs(genId->at(iGen)) != 523 ) {
//          clog << "Not a bu, it's a " << genId->at(iGen) << endl;
//        
//          printDecayChain(iGen);
//          clog << "Bu(" << iSvt <<") -> " << endl << 
//            "Psi2S(" << iPsi2S << ") + k(" << k_idx[0] << ") -> " << endl << 
//            "JPsi(" << iJPsi << ") + pi(" << pi_idx[0] << ") + pi(" << pi_idx[1] << ") + k(" << k_idx[0] << ") -> " <<  endl << 
//            "mu(" << mu_idx[0] << ") + mu(" << mu_idx[1] << ") + pi(" << pi_idx[0] << ") + pi(" << pi_idx[1] << ") + k(" << k_idx[0] << ")" << endl;
//        }
        stringstream decayStream;
        printCompactDecayChain(iGen, decayStream);
        
        auto decayStr = decayStream.str();
        if (decays.count(decayStr) == 0) {
          auto index = decays.size();
          decays[decayStr] = index;
          (*decayOut) << index << " " << decayStr << endl;
        }
        tWriter->svtDecayId->back() = decays.at(decayStr);
      }
      
      
      ROOT::Math::XYZVector svtPos(svtX->at(iSvt), svtY->at(iSvt), svtZ->at(iSvt));
      
      auto iPV = findPV(svtPos, svtP);
      
      if (iPV > 0) {
    
        ROOT::Math::XYZVector pvPos(pvtX->at(iPV),pvtY->at(iPV),pvtZ->at(iPV));
        
        auto BMASS = 5.27932;
        
        auto PVSVdistGV = svtPos - pvPos;
        
        TVector3 PVSVdist(PVSVdistGV.X(), PVSVdistGV.Y(), 0);
        
        TVector3 BsMomentum(svtP.X(), svtP.Y(), 0);
        
        tWriter->svtCt->push_back(BMASS*PVSVdist.Dot(BsMomentum)/BsMomentum.Mag2());
        
        TMatrixF covSV(3,3);
        float covSVArray[]={
          svtSxx->at(iSvt),svtSxy->at(iSvt),svtSxz->at(iSvt),
          svtSxy->at(iSvt),svtSyy->at(iSvt),svtSyz->at(iSvt),
          svtSxz->at(iSvt),svtSyz->at(iSvt),svtSzz->at(iSvt)
        };
        covSV.SetMatrixArray(covSVArray);
          
        TMatrixF covPV(3,3);
        float covPVArray[]={
          pvtSxx->at(iPV),pvtSxy->at(iPV),pvtSxz->at(iPV),
          pvtSxy->at(iPV),pvtSyy->at(iPV),pvtSyz->at(iPV),
          pvtSxz->at(iPV),pvtSyz->at(iPV),pvtSzz->at(iPV)
        };
        covPV.SetMatrixArray(covPVArray);
          
        TMatrixF covTot= covSV+covPV;
        
        double diff2DA[3] = {PVSVdistGV.X(), PVSVdistGV.Y(), 0};
        TVectorD diff2D(3, diff2DA);
        
        auto ctErr = TMath::Abs(BMASS/BsMomentum.Mag()*sqrt(covTot.Similarity(diff2D))/sqrt(diff2D.Norm2Sqr())*cos(PVSVdist.Angle(BsMomentum)));
        
        tWriter->svtCtErr->push_back(ctErr);
      }
      else {
        tWriter->svtCt->push_back(-1);
        tWriter->svtCtErr->push_back(-1);
      }
    }
    else if (type == PDEnumString::svtPsi2S) {
      tWriter->svtCt->push_back(-1);
      tWriter->svtCtErr->push_back(-1);
      tWriter->svtKIndex->push_back(-1);
      tWriter->svtPsi2SIndex->push_back(-1);
      
      int iJPsi  = (subVtxFromSV(iSvt)).at(0);
      if(svtType->at(iJPsi) != PDEnumString::svtJPsi){
          clog<<"(subVtxFromSV(iPsi2S)).at(0) NOT A JPSI but "<<svtType->at(iJPsi)<<endl;
          return false;
      }
      tWriter->svtJPsiIndex->push_back(vertexMap[iJPsi]);
      const vector<int>& mu_idx = tracksFromSV(iJPsi);
      vector<int> pi_idx  = tracksFromSV(iSvt);
      
      pi_idx.erase(remove_if(pi_idx.begin(), pi_idx.end(), [&mu_idx](int trkId){return trkId == mu_idx[0] or trkId == mu_idx[1];}), pi_idx.end());
      if (pi_idx.size() != 2) {
        clog << "pi collection size in psi2S is not 2 but " << pi_idx.size() << endl;
        return false;
      }
      
      tWriter->svtPi1Index->push_back(pi_idx[0]);
      tWriter->svtPi2Index->push_back(pi_idx[1]);
      tWriter->svtMu1Index->push_back(mu_idx[0]);
      tWriter->svtMu2Index->push_back(mu_idx[1]);
      
//      psi2SOut << min(trkPt->at(mu_idx[0]), trkPt->at(mu_idx[1])) << " " << max(trkPt->at(mu_idx[0]), trkPt->at(mu_idx[1])) << " "
//        << min(trkPt->at(pi_idx[0]), trkPt->at(pi_idx[1])) << " " << max(trkPt->at(pi_idx[0]), trkPt->at(pi_idx[1])) << " | ";
//      psi2SOut << svtP.Rho() << " " << svtP.Eta() << " " << svtP.Phi() << " "
//        << svtMass->at(iSvt) << " "
//        << TMath::Prob(svtChi2->at(iSvt), svtNDOF->at(iSvt)) << " | ";
//        
////      if (iGen >= 0) {
////        tWriter->svtGenId->back() = abs(genId->at(iGen));
////        printCompactDecayChain(iGen, psi2SOut);
////      }
////      else {
////        psi2SOut << "No gen available "
////      }
//      psi2SOut << endl;
    }
    else if (type == PDEnumString::svtJPsi) {
      tWriter->svtCt->push_back(-1);
      tWriter->svtCtErr->push_back(-1);
      tWriter->svtPi1Index->push_back(-1);
      tWriter->svtPi2Index->push_back(-1);
      tWriter->svtKIndex->push_back(-1);
      tWriter->svtJPsiIndex->push_back(-1);
      tWriter->svtPsi2SIndex->push_back(-1);
      
      const vector<int>& mu_idx = tracksFromSV(iSvt);
      
      tWriter->svtMu1Index->push_back(mu_idx[0]);
      tWriter->svtMu2Index->push_back(mu_idx[1]);
      
//      jPsiOut << min(trkPt->at(mu_idx[0]), trkPt->at(mu_idx[1])) << " " << max(trkPt->at(mu_idx[0]), trkPt->at(mu_idx[1])) << " | ";
//      jPsiOut << svtP.Rho() << " " << svtP.Eta() << " " << svtP.Phi() << " "
//        << svtMass->at(iSvt) << " "
//        << TMath::Prob(svtChi2->at(iSvt), svtNDOF->at(iSvt)) << " | ";
//        
////      if (iGen >= 0) {
////        tWriter->svtGenId->back() = abs(genId->at(iGen));
////        printCompactDecayChain(iGen, jPsiSOut);
////      }
////      else {
////        jPsiSOut << "No gen available "
////      }
//      jPsiOut << endl;
    }
    else {
      tWriter->svtCt->push_back(-1);
      tWriter->svtCtErr->push_back(-1);
      tWriter->svtPi1Index->push_back(-1);
      tWriter->svtPi2Index->push_back(-1);
      tWriter->svtMu1Index->push_back(-1);
      tWriter->svtMu2Index->push_back(-1);
      tWriter->svtKIndex->push_back(-1);
      tWriter->svtJPsiIndex->push_back(-1);
      tWriter->svtPsi2SIndex->push_back(-1);
    }
  }
  
  sort(bus.begin(), bus.end(), [this, &vertexMap](int iSvt1, int iSvt2){
    int iStoredSv1 = vertexMap[iSvt1];
    int iStoredSv2 = vertexMap[iSvt2];
    
    int mu1Index1 = this->tWriter->svtMu1Index->at(iStoredSv1);
    int mu1Index2 = this->tWriter->svtMu1Index->at(iStoredSv2);
    
    int mu2Index1 = this->tWriter->svtMu2Index->at(iStoredSv1);
    int mu2Index2 = this->tWriter->svtMu2Index->at(iStoredSv2);
    
    int muMaxIndex1 = mu1Index1;
    int muMinIndex1 = mu2Index1;
    int muMaxIndex2 = mu1Index2;
    int muMinIndex2 = mu2Index2;
    
    if (this->trkPt->at(mu1Index1) < this->trkPt->at(mu2Index1)) {
      muMaxIndex1 = mu2Index1;
      muMinIndex1 = mu1Index1;
    }
    if (this->trkPt->at(mu1Index2) < this->trkPt->at(mu2Index2)) {
      muMaxIndex2 = mu2Index2;
      muMinIndex2 = mu1Index2;
    }
    
    int pi1Index1 = this->tWriter->svtPi1Index->at(iStoredSv1);
    int pi1Index2 = this->tWriter->svtPi1Index->at(iStoredSv2);
    
    int pi2Index1 = this->tWriter->svtPi2Index->at(iStoredSv1);
    int pi2Index2 = this->tWriter->svtPi2Index->at(iStoredSv2);
    
    int piMaxIndex1 = pi1Index1;
    int piMinIndex1 = pi2Index1;
    int piMaxIndex2 = pi1Index2;
    int piMinIndex2 = pi2Index2;
    
    if (this->trkPt->at(pi1Index1) < this->trkPt->at(pi2Index1)) {
      piMaxIndex1 = pi2Index1;
      piMinIndex1 = pi1Index1;
    }
    if (this->trkPt->at(pi1Index2) < this->trkPt->at(pi2Index2)) {
      piMaxIndex2 = pi2Index2;
      piMinIndex2 = pi1Index2;
    }
    
    int kIndex1 = this->tWriter->svtKIndex->at(iStoredSv1);
    int kIndex2 = this->tWriter->svtKIndex->at(iStoredSv2);
    
    auto tuple1 = make_tuple(this->trkPt->at(muMinIndex1), this->trkPt->at(muMaxIndex1),this->trkPt->at(piMinIndex1),this->trkPt->at(piMaxIndex1),this->trkPt->at(kIndex1));
    auto tuple2 = make_tuple(this->trkPt->at(muMinIndex2), this->trkPt->at(muMaxIndex2),this->trkPt->at(piMinIndex2),this->trkPt->at(piMaxIndex2),this->trkPt->at(kIndex2));
    
    return tuple1 < tuple2;
    
  });
  
//  for (auto iBu: bus) {
//    auto iStoredBu = vertexMap[iBu];
//    
//    auto mu1 = tWriter->svtMu1Index->at(iStoredBu);
//    auto mu2 = tWriter->svtMu2Index->at(iStoredBu);
//    auto pi1 = tWriter->svtPi1Index->at(iStoredBu);
//    auto pi2 = tWriter->svtPi2Index->at(iStoredBu);
//    auto k = tWriter->svtKIndex->at(iStoredBu);
//    
//    buOut << min(trkPt->at(mu1), trkPt->at(mu2)) << " " << max(trkPt->at(mu1), trkPt->at(mu2)) << " "
//      << min(trkPt->at(pi1), trkPt->at(pi2)) << " " << max(trkPt->at(pi1), trkPt->at(pi2)) << " "
//      << trkPt->at(k) << " | ";
//    buOut << svtMass->at(iBu) << " "
//      << svtChi2->at(iBu) << " " << svtNDOF->at(iBu) << " "
//      << TMath::Prob(svtChi2->at(iBu), svtNDOF->at(iBu)) << " | ";
//      
//    auto iGen = GetClosestGen( tWriter->svtPt->at(iStoredBu), tWriter->svtEta->at(iStoredBu), tWriter->svtPhi->at(iStoredBu) );
//      
//    if (iGen >= 0) {
//      printCompactDecayChain(iGen, buOut);
//    }
//    else {
//      buOut << "No gen available ";
//    }
//    buOut << endl;
//  }

  vector<array<int, 6>> goodBus;
  vector<array<int, 5>> goodPsi2Ss;
  vector<array<int, 3>> goodJPsis;
  
  vector<int> goodBuMatch;
  
  for (int iGen = 0; iGen < nGenP; iGen++) {
    auto id = abs(genId->at(iGen));
    
    if (id != 521 and id != 443 and id != 100443) continue;
    
    tWriter->genPt->push_back(genPt->at(iGen));
    tWriter->genEta->push_back(genEta->at(iGen));
    tWriter->genPhi->push_back(genPhi->at(iGen));
    tWriter->genMass->push_back(genMass->at(iGen));
    tWriter->genCt->push_back( getCt( iGen ));
    tWriter->genCharge->push_back(genCharge->at(iGen));
    tWriter->genId->push_back(genId->at(iGen));
    
    if (id == 521) {
      auto decayTracks = getBuVertexComponents(iGen);
      if (decayTracks) {
        goodBus.push_back(*decayTracks);
        goodBuMatch.push_back(0);
      }
    }
    else if (id == 100443) {
      auto decayTracks = getPsi2SVertexComponents(iGen);
      if (decayTracks) {
        goodPsi2Ss.push_back(*decayTracks);
      }
    }
    else if (id == 443) {
      auto decayTracks = getJPsiVertexComponents(iGen);
      if (decayTracks) {
        goodJPsis.push_back(*decayTracks);
      }
    }
  }
  
  for (auto iGen: genBDau) {
    tWriter->genPt->push_back(genPt->at(iGen));
    tWriter->genEta->push_back(genEta->at(iGen));
    tWriter->genPhi->push_back(genPhi->at(iGen));
    tWriter->genMass->push_back(genMass->at(iGen));
    tWriter->genCt->push_back( getCt( iGen ));
    tWriter->genCharge->push_back(genCharge->at(iGen));
    tWriter->genId->push_back(genId->at(iGen));
  }
  
  vector<int> trackGenMatch(nTracks);
  unordered_map<int, int> genTrackMatch;
  
  for (int iTrk = 0; iTrk < nTracks; iTrk++) {
    int pv = trkPVtx->at(iTrk);
    if (pv >= nPVertices) pv -= nPVertices;
  
    tWriter->trkPt->push_back(trkPt->at(iTrk));
    tWriter->trkEta->push_back(trkEta->at(iTrk));
    tWriter->trkPhi->push_back(trkPhi->at(iTrk));
    tWriter->trkDxy->push_back(trkDxy->at(iTrk));
    tWriter->trkDz->push_back(dZ(iTrk, pvtX->at(pv), pvtY->at(pv), pvtZ->at(pv)));
    tWriter->trkExy->push_back(trkExy->at(iTrk));
    tWriter->trkEz->push_back(trkEz->at(iTrk));
    tWriter->trkCharge->push_back(trkCharge->at(iTrk));
    
    tWriter->trkPfcPt->push_back(0);
    tWriter->trkPfcEta->push_back(0);
    tWriter->trkPfcPhi->push_back(0);
    if (trkPFC->at(iTrk) >= 0) {
        tWriter->trkPfcPt->back() = pfcPt->at(trkPFC->at(iTrk));
        tWriter->trkPfcEta->back() = pfcEta->at(trkPFC->at(iTrk));
        tWriter->trkPfcPhi->back() = pfcPhi->at(trkPFC->at(iTrk));
    }
    
    auto hitPattern = trkHitPattern->at(iTrk);
    tWriter->trkNVHMuon->push_back(hitPattern / 1000000);
    tWriter->trkNVHPixel->push_back((hitPattern % 1000000) / 10000);
    tWriter->trkNVHTracker->push_back((hitPattern % 10000) / 100);
    tWriter->trkNVHAll->push_back(hitPattern % 100);
    
    auto layPattern = trkLayPattern->at(iTrk);
    tWriter->trkLPPixel->push_back(layPattern / 10000);
    tWriter->trkLPStrips->push_back((layPattern % 10000) / 100);
    tWriter->trkLPTracker->push_back(layPattern % 100);
    
    tWriter->trkMIH->push_back(trkMissingInnerHits->at(iTrk));
    
    tWriter->trkIsMu->push_back(trkType->at(iTrk) >= 1024);
    
    int trkKind = 2;
    auto trkPf = trkPFC->at(iTrk);
    if (trkPf >= 0) {
        trkKind = pfcHasTrackDetails->at(trkPf);
    }
    
    tWriter->trkKind->push_back(trkKind);
    
    tWriter->trkQOverPError->push_back(trkQOverPError->at(iTrk));
    tWriter->trkLambdaError->push_back(trkLambdaError->at(iTrk));
    tWriter->trkPhiError->push_back(trkPhiError->at(iTrk));
    tWriter->trkDszError->push_back(trkDszError->at(iTrk));
    tWriter->trkEtaError->push_back(trkEtaError->at(iTrk));
    tWriter->trkPtError->push_back(trkPtError->at(iTrk));
    
    tWriter->trkDxyDszCov->push_back(trkDxyDszCov->at(iTrk));
    tWriter->trkLambdaDszCov->push_back(trkLambdaDszCov->at(iTrk));
    tWriter->trkPhiDxyCov->push_back(trkPhiDxyCov->at(iTrk));
    
    tWriter->trkGenId->push_back(0);
    tWriter->trkGenPt->push_back(0);
    tWriter->trkGenEta->push_back(0);
    tWriter->trkGenPhi->push_back(0);
    
//    auto pxAtVtx = trkVtxPx->at(iTrk);
//    auto pyAtVtx = trkVtxPy->at(iTrk);
//    auto pzAtVtx = trkVtxPz->at(iTrk);
//    
//    auto ptAtVtx  = hypot(pxAtVtx, pyAtVtx);
//    auto etaAtVtx = asinh(pzAtVtx/ptAtVtx);
//    auto phiAtVtx = atan2(pyAtVtx, pxAtVtx);
    
    auto [iTrkGen, minDr] = GetClosestTrackGenFromColl(trkEta->at(iTrk), trkPhi->at(iTrk), trkPt->at(iTrk), trkCharge->at(iTrk), genBDau);
    //auto [iTrkGen, minDr] = GetClosestTrackGen(trkEta->at(iTrk), trkPhi->at(iTrk), trkPt->at(iTrk), trkCharge->at(iTrk));
    
    
    if (iTrkGen >= 0) {
      tWriter->trkGenId->back() = abs(genId->at(iTrkGen));
      tWriter->trkGenPt->back() = genPt->at(iTrkGen);
      tWriter->trkGenEta->back() = genEta->at(iTrkGen);
      tWriter->trkGenPhi->back() = genPhi->at(iTrkGen);
    }
    
    tWriter->trkGenMatchRadius->push_back(minDr);
    
    trackGenMatch[iTrk] = iTrkGen;
    genTrackMatch[iTrkGen] += 1;
  }
  
  auto crosscomp = [](const int* v1, const int* v2) -> bool {
    return (v1[0] == v2[0] and v1[1] == v2[1]) or (v1[0] == v2[1] and v1[1] == v2[0]);
  };
  
  for (size_t iSvt = 0; iSvt < tWriter->svtPt->size(); iSvt++) {
    if ( tWriter->svtType->at(iSvt) != 212) {
      continue;
    }
    
    nBu += 1;
    
    array<int, 5> svtTracks = {{
      tWriter->svtMu1Index->at(iSvt), 
      tWriter->svtMu2Index->at(iSvt), 
      tWriter->svtPi1Index->at(iSvt), 
      tWriter->svtPi2Index->at(iSvt), 
      tWriter->svtKIndex->at(iSvt)
    }};
    for (int i = 0; i < 5; i++) {
      svtTracks[i] = trackGenMatch[svtTracks[i]];
    }
    
    auto iPsi2SSV = tWriter->svtPsi2SIndex->at(iSvt);
    auto iJPsiSV = tWriter->svtJPsiIndex->at(iSvt);
    
    int isJPsiGood = 0;
    if (svtTracks[0] >= 0 and svtTracks[1] >= 0) {
      isJPsiGood = 1;
    }
    else {
      continue;
    }
    if (abs(genId->at(svtTracks[0])) == 13 and abs(genId->at(svtTracks[1])) == 13) {
      isJPsiGood = 2;
    }
    if (isJPsiGood == 2) {
      for (size_t iGoodJPsi = 0; iGoodJPsi < goodJPsis.size(); iGoodJPsi++) {
        if (crosscomp(svtTracks.data(), goodJPsis[iGoodJPsi].data())) {
          if (goodJPsis[iGoodJPsi][2] == 0) {
            isJPsiGood = 4;
          }
          else {
            isJPsiGood = 3;
          }
        }
      }
    }
    tWriter->svtIsGoodMatch->at(iJPsiSV) = isJPsiGood;
    if (isJPsiGood == 0) continue;
    
    int isPsi2SGood = 0;
    if (svtTracks[2] >= 0 and svtTracks[3] >= 0) {
      isPsi2SGood = 1;
    }
    else {
      continue;
    }
    if (isJPsiGood >= 2 and abs(genId->at(svtTracks[2])) == 211 and abs(genId->at(svtTracks[3])) == 211) {
      isPsi2SGood = 2;
    }
    if (isJPsiGood >= 3 and isPsi2SGood == 2) {
      for (size_t iGoodPsi2S = 0; iGoodPsi2S < goodPsi2Ss.size(); iGoodPsi2S++) {
        if (crosscomp(svtTracks.data(), goodPsi2Ss[iGoodPsi2S].data()) and crosscomp(svtTracks.data() + 2, goodPsi2Ss[iGoodPsi2S].data() + 2)) {
          if (goodPsi2Ss[iGoodPsi2S][4] == 0) {
            isPsi2SGood = 4;
          }
          else {
            isPsi2SGood = 3;
          }
        }
      }
    }
    tWriter->svtIsGoodMatch->at(iPsi2SSV) = isPsi2SGood;
    if (isPsi2SGood == 0) continue;
    
    int isBuGood = 0;
    if (svtTracks[4] >= 0) {
      isBuGood = 1;
    }
    else {
      continue;
    }
    if (isPsi2SGood >= 2 and abs(genId->at(svtTracks[4])) == 321) {
      isBuGood = 2;
    }
    if (isPsi2SGood >= 3 and isBuGood == 2) {
      for (size_t iGoodBu = 0; iGoodBu < goodBus.size(); iGoodBu++) {
        if (crosscomp(svtTracks.data(), goodBus[iGoodBu].data()) and 
          crosscomp(svtTracks.data() + 2, goodBus[iGoodBu].data() + 2) and 
          svtTracks[4] == goodBus[iGoodBu][4]) 
        {
          if (goodBus[iGoodBu][5] == 0) {
            isBuGood = 4;
          }
          else {
            isBuGood = 3;
          }
          goodBuMatch[iGoodBu] += 1;
          nMatchedBuRecoToSim += 1;
        }
      }
    }
    tWriter->svtIsGoodMatch->at(iSvt) = isBuGood;
  }
  
  nGenBu += goodBus.size();
  
  for (size_t iGoodBu = 0; iGoodBu < goodBus.size(); iGoodBu++) {
    if (goodBuMatch[iGoodBu] == 0) {
//      clog << "Unmatched Bu: ";
//      clog << boolalpha;
//      clog << "mu1: " << genTrackMatch[goodBus[iGoodBu][0]] << " "
//           << "mu2: " << genTrackMatch[goodBus[iGoodBu][1]] << " "
//           << "pi1: " << genTrackMatch[goodBus[iGoodBu][2]] << " "
//           << "pi2: " << genTrackMatch[goodBus[iGoodBu][3]] << " "
//           << "k: "   << genTrackMatch[goodBus[iGoodBu][4]] << " ";
      for (int i = 0; i < 5; i++) {
        if (genTrackMatch[goodBus[iGoodBu][i]] == 0) {
//          clog << "Missing track";
          nUnmatchedBuForTrkReasons += 1;
          break;
        }
      }
//      clog << endl;
      continue;
    }
    
    nMatchedBuSimToReco += 1;
    
    if (goodBuMatch[iGoodBu] > 1) {
//      clog << "Bu matched multiple times: ";
//      clog << "mu1: " << genTrackMatch[goodBus[iGoodBu][0]] << " "
//           << "mu2: " << genTrackMatch[goodBus[iGoodBu][1]] << " "
//           << "pi1: " << genTrackMatch[goodBus[iGoodBu][2]] << " "
//           << "pi2: " << genTrackMatch[goodBus[iGoodBu][3]] << " "
//           << "k: "   << genTrackMatch[goodBus[iGoodBu][4]] << " ";
//      clog << endl;
      nDupBu += goodBuMatch[iGoodBu];
    }
  }
  
  tWriter->fill();

  return true;

}


void PDAnalyzer::endJob() {
  tWriter->close();

  clog << "eff: " << nMatchedBuSimToReco*1./nGenBu << " ";
  clog << "trk unmatch: " << nUnmatchedBuForTrkReasons*1./nGenBu << " ";
  clog << "signal fraction: " << nMatchedBuRecoToSim*1./nBu << " ";
  clog << "fake rate: " << 1 - nMatchedBuRecoToSim*1./nBu << " ";
  clog << "dup rate: " << nDupBu*1./nBu << endl;
  
  return;
}

float PDAnalyzer::getCt( int genIndex ) 
{

	const vector <int>& aD = allDaughters(genIndex);
	if( aD.size() == 0 ) 
		return -1;

	unsigned int mthIndex = aD[0];

	// pateracchio 
	if( genId->at( genIndex ) == - genId->at(genMother->at(genIndex)) )
	{
		mthIndex = genMother->at(genIndex );
	}

	ROOT::Math::PxPyPzMVector pGen( 
		(double) genPt->at(genIndex),
		(double) genEta->at(genIndex),
		(double) genPhi->at(genIndex),
		(double) genMass->at(genIndex)
	);

	float dx = genVx->at(genIndex)-genVx->at(mthIndex);
	float dy = genVy->at(genIndex)-genVy->at(mthIndex);
	float dz = genVz->at(genIndex)-genVz->at(mthIndex);

	float ct = sqrt( dx*dx+dy*dy+dz*dz )/pGen.Beta()/pGen.Gamma();
	
	return ct;
 
}

int PDAnalyzer::findPV(ROOT::Math::XYZVector sv, ROOT::Math::XYZVector bP) {
  
  float tmpCos = -999.;
  int bestPV = -1;

  int iPV;
  for ( iPV = 0; iPV < nPVertices; ++iPV ) {
    ROOT::Math::XYZVector tmpPV(pvtX->at(iPV),pvtY->at(iPV),pvtZ->at(iPV));
    auto diff= sv-tmpPV; //vector pointing from PV to SV
    // if (abs(sv.Z()-tmpPV.Z())>0.5) continue; //if the PV has a z distance more than 0.5 cm just skip the vertex
    //diff.SetZ(0.);
    //BsP3.SetZ(0.);
    float cosP=cos(diff.Unit().Dot(bP.Unit()));
    if (cosP > tmpCos){
      tmpCos=cosP;
      bestPV=iPV;
    }
  }  

  return bestPV;
}

void PDAnalyzer::printDecayChain(int iGen, const string& pre) {
    cout<<genId->at(iGen)<<endl;

    const vector <int>& vD = allDaughters(iGen);
    uint ndau = vD.size();
    if(ndau == 0) return;

    bool lastLevel = true;
    for(uint id =0; id<ndau; ++id){
        if ( hasDaughter( vD[id] ) ) {
            lastLevel = false;
            break;
        }
    }

    if( lastLevel ){
        cout<<pre<< "+-> ";
        for( uint id=0; id<ndau; ++id ) {
            int d = vD[id];
            cout<<genId->at(d)<<" ";
        }
        cout<<endl;
        return;
    }

    for( uint id=0; id<ndau; ++id ) {
        int d = vD[id];
        cout<<pre<< "+-> ";
        string prepre(pre);
        if ( id == ndau - 1 ) prepre += "    ";
        else prepre += "|   ";
        printDecayChain( d, prepre );
    }
}

void PDAnalyzer::printCompactDecayChain(int iGen, std::ostream& out) {
    out << abs(genId->at(iGen));

    const vector <int>& vD = allDaughters(iGen);
    uint ndau = vD.size();
    if(ndau == 0) return;
    
    out << "( ";

    for( uint id=0; id<ndau; ++id ) {
        int d = vD[id];
        printCompactDecayChain( d, out );
        out << " ";
    }
    out << ")";
}

bool PDAnalyzer::hasDaughter(int iGen)
{
    const vector <int>& vD = allDaughters(iGen);
    return vD.size()>0 ? true : false;
}

int PDAnalyzer::GetClosestGen( float eta, float phi, float pt ) 
{
    double drb = 0.12;
    double dpb = 0.3; 
    int best = -1;
    
    for( uint i = 0; i < genId->size(); ++i ){
       int id = abs(genId->at(i));
       if( not ( (id > 500 and id < 600) or (id > 5000 and id < 6000)) ) continue;
       float dr = deltaR(eta, phi, genEta->at(i), genPhi->at(i));
       float dpt = fabs(genPt->at(i) - pt)/genPt->at(i);

       if( dr > drb ) continue;
       if( dpt > dpb) continue;

       best = (int) i;
       drb = dr;
    } 

    return best;
}

pair<int, float> PDAnalyzer::GetClosestTrackGen( float eta, float phi, float pt, int charge) 
{
    double drb = 0.5;
    //double dpb = 0.1; 
    int best = -1;
    
    for( uint i = 0; i < genId->size(); ++i ){
       int id = abs(genId->at(i));
       if (genStatus->at(i) != 1) continue;
       if( id != 211 and id != 321 and id != 11 and id != 13 and id != 2212 ) continue;
       if (genCharge->at(i) != charge) continue;
       float dr = deltaR(eta, phi, genEta->at(i), genPhi->at(i));
       //float dpt = fabs(genPt->at(i) - pt)/genPt->at(i);

       if( dr > drb ) continue;
       //if( dpt > dpb) continue;

       best = (int) i;
       drb = dr;
    } 

    return {best, drb};
}

pair<int, float> PDAnalyzer::GetClosestTrackGenFromColl( float eta, float phi, float pt, int charge, const vector<int>& coll) 
{
    double drb = 0.5;
    //double dpb = 0.1; 
    int best = -1;
    
    for( auto i: coll ){
       //if (genCharge->at(i) != charge) continue;
       float dr = deltaR(eta, phi, genEta->at(i), genPhi->at(i));
       //float dpt = fabs(genPt->at(i) - pt)/genPt->at(i);

       if( dr > drb ) continue;
       //if( dpt > dpb) continue;

       best = (int) i;
       drb = dr;
    } 

    return {best, drb};
}

std::optional<std::array<int, 6>> PDAnalyzer::getBuVertexComponents(int iBuGen) {
  std::array<int, 6> ret = {{-1, -1, -1, -1, -1, 0}};
  
  if (abs(genId->at(iBuGen)) != 521) {
    return nullopt;
  }
  
  const vector<int>& vBuDau = allDaughters(iBuGen);
  if (vBuDau.size() < 2) {
    return nullopt;
  }
  int iPsi2SGen = vBuDau[0];
  if (abs(genId->at(iPsi2SGen)) != 100443) {
    return nullopt;
  }
  int iKGen = vBuDau[1];
  if (abs(genId->at(iKGen)) != 321) {
    return nullopt;
  }
  ret[4] = iKGen;
  for (size_t iBuDau = 2; iBuDau < vBuDau.size(); iBuDau++) {
    auto buDau = vBuDau[iBuDau];
    if (abs(genId->at(buDau)) != 22) {
      return nullopt;
    }
    ret[5] += 1;
  }
  
  const vector<int>& vPsi2SDau = allDaughters(iPsi2SGen);
  if (vPsi2SDau.size() < 3) {
    return nullopt;
  }
  int iJPsiGen = vPsi2SDau[0];
  if (abs(genId->at(iJPsiGen)) != 443) {
    return nullopt;
  }
  int iPi1Gen = vPsi2SDau[1];
  if (abs(genId->at(iPi1Gen)) != 211) {
    return nullopt;
  }
  ret[2] = iPi1Gen;
  int iPi2Gen = vPsi2SDau[2];
  if (abs(genId->at(iPi2Gen)) != 211) {
    return nullopt;
  }
  ret[3] = iPi2Gen;
  for (size_t iPsi2SDau = 3; iPsi2SDau < vPsi2SDau.size(); iPsi2SDau++) {
    auto psi2SDau = vPsi2SDau[iPsi2SDau];
    if (abs(genId->at(psi2SDau)) != 22) {
      return nullopt;
    }
    ret[5] += 1;
  }
  
  const vector<int>& vJPsiSDau = allDaughters(iJPsiGen);
  if (vJPsiSDau.size() < 2) {
    return nullopt;
  }
  int iMu1Gen = vJPsiSDau[0];
  if (abs(genId->at(iMu1Gen)) != 13) {
    return nullopt;
  }
  ret[0] = iMu1Gen;
  int iMu2Gen = vJPsiSDau[1];
  if (abs(genId->at(iMu2Gen)) != 13) {
    return nullopt;
  }
  ret[1] = iMu2Gen;
  for (size_t iJPsiDau = 2; iJPsiDau < vJPsiSDau.size(); iJPsiDau++) {
    auto jPsiDau = vJPsiSDau[iJPsiDau];
    if (abs(genId->at(jPsiDau)) != 22) {
      return nullopt;
    }
    ret[5] += 1;
  }
  
  return ret;
}

std::optional<std::array<int, 5>> PDAnalyzer::getPsi2SVertexComponents(int iPsi2SGen) {
  std::array<int, 5> ret = {{-1, -1, -1, -1, 0}};
  
  if (abs(genId->at(iPsi2SGen)) != 100443) {
    return nullopt;
  }
  
  const vector<int>& vPsi2SDau = allDaughters(iPsi2SGen);
  if (vPsi2SDau.size() < 3) {
    return nullopt;
  }
  int iJPsiGen = vPsi2SDau[0];
  if (abs(genId->at(iJPsiGen)) != 443) {
    return nullopt;
  }
  int iPi1Gen = vPsi2SDau[1];
  if (abs(genId->at(iPi1Gen)) != 211) {
    return nullopt;
  }
  ret[2] = iPi1Gen;
  int iPi2Gen = vPsi2SDau[2];
  if (abs(genId->at(iPi2Gen)) != 211) {
    return nullopt;
  }
  ret[3] = iPi2Gen;
  for (size_t iPsi2SDau = 3; iPsi2SDau < vPsi2SDau.size(); iPsi2SDau++) {
    auto psi2SDau = vPsi2SDau[iPsi2SDau];
    if (abs(genId->at(psi2SDau)) != 22) {
      return nullopt;
    }
    ret[4] += 1;
  }
  
  const vector<int>& vJPsiSDau = allDaughters(iJPsiGen);
  if (vJPsiSDau.size() < 2) {
    return nullopt;
  }
  int iMu1Gen = vJPsiSDau[0];
  if (abs(genId->at(iMu1Gen)) != 13) {
    return nullopt;
  }
  ret[0] = iMu1Gen;
  int iMu2Gen = vJPsiSDau[1];
  if (abs(genId->at(iMu2Gen)) != 13) {
    return nullopt;
  }
  ret[1] = iMu2Gen;
  for (size_t iJPsiDau = 2; iJPsiDau < vJPsiSDau.size(); iJPsiDau++) {
    auto jPsiDau = vJPsiSDau[iJPsiDau];
    if (abs(genId->at(jPsiDau)) != 22) {
      return nullopt;
    }
    ret[4] += 1;
  }
  
  return ret;
}

std::optional<std::array<int, 3>> PDAnalyzer::getJPsiVertexComponents(int iJPsiGen) {
  std::array<int, 3> ret = {{-1, -1, 0}};
  
  if (abs(genId->at(iJPsiGen)) != 443) {
    return nullopt;
  }
  
  const vector<int>& vJPsiSDau = allDaughters(iJPsiGen);
  if (vJPsiSDau.size() < 2) {
    return nullopt;
  }
  int iMu1Gen = vJPsiSDau[0];
  if (abs(genId->at(iMu1Gen)) != 13) {
    return nullopt;
  }
  ret[0] = iMu1Gen;
  int iMu2Gen = vJPsiSDau[1];
  if (abs(genId->at(iMu2Gen)) != 13) {
    return nullopt;
  }
  ret[1] = iMu2Gen;
  for (size_t iJPsiDau = 2; iJPsiDau < vJPsiSDau.size(); iJPsiDau++) {
    auto jPsiDau = vJPsiSDau[iJPsiDau];
    if (abs(genId->at(jPsiDau)) != 22) {
      return nullopt;
    }
    ret[2] += 1;
  }
  
  return ret;
}



