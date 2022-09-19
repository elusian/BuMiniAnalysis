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
  
  convSpheCart(muoPt, muoEta, muoPhi, muoPx, muoPy, muoPz);
  convSpheCart(trkPt, trkEta, trkPhi, trkPx, trkPy, trkPz);
  convSpheCart(pfcPt, pfcEta, pfcPhi, pfcPx, pfcPy, pfcPz);
  
  map <int, int> vertexMap;
  size_t saved = 0;
  for (int iSvt = 0; iSvt < nSVertices; iSvt++) {
  
    auto type = svtType->at(iSvt);
    if (type != PDEnumString::svtBuPsi2SK and type != PDEnumString::svtPsi2S and type != PDEnumString::svtJPsi) continue;
    
    vertexMap[iSvt] = saved;
    saved++;
  }
  
  std::vector<int> jpsiVerts;
  std::unordered_map<int, int> jpsiAssoc;
  for (int iSvt = 0; iSvt < nSVertices; iSvt++) {
    if (type != PDEnumString::svtJPsi) {
        continue;
    }
    jpsiVerts.push_back(iSvt);
    jpsiAssoc[iSvt] = jpsiVerts.size() - 1;
  }
  
  std::vector<int> psi2SVerts(jpsiVerts.size(), -1);
  std::vector<double> psi2SProb(jpsiVerts.size(), 0);
  std::unordered_map<int, int> psi2SAssoc;
  for (int iSvt = 0; iSvt < nSVertices; iSvt++) {
    if (type != PDEnumString::svtPsi2S) {
        continue;
    }
    int iJPsi  = (subVtxFromSV(iSvt)).at(0);
    auto slot = jpsiAssoc[iJPsi];
    
    auto newProb = TMath::Prob(svtChi2->at(iSvt), svtNDOF->at(iSvt));
    auto bestProb = psi2SProb[slot];
    if (newProb > bestProb) {
        psi2SVerts[slot] = iSvt;
        psi2SProb[slot] = newProb;
        psi2SAssoc[iSvt] = slot;
    }
  }
  
  std::vector<int> buVerts(jpsiVerts.size(), -1);
  std::vector<double> buProb(jpsiVerts.size(), 0);
  for (int iSvt = 0; iSvt < nSVertices; iSvt++) {
    if (type != PDEnumString::svtBuPsi2SK) {
        continue;
    }
    int iPsi2S  = (subVtxFromSV(iSvt)).at(0);
    auto slot = psi2SAssoc[iPsi2S];
    
    auto newProb = TMath::Prob(svtChi2->at(iSvt), svtNDOF->at(iSvt));
    auto bestProb = buProb[slot];
    if (newProb > bestProb) {
        buVerts[slot] = iSvt;
        buProb[slot] = newProb;
    }
  }
  
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
    
    if (type == PDEnumString::svtBuPsi2SK) {
//      clog << subVtxFromSV(iSvt).size() << endl;
//      for (unsigned int i = 0; i < subVtxFromSV(iSvt).size(); i++) {
//        int iPsi2S = (subVtxFromSV(iSvt)).at(i);
//        clog << svtType->at(iPsi2S) << " ";
//      }
//      clog << endl;
    
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
        throw 1;
      }
      
      k_idx.erase(remove_if(k_idx.begin(), k_idx.end(), [&pi_idx, &mu_idx](int trkId){return trkId == mu_idx[0] or trkId == mu_idx[1] or trkId == pi_idx[0] or trkId == pi_idx[1];}), k_idx.end());
      if (k_idx.size() != 1) {
        throw 2;
      }
      
      tWriter->svtPi1Index->push_back(pi_idx[0]);
      tWriter->svtPi2Index->push_back(pi_idx[1]);
      tWriter->svtMu1Index->push_back(mu_idx[0]);
      tWriter->svtMu2Index->push_back(mu_idx[1]);
      tWriter->svtKIndex->push_back(k_idx[0]);
      
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
  
  for (int iGen = 0; iGen < nGenP; iGen++) {
    auto id = abs(genId->at(iGen));
    
    if (id != 521 and id != 443 and id != 30443) continue;
    
    tWriter->genPt->push_back(genPt->at(iGen));
    tWriter->genEta->push_back(genEta->at(iGen));
    tWriter->genPhi->push_back(genPhi->at(iGen));
    tWriter->genMass->push_back(genMass->at(iGen));
    tWriter->genCt->push_back( getCt( iGen ));
    tWriter->genCharge->push_back(genCharge->at(iGen));
    tWriter->genId->push_back(genId->at(iGen));
  }
  
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
  }
  
  tWriter->fill();

  return true;

}


void PDAnalyzer::endJob() {
  tWriter->close();
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


