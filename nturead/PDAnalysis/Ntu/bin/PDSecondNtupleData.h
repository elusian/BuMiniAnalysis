#ifndef PDSecondNtupleData_h
#define PDSecondNtupleData_h
#include <vector>
#include "NtuTool/Common/interface/TreeWrapper.h"
using namespace std;

class PDSecondNtupleData: public virtual TreeWrapper {

public:

void Reset() { autoReset(); }

PDSecondNtupleData() {

    trkPt           = new vector <float>;
    trkEta          = new vector <float>;
    trkPhi          = new vector <float>;
    trkCharge       = new vector <int>;
    trkDxy          = new vector <float>;
    trkDz           = new vector <float>;
    trkExy          = new vector <float>;
    trkEz           = new vector <float>;
    trkPfcPt        = new vector <float>;
    trkPfcEta       = new vector <float>;
    trkPfcPhi       = new vector <float>;
    trkGenPt        = new vector <float>;
    trkGenEta       = new vector <float>;
    trkGenPhi       = new vector <float>;
    trkNVHMuon      = new vector <int>;
    trkNVHPixel     = new vector <int>;
    trkNVHTracker   = new vector <int>;
    trkNVHAll       = new vector <int>;
    trkLPPixel      = new vector <int>;
    trkLPStrips     = new vector <int>;
    trkLPTracker    = new vector <int>;
    trkMIH          = new vector <int>;
    trkKind         = new vector <int>;
    trkQOverPError  = new vector <float>;
    trkLambdaError  = new vector <float>;
    trkPhiError     = new vector <float>;
    trkDszError     = new vector <float>;
    trkEtaError     = new vector <float>;
    trkPtError      = new vector <float>;
    trkDxyDszCov    = new vector <float>;
    trkLambdaDszCov = new vector <float>;
    trkPhiDxyCov    = new vector <float>;
    trkGenId        = new vector <int>;
    trkGenMatchRadius = new vector <float>;
    trkIsMu         = new vector <int>;
    
    
    svtPt = new vector<float>;
    svtEta = new vector<float>;
    svtPhi = new vector<float>;
    svtMass = new vector<float>;
    svtChi2 = new vector<float>;
    svtNdof = new vector<int>;
    svtDist2D = new vector<float>;
    svtSigma2D = new vector<float>;
    svtDist3D = new vector<float>;
    svtSigma3D = new vector<float>;
    svtType = new vector<int>;
    svtCt = new vector<float>;
    svtCtErr = new vector<float>;
    svtDecayCosTheta = new vector<float>;
    svtDecayCosPsi = new vector<float>;
    svtDecayPhi = new vector<float>;
    svtPsi2SIndex = new vector<int>;
    svtJPsiIndex = new vector<int>;
    svtPi1Index = new vector<int>;
    svtPi2Index = new vector<int>;
    svtMu1Index = new vector<int>;
    svtMu2Index = new vector<int>;
    svtKIndex = new vector<int>;
    svtGenId = new vector<int>;
    svtDecayId = new vector<int>;
    svtIsGoodMatch = new vector<int>;
    svtGenMatchRadius = new vector<float>;
    svtGenMatchDpt = new vector<float>;
    svtGenMatchIndex = new vector<int>;
    
    
    genPt = new vector<float>;
    genEta = new vector<float>;
    genPhi = new vector<float>;
    genMass = new vector<float>;
    genCt = new vector<float>;
    genId = new vector<int>;
    genCharge = new vector<int>;

}
virtual ~PDSecondNtupleData() {
}

void initTree() {
    treeName = "PDsecondTree";

    setBranch( "svtPt", &svtPt , 8192, 99, &b_svtPt );
    setBranch( "svtEta", &svtEta , 8192, 99, &b_svtEta );
    setBranch( "svtPhi", &svtPhi , 8192, 99, &b_svtPhi );
    setBranch( "svtMass", &svtMass , 8192, 99, &b_svtMass );
    setBranch( "svtChi2", &svtChi2 , 8192, 99, &b_svtChi2 );
    setBranch( "svtNdof", &svtNdof , 8192, 99, &b_svtNdof );
    setBranch( "svtDist2D", &svtDist2D , 8192, 99, &b_svtDist2D );
    setBranch( "svtSigma2D", &svtSigma2D , 8192, 99, &b_svtSigma2D );
    setBranch( "svtDist3D", &svtDist3D , 8192, 99, &b_svtDist3D );
    setBranch( "svtSigma3D", &svtSigma3D , 8192, 99, &b_svtSigma3D );
    setBranch( "svtCt", &svtCt , 8192, 99, &b_svtCt );
    setBranch( "svtCtErr", &svtCtErr , 8192, 99, &b_svtCtErr );
    setBranch( "svtType", &svtType , 8192, 99, &b_svtType );
    setBranch( "svtPsi2SIndex", &svtPsi2SIndex, 8192, 99, &b_svtPsi2SIndex);
    setBranch( "svtJPsiIndex", &svtJPsiIndex, 8192, 99, &b_svtJPsiIndex);
    setBranch( "svtPi1Index", &svtPi1Index, 8192, 99, &b_svtPi1Index);
    setBranch( "svtPi2Index", &svtPi2Index, 8192, 99, &b_svtPi2Index);
    setBranch( "svtMu1Index", &svtMu1Index, 8192, 99, &b_svtMu1Index);
    setBranch( "svtMu2Index", &svtMu2Index, 8192, 99, &b_svtMu2Index);
    setBranch( "svtKIndex", &svtKIndex, 8192, 99, &b_svtKIndex);
    setBranch( "svtGenId", &svtGenId, 8192, 99, &b_svtGenId);
    setBranch( "svtDecayId", &svtDecayId, 8192, 99, &b_svtDecayId);
    setBranch( "svtIsGoodMatch", &svtIsGoodMatch, 8192, 99, &b_svtIsGoodMatch);
    setBranch( "svtGenMatchRadius", &svtGenMatchRadius, 8192, 99, &b_svtGenMatchRadius);
    setBranch( "svtGenMatchDpt", &svtGenMatchDpt, 8192, 99, &b_svtGenMatchDpt);
    setBranch( "svtGenMatchIndex", &svtGenMatchIndex, 8192, 99, &b_svtGenMatchIndex);
    

    setBranch( "evtNumber", &evtNumber, "evtNumber/I", &b_evtNumber );

    setBranch( "trkPt", &trkPt , 8192, 99, &b_trkPt );
    setBranch( "trkEta", &trkEta , 8192, 99, &b_trkEta );
    setBranch( "trkPhi", &trkPhi , 8192, 99, &b_trkPhi );
    setBranch( "trkCharge", &trkCharge , 8192, 99, &b_trkCharge );
    setBranch( "trkDxy", &trkDxy , 8192, 99, &b_trkDxy );
    setBranch( "trkDz", &trkDz , 8192, 99, &b_trkDz );
    setBranch( "trkExy", &trkExy , 8192, 99, &b_trkExy );
    setBranch( "trkEz", &trkEz , 8192, 99, &b_trkEz );
    setBranch( "trkPfcPt", &trkPfcPt , 8192, 99, &b_trkPfcPt );
    setBranch( "trkPfcEta", &trkPfcEta , 8192, 99, &b_trkPfcEta );
    setBranch( "trkPfcPhi", &trkPfcPhi , 8192, 99, &b_trkPfcPhi );
    setBranch( "trkGenPt", &trkGenPt , 8192, 99, &b_trkGenPt );
    setBranch( "trkGenEta", &trkGenEta , 8192, 99, &b_trkGenEta );
    setBranch( "trkGenPhi", &trkGenPhi , 8192, 99, &b_trkGenPhi );
    setBranch( "trkNVHMuon", &trkNVHMuon, 8192, 99, &b_trkNVHMuon);
    setBranch( "trkNVHPixel", &trkNVHPixel, 8192, 99, &b_trkNVHPixel);
    setBranch( "trkNVHTracker", &trkNVHTracker, 8192, 99, &b_trkNVHTracker);
    setBranch( "trkNVHAll", &trkNVHAll, 8192, 99, &b_trkNVHAll);
    setBranch( "trkLPPixel", &trkLPPixel, 8192, 99, &b_trkLPPixel);
    setBranch( "trkLPStrips", &trkLPStrips, 8192, 99, &b_trkLPStrips);
    setBranch( "trkLPTracker", &trkLPTracker, 8192, 99, &b_trkLPTracker);
    setBranch( "trkMIH", &trkMIH, 8192, 99, &b_trkMIH);
    setBranch( "trkKind", &trkKind, 8192, 99, &b_trkKind);
    setBranch( "trkQOverPError", &trkQOverPError, 8192, 99, &b_trkQOverPError);
    setBranch( "trkLambdaError", &trkLambdaError, 8192, 99, &b_trkLambdaError);
    setBranch( "trkPhiError", &trkPhiError, 8192, 99, &b_trkPhiError);
    setBranch( "trkDszError", &trkDszError, 8192, 99, &b_trkDszError);
    setBranch( "trkEtaError", &trkEtaError, 8192, 99, &b_trkEtaError);
    setBranch( "trkPtError", &trkPtError, 8192, 99, &b_trkPtError);
    setBranch( "trkDxyDszCov", &trkDxyDszCov, 8192, 99, &b_trkDxyDszCov);
    setBranch( "trkLambdaDszCov", &trkLambdaDszCov, 8192, 99, &b_trkLambdaDszCov);
    setBranch( "trkPhiDxyCov", &trkPhiDxyCov, 8192, 99, &b_trkPhiDxyCov);
    setBranch( "trkGenId", &trkGenId, 8192, 99, &b_trkGenId);
    setBranch( "trkGenMatchRadius", &trkGenMatchRadius, 8192, 99, &b_trkGenMatchRadius);
    setBranch( "trkIsMu", &trkIsMu, 8192, 99, &b_trkIsMu);
    
    setBranch( "genPt", &genPt , 8192, 99, &b_genPt );
    setBranch( "genEta", &genEta , 8192, 99, &b_genEta );
    setBranch( "genPhi", &genPhi , 8192, 99, &b_genPhi );
    setBranch( "genMass", &genMass , 8192, 99, &b_genMass );
    setBranch( "genCt", &genCt , 8192, 99, &b_genCt );
    setBranch( "genId", &genId , 8192, 99, &b_genId );
    setBranch( "genCharge", &genCharge , 8192, 99, &b_genCharge );

}

int evtNumber;
TBranch *b_evtNumber;

vector <float> *trkPt, *trkEta, *trkPhi, *trkPfcPt, *trkPfcEta, *trkPfcPhi, *trkGenPt, *trkGenEta, *trkGenPhi;
vector <float> *trkDxy, *trkDz, *trkExy, *trkEz;
vector <int> *trkCharge;
vector <int> *trkNVHMuon, *trkNVHPixel, *trkNVHTracker, *trkNVHAll, *trkLPPixel, *trkLPStrips, *trkLPTracker, *trkMIH, *trkKind;
vector <float> *trkQOverPError, *trkLambdaError, *trkPhiError, *trkDszError, *trkEtaError, *trkPtError;
vector <float> *trkDxyDszCov, *trkLambdaDszCov, *trkPhiDxyCov;
vector <int> *trkGenId, *trkIsMu;
vector <float> *trkGenMatchRadius;
TBranch *b_trkPt, *b_trkEta, *b_trkPhi, *b_trkPfcPt, *b_trkPfcEta, *b_trkPfcPhi, *b_trkGenPt, *b_trkGenEta, *b_trkGenPhi;
TBranch *b_trkDxy, *b_trkDz, *b_trkExy, *b_trkEz;
TBranch *b_trkCharge;
TBranch *b_trkNVHMuon, *b_trkNVHPixel, *b_trkNVHTracker, *b_trkNVHAll, *b_trkLPPixel, *b_trkLPStrips, *b_trkLPTracker, *b_trkMIH, *b_trkKind;
TBranch *b_trkQOverPError, *b_trkLambdaError, *b_trkPhiError, *b_trkDszError, *b_trkEtaError, *b_trkPtError;
TBranch *b_trkDxyDszCov, *b_trkLambdaDszCov, *b_trkPhiDxyCov;
TBranch *b_trkGenId, *b_trkGenMatchRadius, *b_trkIsMu;


vector <float> *svtPt, *svtEta, *svtPhi, *svtMass, *svtGenMatchRadius, *svtGenMatchDpt;
vector <float> *svtChi2, *svtDist2D, *svtSigma2D, *svtDist3D, *svtSigma3D, *svtCt, *svtCtErr;
vector <float> *svtDecayCosTheta, *svtDecayCosPsi, *svtDecayPhi;
vector <int> *svtPsi2SIndex, *svtJPsiIndex;
vector <int> *svtPi1Index, *svtPi2Index, *svtMu1Index, *svtMu2Index, *svtKIndex;
vector <int> *svtNdof, *svtType;
vector <int> *svtGenId, *svtDecayId, *svtIsGoodMatch, *svtGenMatchIndex;
TBranch *b_svtPt, *b_svtEta, *b_svtPhi, *b_svtMass, *b_svtGenMatchRadius, *b_svtGenMatchDpt;
TBranch *b_svtChi2, *b_svtDist2D, *b_svtSigma2D, *b_svtDist3D, *b_svtSigma3D, *b_svtCt, *b_svtCtErr;
TBranch *b_svtDecayCosTheta, *b_svtDecayCosPsi, *b_svtDecayPhi;
TBranch *b_svtPsi2SIndex, *b_svtJPsiIndex;
TBranch *b_svtPi1Index, *b_svtPi2Index, *b_svtMu1Index, *b_svtMu2Index, *b_svtKIndex;
TBranch *b_svtNdof, *b_svtType;
TBranch *b_svtGenId, *b_svtDecayId, *b_svtIsGoodMatch, *b_svtGenMatchIndex;

vector <float> *genPt, *genEta, *genPhi, *genMass;
vector <float> *genCt;
vector <int> *genId, *genCharge;
TBranch *b_genPt, *b_genEta, *b_genPhi, *b_genMass;
TBranch *b_genCt;
TBranch *b_genId, *b_genCharge;

private:

PDSecondNtupleData ( const PDSecondNtupleData& a );
PDSecondNtupleData& operator=( const PDSecondNtupleData& a );

};

#endif

