
#include <TH1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TStyle.h>
#include <TMath.h>
#include <TVector3.h>

#include <RooMsgService.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooJohnson.h>
#include <RooExponential.h>
#include <RooGaussian.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooCategory.h>
#include <RooFitResult.h>

#include <map>

#include "pulls.C"

using namespace std;

map<string, int> colors = [](){
    map<string, int> m;
    m["no mkFit"] = kBlue;
    m["mkFit"] = kBlack;
    
    return m;
}();

map<string, string> xaxis, yaxis;

void drawAll(TString imgName, map<string, TH1*> hists, TString baseLabel, bool log = false) {
    TCanvas c;
    TLegend l(0.72, 0.76, 0.92, 0.91);
    bool setHist = false;
    TH1* mainHist = nullptr;
    double maxVal = 0;
    for (auto&& [name, hist]: hists) {
        hist->SetLineColor(colors[name]);
        hist->SetLineWidth(2);
        hist->SetYTitle(yaxis[imgName.Data()].c_str());
        hist->SetXTitle(xaxis[imgName.Data()].c_str());
//        if (not mainHist) {
//            hist->Draw("HIST");
//            mainHist = hist;
//        }
//        else {
//            hist->Draw("HIST SAME");
//        }
//        l.AddEntry(hist, name.c_str());
//        auto newMax = hist->GetMaximum();
//        if (newMax > maxVal) {
//            maxVal = newMax;
//        }
    }
    auto ratio = new TRatioPlot(hists["mkFit"], hists["no mkFit"]);
    ratio->Draw();
    l.AddEntry(hists["mkFit"], "mkFit");
    l.AddEntry(hists["no mkFit"], "no mkFit");
    
    ratio->GetLowerRefYaxis()->SetRangeUser(0.2, 1.8);
    
//    if (not log) mainHist->GetYaxis()->SetRangeUser(0, maxVal*1.2);
//    else mainHist->GetYaxis()->SetRangeUser(maxVal*0.001, maxVal*1.2);
    //gPad->BuildLegend(0.72, 0.76, 0.92, 0.91);
    l.Draw();
    if (log) ratio->GetUpperPad()->SetLogy();
    c.SaveAs("plotSvt/" + baseLabel + "/" + imgName + ".png");
}

void sbSub(TH1* h, TH1* sb, double sbWeight, double sigRescale) {
    h->Add(sb, -sbWeight);
    h->Scale(sigRescale);
}

void drawSvtVars(TString baselabel) {
    gStyle->SetOptTitle(false);
    gStyle->SetOptStat(false);
    gStyle->SetHistMinimumZero();
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    gErrorIgnoreLevel = kWarning;
    
    TString labelRef = baselabel + "_base";
    TString labelMkFit = baselabel + "_mkFit";
    
    map<string, const char*> names;
    names[labelRef.Data()] = "no mkFit";
    names[labelMkFit.Data()] = "mkFit";
    
    vector<pair<string, const char*>> orderedNames;
    
    for (string name: {labelRef.Data(), labelMkFit.Data()}) {
        orderedNames.push_back(make_pair(name, names[name]));
    }
    
    map<string, TFile*> files;
    map<string, TTree*> trees;
    
    for (auto&& [name, label]: orderedNames) {
        files[name] = TFile::Open(Form("%s/second.root", name.c_str()));
        
        trees[name] = (TTree*)files[name]->Get("PDsecondTree");
    }
    
    map<string, pair<double, double>> muVal;
    map<string, pair<double, double>> s1Val;
    map<string, pair<double, double>> s2Val;
    map<string, pair<double, double>> lambdaVal;
    map<string, pair<double, double>> gammaVal;
    map<string, pair<double, double>> deltaVal;
    map<string, pair<double, double>> yieldVal;
    map<string, pair<double, double>> bkgyieldVal;
    map<string, pair<double, double>> fracVal;
    map<string, pair<double, double>> coreFracs;
    
    map<string, map<string, RooRealVar*>> vars;
    map<string, RooFitResult*> res;
    
    double massLimitLow = 5.16;
    double massSbLimitLow = 5.18;
    double massCentLimitLow = 5.25;
    double massCentLimitHigh = 5.31;
    double massSbLimitHigh = 5.35;
    double massLimitHigh = 5.41;
    
    float kRadiusCut = 0.04;
    float kPtDiff = 0.3;
//    float piRadiusCut = 0.06;
//    float piPtDiff = 0.15;
    float piRadiusCut = 0.3;
    float piPtDiff = 0.3;
    float muRadiusCut = 0.04;
    float muPtDiff = 0.15;
    
    map<string, TH1*> hsvtMass;
    
    map<string, TH1*> hsvtPt;
    map<string, TH1*> hsvtEta;
    
    map<string, TH1*> hsvtCt;
    map<string, TH1*> hsvtCtErr;
    
    map<string, TH1*> hsvtMassJPsi;
    map<string, TH1*> hsvtPtJPsi;
    map<string, TH1*> hsvtEtaJPsi;
    map<string, TH1*> hsvtMassPsi2S;
    map<string, TH1*> hsvtPtPsi2S;
    map<string, TH1*> hsvtEtaPsi2S;
    
    map<string, TH1*> hsvtDist3D;
    map<string, TH1*> hsvtSign3D;
    map<string, TH1*> hsvtSigma3D;
    
    map<string, TH1*> hsvtProb;
    
#define trkHist(varname, part)\
    map<string, TH1*> h##varname##part##All;\
    map<string, TH1*> h##varname##part##Min;\
    map<string, TH1*> h##varname##part##Max;
    
    trkHist(trkPt, Pi)
    trkHist(trkEta, Pi)
    trkHist(trkDxy, Pi)
    trkHist(trkDz, Pi)
    trkHist(trkExy, Pi)
    trkHist(trkEz, Pi)
    trkHist(trkNVHAll, Pi)
    trkHist(trkNVHPixel, Pi)
    
    map<string, TH1*> htrkPtK;
    map<string, TH1*> htrkEtaK;
    map<string, TH1*> htrkDxyK;
    map<string, TH1*> htrkDzK;
    map<string, TH1*> htrkExyK;
    map<string, TH1*> htrkEzK;
    map<string, TH1*> htrkNVHAllK;
    map<string, TH1*> htrkNVHPixelK;
    
    TH1* hgenPt = new TH1D("hgenPt", "", 50, 0, 50);
    TH1* hgenEta = new TH1D("hgenEta", "", 50, -8, 8);
    
    map<string, TH1*> hPiAngle;
        
    RooRealVar mass("mass", "M(#pi^{+}#pi^{-}#mu^{+}#mu^{-}K^{#pm})", massLimitLow, massLimitHigh, "GeV");
    mass.setRange("sbLeft", massLimitLow, massSbLimitLow);
    mass.setRange("sbRight", massSbLimitHigh, massLimitHigh);
    mass.setRange("center", massSbLimitLow, massSbLimitHigh);
    mass.setRange("core", massCentLimitLow, massCentLimitHigh);
    
    RooCategory genMatchState("genMatchState", "");
    genMatchState.defineType("matched", 1);
    genMatchState.defineType("unmatched", 0);
    
    auto comparisonFrame = mass.frame();
    auto comparisonFrameNorm = mass.frame();
    
    map<string, float> discardedBkgYield;
    
    unordered_multimap<int, size_t> passedEvents;
    size_t nEvents = trees[labelRef.Data()]->GetEntries();
    //size_t eventRadius = 100;
    
    for (auto&& [name, label]: orderedNames) {
        clog << "Processing " << label << endl;
#define initSvtHist(varname, bins, low, high, xaxisLabel, udm)\
        h##varname[label] = new TH1D(Form("h" #varname "%s", name.c_str()), label, bins, low, high);\
        h##varname[label]->Sumw2();\
        auto h##varname##SB = new TH1D(Form("h" #varname "%sSB", name.c_str()), label, bins, low, high);\
        h##varname##SB->Sumw2();\
        xaxis[#varname] = xaxisLabel;\
        yaxis[#varname] = Form("Events / (%g %s)", (high-low)/bins, udm);
        
        initSvtHist(svtMass, 50, massLimitLow, massLimitHigh, "Bu mass (GeV)", "GeV");
        
        initSvtHist(svtPt, 40, 0., 50, "Bu p_{T} (GeV)", "GeV");
        initSvtHist(svtEta, 25, -3., 3, "Bu #eta", "");
        
        initSvtHist(svtCt, 25, -0.01, 0.4, "Bu ct (cm)", "cm");
        initSvtHist(svtCtErr, 25, 0., 0.005, "Bu #sigma(ct) (cm)", "cm");
        
        initSvtHist(svtMassPsi2S, 50, 3.6, 3.75, "#psi' mass (GeV)", "GeV");
        initSvtHist(svtPtPsi2S, 40, 0, 50., "#psi' p_{T} (GeV)", "GeV");
        initSvtHist(svtEtaPsi2S, 25, -3, 3., "#psi' #eta", "");
        initSvtHist(svtMassJPsi, 50, 2.9, 3.3, "J/#psi mass (GeV)", "GeV");
        initSvtHist(svtPtJPsi, 40, 0, 50., "J/#psi p_{T} (GeV)", "GeV");
        initSvtHist(svtEtaJPsi, 25, -3, 3., "J/#psi #eta", "");
        
        initSvtHist(svtDist3D, 25, 0, 1., "Bu PV-SV distance (cm)", "cm");
        initSvtHist(svtSigma3D, 25, 0, 0.2, "Bu PV-SV sigma (cm)", "cm");
        initSvtHist(svtSign3D, 100, 0, 50., "Bu PV-SV significance", "");
        
        initSvtHist(svtProb, 50, 0, 1., "Bu SV probability", "");
        

#define initTrkHist(varname, part, bins, low, high, xaxisLabel, udm)\
        h##varname##part##All[label] = new TH1D(Form("h" #varname #part "All%s", name.c_str()), label, bins, low, high);\
        h##varname##part##All[label]->Sumw2();\
        auto h##varname##part##AllSB = new TH1D(Form("h" #varname #part "All%sSB", name.c_str()), label, bins, low, high);\
        h##varname##part##AllSB->Sumw2();\
        if (udm != TString("")) xaxis[#varname #part "All"] = Form("%s, all (%s)", xaxisLabel, udm);\
        else xaxis[#varname #part "All"] = Form("%s, all", xaxisLabel);\
        yaxis[#varname #part "All"] = Form("Events / (%g %s)", (high-low)/bins, udm);\
        h##varname##part##Min[label] = new TH1D(Form("h" #varname #part "Min%s", name.c_str()), label, bins, low, high);\
        h##varname##part##Min[label]->Sumw2();\
        auto h##varname##part##MinSB = new TH1D(Form("h" #varname #part "Min%sSB", name.c_str()), label, bins, low, high);\
        h##varname##part##MinSB->Sumw2();\
        if (udm != TString("")) xaxis[#varname #part "Min"] = Form("%s, lowest p_{T} (%s)", xaxisLabel, udm);\
        else xaxis[#varname #part "Min"] = Form("%s, lowest p_{T}", xaxisLabel);\
        yaxis[#varname #part "Min"] = Form("Events / (%g %s)", (high-low)/bins, udm);\
        h##varname##part##Max[label] = new TH1D(Form("h" #varname #part "Max%s", name.c_str()), label, bins, low, high);\
        h##varname##part##Max[label]->Sumw2();\
        auto h##varname##part##MaxSB = new TH1D(Form("h" #varname #part "Max%sSB", name.c_str()), label, bins, low, high);\
        h##varname##part##MaxSB->Sumw2();\
        if (udm != TString("")) xaxis[#varname #part "Max"] = Form("%s, highest p_{T} (%s)", xaxisLabel, udm);\
        else xaxis[#varname #part "Max"] = Form("%s, highest p_{T}", xaxisLabel);\
        yaxis[#varname #part "Max"] = Form("Events / (%g %s)", (high-low)/bins, udm);
        
        initTrkHist(trkPt, Pi, 25, 0, 5., "Pion p_{T}", "GeV")
        initTrkHist(trkEta, Pi, 25, -3, 3., "Pion #eta", "")
        initTrkHist(trkDxy, Pi, 25, -0.1, 0.1, "Pion d_{xy}", "cm")
        initTrkHist(trkExy, Pi, 25, 0, 0.2, "Pion #sigma(d_{xy})", "cm")
        initTrkHist(trkDz, Pi, 25, -0.1, 0.1, "Pion d_{z}", "cm")
        initTrkHist(trkEz, Pi, 25, 0, 0.2, "Pion #sigma(d_{z})", "cm")
        initTrkHist(trkNVHAll, Pi, 51, -0.5, 50.5, "Pion #hits", "")
        initTrkHist(trkNVHPixel, Pi, 31, -0.5, 30.5, "Pion #px hits", "")
        
        initSvtHist(trkPtK, 25, 0, 10., "Kaon p_{T} (GeV)", "GeV")
        initSvtHist(trkEtaK, 25, -3, 3., "Kaon #eta", "")
        initSvtHist(trkDxyK, 25, -0.1, 0.1, "Kaon d_{xy} (cm)", "cm")
        initSvtHist(trkExyK, 25, 0, 0.04, "Kaon #sigma(d_{xy}) (cm)", "cm")
        initSvtHist(trkDzK, 25, -0.1, 0.1, "Kaon d_{z} (cm)", "cm")
        initSvtHist(trkEzK, 25, 0, 0.05, "Kaon #sigma(d_{z}) (cm)", "cm")
        initSvtHist(trkNVHAllK, 51, -0.5, 50.5, "Kaon #hits", "")
        initSvtHist(trkNVHPixelK, 31, -0.5, 30.5, "Kaon #px hits", "")
        
        initSvtHist(PiAngle, 20, 0, 0.6, "#pi#pi angle (rad)", "rad")
    
        TTreeReader fReader(trees[name]);
    
        TTreeReaderArray<float> svtPt = {fReader, "svtPt"};
        TTreeReaderArray<float> svtEta = {fReader, "svtEta"};
        TTreeReaderArray<float> svtPhi = {fReader, "svtPhi"};
        TTreeReaderArray<float> svtMass = {fReader, "svtMass"};
        TTreeReaderArray<float> svtChi2 = {fReader, "svtChi2"};
        TTreeReaderArray<int> svtNdof = {fReader, "svtNdof"};
        TTreeReaderArray<float> svtDist2D = {fReader, "svtDist2D"};
        TTreeReaderArray<float> svtSigma2D = {fReader, "svtSigma2D"};
        TTreeReaderArray<float> svtDist3D = {fReader, "svtDist3D"};
        TTreeReaderArray<float> svtSigma3D = {fReader, "svtSigma3D"};
        TTreeReaderArray<float> svtCt = {fReader, "svtCt"};
        TTreeReaderArray<float> svtCtErr = {fReader, "svtCtErr"};
        TTreeReaderArray<int> svtType = {fReader, "svtType"};
        TTreeReaderArray<int> svtPsi2SIndex = {fReader, "svtPsi2SIndex"};
        TTreeReaderArray<int> svtJPsiIndex = {fReader, "svtJPsiIndex"};
        TTreeReaderArray<int> svtIsGoodMatch = {fReader, "svtIsGoodMatch"};
        TTreeReaderArray<float> svtGenMatchRadius = {fReader, "svtGenMatchRadius"};
        TTreeReaderArray<float> svtGenMatchDpt = {fReader, "svtGenMatchDpt"};
        TTreeReaderArray<int> svtPi1Index = {fReader, "svtPi1Index"};
        TTreeReaderArray<int> svtPi2Index = {fReader, "svtPi2Index"};
        TTreeReaderArray<int> svtMu1Index = {fReader, "svtMu1Index"};
        TTreeReaderArray<int> svtMu2Index = {fReader, "svtMu2Index"};
        TTreeReaderArray<int> svtKIndex = {fReader, "svtKIndex"};
        TTreeReaderArray<int> svtGenId = {fReader, "svtGenId"};
        TTreeReaderValue<Int_t> evtNumber = {fReader, "evtNumber"};
        TTreeReaderArray<float> trkPt = {fReader, "trkPt"};
        TTreeReaderArray<float> trkEta = {fReader, "trkEta"};
        TTreeReaderArray<float> trkPhi = {fReader, "trkPhi"};
        TTreeReaderArray<int> trkCharge = {fReader, "trkCharge"};
        TTreeReaderArray<float> trkDxy = {fReader, "trkDxy"};
        TTreeReaderArray<float> trkDz = {fReader, "trkDz"};
        TTreeReaderArray<float> trkExy = {fReader, "trkExy"};
        TTreeReaderArray<float> trkEz = {fReader, "trkEz"};
        TTreeReaderArray<int> trkNVHMuon = {fReader, "trkNVHMuon"};
        TTreeReaderArray<int> trkNVHPixel = {fReader, "trkNVHPixel"};
        TTreeReaderArray<int> trkNVHTracker = {fReader, "trkNVHTracker"};
        TTreeReaderArray<int> trkNVHAll = {fReader, "trkNVHAll"};
        TTreeReaderArray<int> trkLPPixel = {fReader, "trkLPPixel"};
        TTreeReaderArray<int> trkLPStrips = {fReader, "trkLPStrips"};
        TTreeReaderArray<int> trkLPTracker = {fReader, "trkLPTracker"};
        TTreeReaderArray<int> trkMIH = {fReader, "trkMIH"};
        TTreeReaderArray<float> trkGenMatchRadius = {fReader, "trkGenMatchRadius"};
        TTreeReaderArray<float> genPt = {fReader, "genPt"};
        TTreeReaderArray<float> genEta = {fReader, "genEta"};
        TTreeReaderArray<float> genPhi = {fReader, "genPhi"};
        TTreeReaderArray<float> genMass = {fReader, "genMass"};
        TTreeReaderArray<float> genCt = {fReader, "genCt"};
        TTreeReaderArray<int> genId = {fReader, "genId"};
        TTreeReaderArray<int> genCharge = {fReader, "genCharge"};
        
        RooArgSet row(mass, genMatchState);
        RooDataSet massData("massData", "", row);
        
        for (auto entry: fReader) {
            if (entry% 5000 == 0) {
                clog << "entry " << entry << '/' << nEvents << endl;
            }
            
            auto nSV = svtPt.GetSize();

            for (int iSV = 0; iSV < nSV; iSV++) {
                if (svtType[iSV] != 212) {
                    continue;
                }
                
                if (svtMass[svtJPsiIndex[iSV]] > 3.3 or svtMass[svtJPsiIndex[iSV]] < 2.9 ) {
                    continue;
                }
                
                if (svtMass[svtPsi2SIndex[iSV]] > 3.7 or svtMass[svtPsi2SIndex[iSV]] < 3.673 ) {
                    continue;
                }
                
                if (svtMass[iSV] < massLimitLow or svtMass[iSV] > massLimitHigh) {
//                    clog << svtType[iSV] << endl;
//                    clog << massLimitLow << " < " << svtMass[iSV] << " < " << massLimitHigh << endl;
                    continue;
                }
                
//                if (svtIsGoodMatch[iSV] < 3) {
//                    if (svtMass[iSV] > massCentLimitLow and svtMass[iSV] < massCentLimitHigh) {
//                        discardedBkgYield[name] += 1;
//                    }
//                    continue;
//                }
                
//                if (svtGenMatchRadius[svtJPsiIndex[iSV]] > muRadiusCut or svtGenMatchDpt[svtJPsiIndex[iSV]] > muPtDiff) 
//                {
//                    if (svtMass[iSV] > massCentLimitLow and svtMass[iSV] < massCentLimitHigh) {
//                        discardedBkgYield[name] += 1;
//                    }
//                    continue;
//                }
                
//                if (svtGenMatchRadius[svtPsi2SIndex[iSV]] > piRadiusCut or svtGenMatchDpt[svtPsi2SIndex[iSV]] > piPtDiff) 
//                {
//                    if (svtMass[iSV] > massCentLimitLow and svtMass[iSV] < massCentLimitHigh) {
//                        discardedBkgYield[name] += 1;
//                    }
//                    continue;
//                }
//                
//                if (svtGenMatchRadius[iSV] > kRadiusCut or svtGenMatchDpt[iSV] > kPtDiff) 
//                {
//                    if (svtMass[iSV] > massCentLimitLow and svtMass[iSV] < massCentLimitHigh) {
//                        discardedBkgYield[name] += 1;
//                    }
//                    continue;
//                }
                
                //svtMass[svtPsi2SIndex] < 3.74 && svtMass[svtPsi2SIndex] > 3.64 && svtMass[svtJPsiIndex] > 2.9 && svtMass[svtJPsiIndex] < 3.4
                
                hsvtMass[label]->Fill(svtMass[iSV]);
                
                mass.setVal(svtMass[iSV]);
                int matchState = 0;
                if (svtIsGoodMatch[iSV] >= 3){
                    matchState = 1;
                }
                genMatchState.setIndex(matchState);
                massData.add(row);
                
#define fillSvtHist(varname)\
                if (svtMass[iSV] < massSbLimitLow or svtMass[iSV] > massSbLimitHigh) h##varname##SB->Fill(varname[iSV]);\
                else h##varname[label]->Fill(varname[iSV]);
                
                fillSvtHist(svtPt)
                fillSvtHist(svtEta)
                
                fillSvtHist(svtCt)
                fillSvtHist(svtCtErr)
                
                fillSvtHist(svtDist3D);
                fillSvtHist(svtSigma3D);
                
                if (svtMass[iSV] < massSbLimitLow or svtMass[iSV] > massSbLimitHigh) hsvtSign3DSB->Fill(svtDist3D[iSV]/svtSigma3D[iSV]);
                else hsvtSign3D[label]->Fill(svtDist3D[iSV]/svtSigma3D[iSV]);
                
                if (svtMass[iSV] < massSbLimitLow or svtMass[iSV] > massSbLimitHigh) hsvtProbSB->Fill(TMath::Prob(svtChi2[iSV], svtNdof[iSV]));
                else hsvtProb[label]->Fill(TMath::Prob(svtChi2[iSV], svtNdof[iSV]));
                
#define fillSvtHistByProxy(varname, proxy)\
                if (svtMass[iSV] < massSbLimitLow or svtMass[iSV] > massSbLimitHigh) h##varname##proxy##SB->Fill(varname[svt##proxy##Index[iSV]]);\
                if (svtMass[iSV] < massCentLimitHigh and svtMass[iSV] > massCentLimitLow) h##varname##proxy[label]->Fill(varname[svt##proxy##Index[iSV]]);
                
                fillSvtHistByProxy(svtMass, JPsi)
                fillSvtHistByProxy(svtPt, JPsi)
                fillSvtHistByProxy(svtEta, JPsi)
                fillSvtHistByProxy(svtMass, Psi2S)
                fillSvtHistByProxy(svtPt, Psi2S)
                fillSvtHistByProxy(svtEta, Psi2S)
                
#define fillHistByTrkProxy(varname, part)\
                if (svtMass[iSV] < massSbLimitLow or svtMass[iSV] > massSbLimitHigh) {\
                    h##varname##part##AllSB->Fill(varname[svt##part##1Index[iSV]]);\
                    h##varname##part##AllSB->Fill(varname[svt##part##2Index[iSV]]);\
                }\
                if (svtMass[iSV] < massCentLimitHigh and svtMass[iSV] > massCentLimitLow) {\
                    h##varname##part##All[label]->Fill(varname[svt##part##1Index[iSV]]);\
                    h##varname##part##All[label]->Fill(varname[svt##part##2Index[iSV]]);\
                }\
                if (trkPt[svt##part##1Index[iSV]] > trkPt[svt##part##2Index[iSV]]) {\
                    if (svtMass[iSV] < massSbLimitLow or svtMass[iSV] > massSbLimitHigh) {\
                        h##varname##part##MaxSB->Fill(varname[svt##part##1Index[iSV]]);\
                        h##varname##part##MinSB->Fill(varname[svt##part##2Index[iSV]]);\
                    }\
                    if (svtMass[iSV] < massCentLimitHigh and svtMass[iSV] > massCentLimitLow) {\
                        h##varname##part##Max[label]->Fill(varname[svt##part##1Index[iSV]]);\
                        h##varname##part##Min[label]->Fill(varname[svt##part##2Index[iSV]]);\
                    }\
                }\
                else {\
                    if (svtMass[iSV] < massSbLimitLow or svtMass[iSV] > massSbLimitHigh) {\
                        h##varname##part##MinSB->Fill(varname[svt##part##1Index[iSV]]);\
                        h##varname##part##MaxSB->Fill(varname[svt##part##2Index[iSV]]);\
                    }\
                    if (svtMass[iSV] < massCentLimitHigh and svtMass[iSV] > massCentLimitLow) {\
                        h##varname##part##Min[label]->Fill(varname[svt##part##1Index[iSV]]);\
                        h##varname##part##Max[label]->Fill(varname[svt##part##2Index[iSV]]);\
                    }\
                }
                
                fillHistByTrkProxy(trkPt, Pi)
                fillHistByTrkProxy(trkEta, Pi)
                fillHistByTrkProxy(trkDxy, Pi)
                fillHistByTrkProxy(trkExy, Pi)
                fillHistByTrkProxy(trkDz, Pi)
                fillHistByTrkProxy(trkEz, Pi)
                fillHistByTrkProxy(trkNVHAll, Pi)
                fillHistByTrkProxy(trkNVHPixel, Pi)
                
                fillSvtHistByProxy(trkPt, K)
                fillSvtHistByProxy(trkEta, K)
                fillSvtHistByProxy(trkDxy, K)
                fillSvtHistByProxy(trkExy, K)
                fillSvtHistByProxy(trkDz, K)
                fillSvtHistByProxy(trkEz, K)
                fillSvtHistByProxy(trkNVHAll, K)
                fillSvtHistByProxy(trkNVHPixel, K)
                
                TVector3 Pi1P;
                Pi1P.SetPtEtaPhi(trkPt[svtPi1Index[iSV]], trkEta[svtPi1Index[iSV]], trkPhi[svtPi1Index[iSV]]);
                TVector3 Pi2P;
                Pi2P.SetPtEtaPhi(trkPt[svtPi2Index[iSV]], trkEta[svtPi2Index[iSV]], trkPhi[svtPi2Index[iSV]]);
                
                auto piangle = TMath::ACos(Pi1P.Unit().Dot(Pi2P.Unit()));
                
                if (svtMass[iSV] < massSbLimitLow or svtMass[iSV] > massSbLimitHigh) hPiAngleSB->Fill(piangle);
                else hPiAngle[label]->Fill(piangle);
            }
        }
        
        if (massData.sumEntries() == 0) {
            clog << "No Entries" << endl;
            return;
        }
        
        auto& mu = *new RooRealVar("mu", "mu", 5.27932, 5.26, 5.29);
        auto& lambda = *new RooRealVar("lambda", "lambda", 0.014, 0.0001, 0.1);
        auto& gamma = *new RooRealVar("gamma", "gamma", 0, -1, 1);
        auto& delta = *new RooRealVar("delta", "delta", 1, 0.5, 1.5);
        RooJohnson sigPdfJohn("sigPdf", "", mass, mu, lambda, gamma, delta);
        
        auto& mu2 = *new RooRealVar("mass_mu2", "mass_mu", 5.27932, 5.26, 5.29);
        auto& lambda2 = *new RooRealVar("mass_lambda2", "mass_lambda", 0.027, 0, 0.03);
        RooGaussian sigPdf1("sigPdf1", "", mass, mu, lambda);
        RooGaussian sigPdf2("sigPdf2", "", mass, mu2, lambda2);
        auto& sigFraction = *new RooRealVar("sigFrac", "", 0, 1);
        RooAbsPdf* sigPdf;
        sigPdf = &sigPdfJohn;
        
        auto& bkgPeakMu = *new RooRealVar("mass_bkgPeakMu", "", 5.2, 5.14, 5.39);
        auto& bkgPeakSigma = *new RooRealVar("mass_bkgPeakSigma", "", 0.1, 0.08, 0.4);
        RooGaussian bkgPeak("mass_bkgPeak", "", mass, bkgPeakMu, bkgPeakSigma);
        auto& bkgSlope = *new RooRealVar("bkgSlope", "", -5, -10, 10);
        RooExponential bkgExp("mass_bkgExp", "mass bkg", mass, bkgSlope);
        auto& bkgPeakFrac = *new RooRealVar("mass_bkgPeakFrac", "", 0.1, 0.0, 1);
        RooAddPdf bkgPdf("bkgPdf", "", RooArgList(bkgPeak, bkgExp), bkgPeakFrac);
        if (true or name == "original") {
            bkgPeakFrac.setVal(0);
            bkgPeakFrac.setConstant();
            bkgPeakMu.setConstant();
            bkgPeakSigma.setConstant();
        }
        
        auto& nSig = *new RooRealVar("nSig", "Number of Signal Events in SIGNAL MC", 7000, 0., 1.5*hsvtMass[label]->Integral());
        auto& nBkg = *new RooRealVar("nBkg", "Number of BG events in SIGNAL MC", 3./4*hsvtMass[label]->Integral(), 0., 1.5*hsvtMass[label]->Integral());
        
        
        
        RooAddPdf massPdf("massPdf", "Total pdf", RooArgList(*sigPdf, bkgPdf), RooArgList(nSig, nBkg));
        //RooAddPdf massPdf("massPdf", "Total pdf", RooArgList(*sigPdf), RooArgList(nSig));
        
        
        auto sbFitRes = bkgExp.fitTo(massData, RooFit::PrintLevel(-1), RooFit::Save(), RooFit::Range("sbLeft,sbRight"), RooFit::Offset(true));
        sbFitRes->Print("V");
        
//        nBkg.setVal(0);
//        nBkg.setConstant();
//        bkgSlope.setConstant();
        
        
        auto fitRes = massPdf.fitTo(massData, RooFit::PrintLevel(-1), RooFit::Save(), RooFit::Offset(true));
        
        auto frame = mass.frame(RooFit::Bins(50));
        massData.plotOn(frame, RooFit::Name("data"));
        massPdf.plotOn(frame, RooFit::Components("sigPdf"), RooFit::LineColor(kGreen), RooFit::LineStyle(kDashed));
        massPdf.plotOn(frame, RooFit::Components("bkgPdf"), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
        massPdf.plotOn(frame, RooFit::Name("fit"));
        auto lineLow = new TLine(massCentLimitLow, 0, massCentLimitLow, frame->GetMaximum()*0.75);
        lineLow->SetLineStyle(kDashed);
        lineLow->SetLineColor(kRed);
        lineLow->SetLineWidth(2);
        frame->addObject(lineLow, "L");
        auto lineHigh = new TLine(massCentLimitHigh, 0, massCentLimitHigh, frame->GetMaximum()*0.75);
        lineHigh->SetLineStyle(kDashed);
        lineHigh->SetLineColor(kRed);
        lineHigh->SetLineWidth(2);
        frame->addObject(lineHigh, "L");
        massPdf.paramOn(frame, RooFit::Layout(0.6, 0.98, 0.95));
        TCanvas cMassFit;
        roo_pulls(cMassFit, frame, "data", "fit");
        cMassFit.SaveAs(Form("plotSvt/" + baselabel + "/massFit%s.png", label));
        
        fitRes->Print("V");
        if (fitRes->covQual() != 3 or fitRes->statusCodeHistory(0) != 0 or fitRes->statusCodeHistory(1) != 0) {
            throw std::runtime_error("bad fit");
        }

        
        massPdf.plotOn(comparisonFrame, RooFit::LineColor(colors[names[name]]), RooFit::Components("sigPdf"), RooFit::Normalization(massData.sumEntries(), RooAbsReal::NumEvent), RooFit::Name(names[name]));
        dynamic_cast<TNamed*>(comparisonFrame->findObject(names[name]))->SetTitle(names[name]);
        sigPdf->plotOn(comparisonFrameNorm, RooFit::LineColor(colors[names[name]]), RooFit::Name(names[name]));
        dynamic_cast<TNamed*>(comparisonFrameNorm->findObject(names[name]))->SetTitle(names[name]);
        
        yieldVal[name] = make_pair(nSig.getVal(), nSig.getError());
        bkgyieldVal[name] = make_pair(nBkg.getVal(), nBkg.getError());
        muVal[name] = make_pair(mu.getVal(), mu.getError());
        //s1Val[name] = make_pair(lambda.getVal(), lambda.getError());
        s2Val[name] = make_pair(lambda2.getVal(), lambda2.getError());
        lambdaVal[name] = make_pair(lambda.getVal(), lambda.getError());
        gammaVal[name] = make_pair(gamma.getVal(), gamma.getError());
        deltaVal[name] = make_pair(delta.getVal(), delta.getError());
        fracVal[name] = make_pair(sigFraction.getVal(), sigFraction.getError());
        
        for (auto varObj: *massPdf.getParameters(massData)) {
            auto var = dynamic_cast<RooRealVar*>(varObj);
            vars[name][var->GetName()] = var;
        }
        res[name] = fitRes;
        
        auto bkgInSbInt = bkgPdf.createIntegral(mass, mass, "sbLeft,sbRight")->getVal();
        auto bkgInCentInt = bkgPdf.createIntegral(mass, mass, "center")->getVal();
        //clog << bkgInCentInt << endl;
        double bkgInCoreInt = bkgPdf.createIntegral(mass, mass, "core")->getVal();
        double sigInCoreInt = sigPdf->createIntegral(mass, mass, "core")->getVal();
        
        coreFracs[name] = make_pair(sigInCoreInt, bkgInCoreInt);
        
        auto sigInSbInt = sigPdf->createIntegral(mass, mass, "sbLeft,sbRight")->getVal();
        auto sigInCentInt = sigPdf->createIntegral(mass, mass, "center")->getVal();
        //clog << sigInCentInt << endl;
        
        auto bkgSbWeight = bkgInCentInt/bkgInSbInt;
        auto sigRescale = 1./(sigInCentInt-sigInSbInt*bkgInCentInt/bkgInSbInt);//1/(1-sigInSbInt/bkgInSbInt);
        
        
#define subSb(varname)\
        sbSub(h##varname[label], h##varname##SB, bkgSbWeight, sigRescale);
        
//        clog << "before " << hsvtPt[label]->GetEntries() << endl;
//        clog << "sb " << hsvtPtSB->GetEntries() << endl;
//        clog << nSig.getVal() << " " << nBkg.getVal() << " " << nSig.getVal() + nBkg.getVal() << endl;
//        clog << "bkg in sb " << nBkg.getVal()*bkgInSbInt << endl;
//        clog << "bkg in cent " << nBkg.getVal()*(1-bkgInSbInt) << endl;
//        clog << "sig in sb " << nSig.getVal()*sigInSbInt << endl;
//        clog << "sig in cent " << nSig.getVal()*(1-sigInSbInt) << endl;
        
//        subSb(svtPt)
//        subSb(svtEta)
//        
////        clog << "after " << hsvtPt[label]->Integral() << endl;
//        
//        subSb(svtCt)
//        subSb(svtCtErr)
//        
//        subSb(svtMassPsi2S)
//        subSb(svtPtPsi2S)
//        subSb(svtEtaPsi2S)
//        subSb(svtMassJPsi)
//        subSb(svtPtJPsi)
//        subSb(svtEtaJPsi)
//        
//        subSb(svtDist3D)
//        subSb(svtSign3D)
//        subSb(svtSigma3D)
//        
//        subSb(svtProb)
        
#define subSbTrk(varname, part)\
        sbSub(h##varname##part##All[label], h##varname##part##AllSB, bkgSbWeight, sigRescale);\
        sbSub(h##varname##part##Min[label], h##varname##part##MinSB, bkgSbWeight, sigRescale);\
        sbSub(h##varname##part##Max[label], h##varname##part##MaxSB, bkgSbWeight, sigRescale);
//        #define subSbTrk(varname, part)
        
//        clog << "before " << htrkPtPiMin[label]->GetEntries() << endl;
//        clog << "sb " << htrkPtPiMinSB->GetEntries() << endl;
//        subSbTrk(trkPt, Pi)
//        subSbTrk(trkEta, Pi)
//        subSbTrk(trkDxy, Pi)
//        subSbTrk(trkDz, Pi)
//        subSbTrk(trkExy, Pi)
//        subSbTrk(trkEz, Pi)
//        
//        subSb(trkPtK)
//        subSb(trkEtaK)
//        subSb(trkDxyK)
//        subSb(trkDzK)
//        subSb(trkExyK)
//        subSb(trkEzK)
        
//        clog << "after " << htrkPtPiMin[label]->Integral() << endl;
        
//        subSb(PiAngle)
    }
    
#define drawSvt(varname)\
    drawAll(#varname, h##varname, baselabel);

#define drawSvtLog(varname)\
    drawAll(#varname, h##varname, baselabel, true);
    
    drawSvt(svtMass)
    
    drawSvtLog(svtPt)
    drawSvt(svtEta)
    
    drawSvtLog(svtCt)
    drawSvt(svtCtErr)
    
    for (auto&& [name, hist]: hsvtMassJPsi) {
        hist->Scale(1./hist->Integral());
    }
    
    drawSvt(svtMassJPsi)
    drawSvtLog(svtPtJPsi)
    drawSvt(svtEtaJPsi)
    
    for (auto&& [name, hist]: hsvtMassPsi2S) {
        hist->Scale(1./hist->Integral());
    }
    drawSvt(svtMassPsi2S)
    drawSvtLog(svtPtPsi2S)
    drawSvt(svtEtaPsi2S)
    
    drawSvt(svtDist3D)
    drawSvt(svtSign3D)
    drawSvtLog(svtSigma3D)
    
    drawSvt(svtProb)

#define drawTrkLog(varname, part)\
    drawAll(#varname #part "All", h##varname##part##All, baselabel, true);\
    drawAll(#varname #part "Min", h##varname##part##Min, baselabel, true);\
    drawAll(#varname #part "Max", h##varname##part##Max, baselabel, true);

#define drawTrk(varname, part)\
    drawAll(#varname #part "All", h##varname##part##All, baselabel);\
    drawAll(#varname #part "Min", h##varname##part##Min, baselabel);\
    drawAll(#varname #part "Max", h##varname##part##Max, baselabel);
    
    drawTrkLog(trkPt, Pi)
    drawTrk(trkEta, Pi)
    drawTrk(trkDxy, Pi)
    drawTrkLog(trkExy, Pi)
    drawTrk(trkDz, Pi)
    drawTrkLog(trkEz, Pi)
    drawTrk(trkNVHPixel, Pi)
    drawTrk(trkNVHAll, Pi)
    
    drawSvtLog(trkPtK)
    drawSvt(trkEtaK)
    drawSvt(trkDxyK)
    drawSvtLog(trkExyK)
    drawSvt(trkDzK)
    drawSvtLog(trkEzK)
    drawSvt(trkNVHPixelK)
    drawSvt(trkNVHAllK)
    
    drawSvt(PiAngle)
    
    TCanvas cComp;
    comparisonFrame->SetYTitle("");
    comparisonFrame->Draw();
    auto legend = comparisonFrame->BuildLegend();
    legend->Draw();
    cComp.Update();
    legend->SetX1NDC(0.1);
    legend->SetX2NDC(0.4);
    legend->SetY1NDC(0.6);
    legend->SetY2NDC(0.9);
    cComp.Modified();
    cComp.SaveAs("plotSvt/" + baselabel + "/sigComparison.png");
    
    comparisonFrameNorm->Draw();
    comparisonFrameNorm->SetYTitle("");
    auto legendNorm = comparisonFrameNorm->BuildLegend();
    legendNorm->Draw();
    cComp.Update();
    legendNorm->SetX1NDC(0.1);
    legendNorm->SetX2NDC(0.4);
    legendNorm->SetY1NDC(0.6);
    legendNorm->SetY2NDC(0.9);
    cComp.Modified();
    cComp.SaveAs("plotSvt/" + baselabel + "/sigComparisonNorm.png");
    
    for (auto&& [name, label]: orderedNames) {
        clog << name << endl;
        clog << "\traw yields" << endl
            << "\t\tsignal " << *vars[name]["nSig"]->format(0, "NEXUP") << endl
            << "\t\tbkg " << *vars[name]["nBkg"]->format(0, "NEXUP") << endl;
        clog << "\tyields in window" << endl
            << "\t\tsignal " << vars[name]["nSig"]->getVal()*coreFracs[name].first << " \\pm " << vars[name]["nSig"]->getError()*coreFracs[name].first << endl
            << "\t\tbkg " << vars[name]["nBkg"]->getVal()*coreFracs[name].second << " \\pm " << vars[name]["nBkg"]->getError()*coreFracs[name].second << endl;
        
        RooFormulaVar significance("significance", Form("nSig*%f/(nBkg*%f)", coreFracs[name].first, coreFracs[name].second), 
            RooArgList(
                *vars[name]["nSig"],
                *vars[name]["nBkg"]
            )
        );
        RooFormulaVar significanceSqr("significanceSqr", Form("pow(nSig*%1$f, 2)/(nBkg*%2$f + nSig*%1$f)", coreFracs[name].first, coreFracs[name].second), 
            RooArgList(
                *vars[name]["nSig"],
                *vars[name]["nBkg"]
            )
        );
        clog << "\tsignificance" << endl
            << "\t\tS/B " << significance.getVal() << " \\pm " << significance.getPropagatedError(*res[name]) << endl
            << "\t\tS^2/(S+B) " << significanceSqr.getVal() << " \\pm " << significanceSqr.getPropagatedError(*res[name]) << endl;
        
        RooFormulaVar effSigma("effSigma", "lambda*sqrt((exp(pow(delta, -2)) - 1)*(exp(pow(delta, -2))*cosh(2*gamma/delta) + 1)/2)", 
            RooArgList(
                *vars[name]["lambda"],
                *vars[name]["delta"],
                *vars[name]["gamma"]
            )
        );
        clog << "\tsigma" << endl
            << "\t\tlambda " << *vars[name]["lambda"]->format(0, "NEXUP") << endl
            << "\t\teff sigma " << effSigma.getVal() << " \\pm " << effSigma.getPropagatedError(*res[name]) << endl;
        
//        clog << name << " "
//            << yieldVal[name].first << "+/-" << yieldVal[name].second << " " 
//            << bkgyieldVal[name].first << "+/-" << bkgyieldVal[name].second << " "
//            << yieldVal[name].first/bkgyieldVal[name].first << " "
//            << yieldVal[name].first*coreFracs[name].first << "+/-" << yieldVal[name].second*coreFracs[name].first << " " 
//            << bkgyieldVal[name].first*coreFracs[name].second + discardedBkgYield[name] << "+/-" << bkgyieldVal[name].second*coreFracs[name].second << " "
//            << yieldVal[name].first*coreFracs[name].first/(bkgyieldVal[name].first*coreFracs[name].second + discardedBkgYield[name]) << " "
//            << muVal[name].first << "+/-" << muVal[name].second << " " 
//            << lambdaVal[name].first << "+/-" << lambdaVal[name].second << " " 
//            << gammaVal[name].first << "+/-" << gammaVal[name].second << " " 
//            << deltaVal[name].first << "+/-" << deltaVal[name].second << endl;
    }
    
    auto prevDir = gDirectory;
    auto massFile = TFile::Open("massHists.root", "RECREATE");
    for (auto&& [name, hist]: hsvtMass) {
        hist->Write();
    }
    delete massFile;
    prevDir->cd();
}
