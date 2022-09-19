#ifndef PDAnalyzer_H
#define PDAnalyzer_H

#include "TH1.h"
#include "PDAnalyzerUtil.h"
#include "PDAnalysis/Ntu/interface/PDGenHandler.h"
#include "Math/Vector3D.h"

#include <fstream>
#include <unordered_map>
#include <array>
#include <optional>

class PDSecondNtupleWriter;

class PDAnalyzer: public virtual PDAnalyzerUtil, public virtual PDGenHandler {

 public:

  PDAnalyzer();
  virtual ~PDAnalyzer();

  // function called before starting the analysis
  virtual void beginJob();

  // functions to book the histograms
  void book();

  // functions called for each event
  // function to reset class content before reading from file
  virtual void reset();
  // function to do event-by-event analysis,
  // return value "true" for accepted events
  virtual bool analyze( int entry, int event_file, int event_tot );

  // function called at the end of the analysis
  virtual void endJob();

  bool verbose;

 private:

  TH1D* hptmumax;
  TH1D* hptmuall;

// additional features: second ntuple
  PDSecondNtupleWriter* tWriter;
  
  float getCt(int iGen);
  int findPV(ROOT::Math::XYZVector sv, ROOT::Math::XYZVector bP);

  // dummy copy constructor and assignment
  PDAnalyzer           ( const PDAnalyzer& );
  PDAnalyzer& operator=( const PDAnalyzer& );
  
  bool hasDaughter(int iGen);
  void printDecayChain(int iGen, const std::string& pre = "");
  void printCompactDecayChain(int iGen, std::ostream& out);
  
  int GetClosestGen( float eta, float phi, float pt );
  std::pair<int, float> GetClosestTrackGen( float eta, float phi, float pt, int charge );
  std::pair<int, float> GetClosestTrackGenFromColl( float eta, float phi, float pt, int charge, const std::vector<int>& coll);
  std::optional<std::array<int, 6>> getBuVertexComponents(int iBuGen);
  std::optional<std::array<int, 5>> getPsi2SVertexComponents(int iPsi2SGen);
  std::optional<std::array<int, 3>> getJPsiVertexComponents(int iJPsiGen);

  bool isTrkHighPurity(int itk){ return (( trkQuality->at( itk ) >> 2 ) & 1); }
  
  std::ostream* decayOut;
  std::unordered_map<std::string, int> decays;

};


#endif

