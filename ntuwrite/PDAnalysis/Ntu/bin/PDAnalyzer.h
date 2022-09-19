#ifndef PDAnalyzer_H
#define PDAnalyzer_H

#include "TH1.h"
#include "PDAnalyzerUtil.h"

class PDAnalyzer: public virtual PDAnalyzerUtil {

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

  // dummy copy constructor and assignment
  PDAnalyzer           ( const PDAnalyzer& );
  PDAnalyzer& operator=( const PDAnalyzer& );

};


#endif

