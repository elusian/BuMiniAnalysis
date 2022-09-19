#ifndef PDProofAnalyzer_h
#define PDProofAnalyzer_h

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TSelector.h"
#include "TProofServ.h"
#include "TH1.h"
#include <string>
#include <vector>

#include "PDAnalyzer.h"

#include "PDProofSetup.cc"

class ProofLightReaderInterface: public PDAnalyzer {
 public:
  void bJ() {
    this->setupNtuple();
    this->beginJob();
    return;
  }
  bool gE( int ientry ) {
    return this->getEntry( ientry );
  };
  void process( int ientry ) {
    return;
  };
  void iT( TTree* tree ) {
    cout << "openNtuple" << endl;
    this->openTree( tree );
  }
};

class PDProofLightAnalyzer: public TSelector,
                             public ProofLightReaderInterface {

 public:

  PDProofLightAnalyzer( TTree* tree = 0 );
  virtual ~PDProofLightAnalyzer();
  virtual Int_t   Version() const { return 2; }
  virtual void    Begin( TTree *tree );
  virtual void    SlaveBegin( TTree *tree );
  virtual void    Init( TTree *tree );
  virtual Bool_t  Notify();
  virtual Bool_t  Process(Long64_t entry);
  virtual void    SetOption(const char *option) { fOption = option; }
  virtual void    SetObject(TObject *obj) { fObject = obj; }
  virtual void    SetInputList(TList *input) { fInput = input; }
  virtual TList  *GetOutputList() const { return fOutput; }
  virtual void    SlaveTerminate();
  virtual void    Terminate();

  ClassDef(PDProofLightAnalyzer,0);

};

#endif

