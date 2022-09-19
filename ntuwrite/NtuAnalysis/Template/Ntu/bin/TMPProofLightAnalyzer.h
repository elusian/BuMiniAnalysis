#ifndef TMPProofAnalyzer_h
#define TMPProofAnalyzer_h

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TSelector.h"
#include "TProofServ.h"
#include "TH1.h"
#include <string>
#include <vector>

#include "TMPAnalyzer.h"

#include "TMPProofSetup.cc"

class ProofLightReaderInterface: public TMPAnalyzer {
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

class TMPProofLightAnalyzer: public TSelector,
                             public ProofLightReaderInterface {

 public:

  TMPProofLightAnalyzer( TTree* tree = 0 );
  virtual ~TMPProofLightAnalyzer();
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

  ClassDef(TMPProofLightAnalyzer,0);

};

#endif

