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
#include "NtuTool/Read/interface/TreeReader.h"
#include "NtuTool/Read/interface/EDMTreeReader.h"

#include "PDProofSetup.cc"

class PDProofReaderInterface: public PDAnalyzer {
 public:
  void bJ() {
    this->setupNtuple();
    this->beginJob();
    return;
  }
  bool gE( int ientry ) {
    return this->getEntry( ientry );
  };
};

class PDProofReader: public PDProofReaderInterface,
                      public TreeReader {
 protected:
  void iT( TTree* tree ) {
    cout << "initRead" << endl;
    this->initRead( tree );
  };
};

class PDProofEDMReader: public PDProofReaderInterface,
                         public EDMTreeReader {
 public:
  PDProofEDMReader() {
    setProcessName ( processName  );
    setProducerName( producerName );
  }
 protected:
  void iT( TTree* tree ) {
    cout << "initRead" << endl;
    this->initRead( tree );
  };
};


class PDProofAnalyzer: public TSelector,
                        public NtuPReader {

 public:

  PDProofAnalyzer( TTree* tree = 0 );
  virtual ~PDProofAnalyzer();
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

  ClassDef(PDProofAnalyzer,0);

};

#endif

