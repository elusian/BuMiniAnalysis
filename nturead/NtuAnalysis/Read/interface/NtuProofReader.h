#ifndef NtuProofReader_h
#define NtuProofReader_h

//template<class T>
class NtuProofReaderInterface: public NtuAnalyzer {
 public:
  void bJ() {
    this->setupNtuple();
    this->beginJob();
    return;
  }
  bool gE( int ientry ) {
    return this->getEntry( ientry );
  };
  ClassDef(NtuProofReaderInterface,0);
};

//template <class T>
class NtuProofReader: public NtuProofReaderInterface,
                      public TreeReader {
 protected:
  void iT( TTree* tree ) {
    cout << "initRead" << endl;
    this->initRead( tree );
  };
  ClassDef(NtuProofReader,0);
};

//template <class T>
class NtuProofEDMReader: public NtuProofReaderInterface,
                         public EDMTreeReader {
 public:
  NtuProofEDMReader() {
    setProcessName ( processName  );
    setProducerName( producerName );
  }
 protected:
  void iT( TTree* tree ) {
    cout << "initRead" << endl;
    this->initRead( tree );
  };
  ClassDef(NtuProofEDMReader,0);
};


#endif

