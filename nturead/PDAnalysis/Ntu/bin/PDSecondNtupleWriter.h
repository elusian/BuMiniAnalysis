#ifndef PDSecondNtupleWriter_h
#define PDSecondNtupleWriter_h

#include "PDSecondNtupleData.h"
#include "NtuTool/Common/interface/TreeWriter.h"
#include "TFile.h"
#include <string>

class PDSecondNtupleWriter: public PDSecondNtupleData, public TreeWriter {

 public:

  PDSecondNtupleWriter() {
  }
  virtual ~PDSecondNtupleWriter() {
  }

  void open( const std::string& name, const std::string mode = "CREATE" ) {
    TDirectory* current = gDirectory;
    file = new TFile( name.c_str(), mode.c_str() );
    initTree();
    initWrite( file );
    current->cd();
    return;
  }
 
  void close() {
    TreeWriter::close();
    file->Close();
  }

 private:

  TFile* file;

  PDSecondNtupleWriter           ( const PDSecondNtupleWriter& a );
  PDSecondNtupleWriter& operator=( const PDSecondNtupleWriter& a );

};

#endif

