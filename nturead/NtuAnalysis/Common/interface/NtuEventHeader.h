#ifndef NtuEventHeader_H
#define NtuEventHeader_H

class TBranch;

class NtuEventHeader {

 protected:

  NtuEventHeader() { b_runNumber   =
                     b_lumiSection =
                     b_eventNumber = 0; }

  // header
  unsigned int runNumber;
  unsigned int lumiSection;
  unsigned int eventNumber;

  TBranch* b_runNumber;
  TBranch* b_lumiSection;
  TBranch* b_eventNumber;

 private:

  // dummy copy constructor and assignment
  NtuEventHeader           ( const NtuEventHeader& );
  NtuEventHeader& operator=( const NtuEventHeader& );

};

#endif

