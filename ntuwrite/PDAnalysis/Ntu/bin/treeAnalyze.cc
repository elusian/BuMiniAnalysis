#define UTIL_USE FULL
#include "PDAnalyzer.h"
#include "NtuAnalysis/Read/interface/NtuFlexibleAnalyzer.h"
#include "NtuTool/Read/bin/treeAnalyze.cc"
static NtuFlexibleAnalyzer<PDAnalyzer> nfa( "pdAnalysis", "pdAnalyzer" );
