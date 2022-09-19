#define UTIL_USE FULL
#include "TMPAnalyzer.h"
#include "NtuAnalysis/Read/interface/NtuFlexibleAnalyzer.h"
#include "NtuTool/Read/bin/treeAnalyze.cc"
static NtuFlexibleAnalyzer<TMPAnalyzer> nfa( "tmpAnalysis", "tmpAnalyzer" );
