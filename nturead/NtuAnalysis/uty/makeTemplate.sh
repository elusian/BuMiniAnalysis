#!/bin/sh

if [ -z ${CMSSW_BASE} ] ; then
  echo "CMSSW_BASE not defined; run \"cmsenv\" before"
  exit
fi

export CURRENTDIR=`/bin/pwd`

export PREFIX=`echo $1 | awk '{print toupper($0)}'`
export LOWFIX=`echo $1 | awk '{print tolower($0)}'`

export PRETMP="TMP"
export LOWTMP="tmp"
export SUBSYSTEM=${CMSSW_BASE}/src/${PREFIX}Analysis
export SUBSYSTMP=${CMSSW_BASE}/src/NtuAnalysis/Template
#export SUBSYSTMP=${CMSSW_BASE}/src/${PRETMP}Analysis

export AnalyzerUtil=${PREFIX}AnalyzerUtil
export AnalyzerFW=${PREFIX}AnalyzerFW
export Analyzer=${PREFIX}Analyzer
export NtupleData=${PREFIX}NtupleData
export NtupleBranch=${PREFIX}NtupleBranch
export EDMToNtuple=${PREFIX}EDMToNtuple
export Ntuplizer=${PREFIX}Ntuplizer
export EDMNtuplizer=${PREFIX}EDMNtuplizer
export Filter=${PREFIX}Filter

create_file() {
cat $1 |\
    sed s/${PRETMP}tree/${PREFIX}tree/g |\
    sed s/${PRETMP}Analysis/${PREFIX}Analysis/g |\
    sed s/${PRETMP}AnalyzerUtil/${PREFIX}AnalyzerUtil/g |\
    sed s/${PRETMP}AnalyzerFW/${PREFIX}AnalyzerFW/g |\
    sed s/${PRETMP}Analyzer/${PREFIX}Analyzer/g |\
    sed s/${PRETMP}MacroAnalyzer/${PREFIX}MacroAnalyzer/g |\
    sed s/${PRETMP}EDMAnalyzer/${PREFIX}EDMAnalyzer/g |\
    sed s/${PRETMP}NtupleData/${PREFIX}NtupleData/g |\
    sed s/${PRETMP}NtupleBranch/${PREFIX}NtupleBranch/g |\
    sed s/${PRETMP}EDMToNtuple/${PREFIX}EDMToNtuple/g |\
    sed s/${PRETMP}Ntuplizer/${PREFIX}Ntuplizer/g |\
    sed s/${PRETMP}EDMNtuplizer/${PREFIX}EDMNtuplizer/g |\
    sed s/${PRETMP}Proof/${PREFIX}Proof/g |\
    sed s/${PRETMP}Filter/${PREFIX}Filter/g |\
    sed s/${PRETMP}_PACK/${PREFIX}_PACK/g |\
    sed s/${PRETMP}_PATH/${PREFIX}_PATH/g |\
    sed s/${PRETMP}_ROOT_DIR/${PREFIX}_ROOT_DIR/g |\
    sed s/${PRETMP}_EXE_FILE/${PREFIX}_EXE_FILE/g |\
    sed s/${LOWTMP}Macro/${LOWFIX}Macro/g |\
    sed s/${LOWTMP}Proof/${LOWFIX}Proof/g |\
    sed s/${LOWTMP}Analysis/${LOWFIX}Analysis/g |\
    sed s/${LOWTMP}Analyzer/${LOWFIX}Analyzer/g |\
    sed s/${LOWTMP}Select/${LOWFIX}Select/g |\
    sed s/${LOWTMP}Filter/${LOWFIX}Filter/g |\
    sed s/${LOWTMP}FWLite/${LOWFIX}FWLite/g |\
    sed s/${LOWTMP}TreeAnalyze/${LOWFIX}TreeAnalyze/g >\
    $2
}

create_name() {
create_file ${SUBSYSTMP}/$1/${PRETMP}$2 ${SUBSYSTEM}/$1/${PREFIX}$2
}

create_nlow() {
create_file ${SUBSYSTMP}/$1/${LOWTMP}$2 ${SUBSYSTEM}/$1/${LOWFIX}$2
}

create_same() {
create_file ${SUBSYSTMP}/$1 ${SUBSYSTEM}/$1
}

create_temp() {
mv ${SUBSYSTEM}/Ntu/bin/${PREFIX}$1 ${SUBSYSTEM}/Ntu/bin/${PRETMP}$1
sed s#TMPNTUANADIR#${SUBSYSTEM}/Ntu/bin# \
    ${SUBSYSTEM}/Ntu/bin/${PRETMP}$1 >\
    ${SUBSYSTEM}/Ntu/bin/${PREFIX}$1
rm -f ${SUBSYSTEM}/Ntu/bin/${PRETMP}$1
}

create_tlow() {
mv ${SUBSYSTEM}/Ntu/bin/${LOWFIX}$1 ${SUBSYSTEM}/Ntu/bin/${LOWTMP}$1
sed s#TMPNTUANADIR#${SUBSYSTEM}/Ntu/bin# \
    ${SUBSYSTEM}/Ntu/bin/${LOWTMP}$1 >\
    ${SUBSYSTEM}/Ntu/bin/${LOWFIX}$1
rm -f ${SUBSYSTEM}/Ntu/bin/${LOWTMP}$1
}

mkdir -p ${SUBSYSTEM}/Ntu
create_same Ntu/BuildFile.xml

mkdir -p ${SUBSYSTEM}/Ntu/interface
create_name Ntu/interface NtupleData.h
create_name Ntu/interface NtupleBranch.h
create_name Ntu/interface NtupleBranch.hpp

mkdir -p ${SUBSYSTEM}/Ntu/src
create_name Ntu/src NtupleData.cc
create_name Ntu/src AnalyzerFW.cc

mkdir -p ${SUBSYSTEM}/Ntu/bin
create_same Ntu/bin/BuildFile.xml
create_name Ntu/bin Analyzer.h
create_name Ntu/bin AnalyzerUtil.h
create_name Ntu/bin Analyzer.cc
create_name Ntu/bin AnalyzerUtil.cc
create_name Ntu/bin ProofSetup.cc
create_name Ntu/bin ProofAnalyzer.h
create_name Ntu/bin ProofAnalyzer.cc
create_name Ntu/bin ProofLightAnalyzer.h
create_name Ntu/bin ProofLightAnalyzer.cc
create_nlow Ntu/bin Macro.C
create_nlow Ntu/bin Proof.mac
create_same Ntu/bin/treeAnalyze.cc
create_same Ntu/bin/compile.sh
create_same Ntu/bin/export.sh
create_same Ntu/bin/setupMacroAnalysis.sh
cat ${SUBSYSTMP}/Ntu/bin/installProofAnalysis.sh |\
    sed s/PROOFPREFIX/${PREFIX}/ |\
    sed s/PROOFPRETMP/${PRETMP}/ |\
    sed s/${PRETMP}Analysis/${PREFIX}Analysis/ >\
    ${SUBSYSTEM}/Ntu/bin/installProofAnalysis.sh
chmod u+x ${SUBSYSTEM}/Ntu/bin/*.sh
cp ${SUBSYSTMP}/Ntu/bin/*list ${SUBSYSTEM}/Ntu/bin

create_temp ProofAnalyzer.cc
create_temp ProofLightAnalyzer.cc
create_tlow Proof.mac

mkdir -p ${SUBSYSTEM}/EDM
create_same EDM/BuildFile.xml

mkdir -p ${SUBSYSTEM}/EDM/interface
create_name EDM/interface EDMToNtuple.h

mkdir -p ${SUBSYSTEM}/EDM/src
create_name EDM/src EDMToNtuple.cc

mkdir -p ${SUBSYSTEM}/EDM/bin
create_same EDM/bin/BuildFile.xml
create_nlow EDM/bin FWLite.cc
create_same EDM/bin/cfg_fwlite.py
create_same EDM/bin/cfg_fwfull.py
create_same EDM/bin/cfg_fwfEDM.py
create_same EDM/bin/cfg_filter.py

mkdir -p ${SUBSYSTEM}/EDM/plugins
create_same EDM/plugins/BuildFile.xml
create_name EDM/plugins Ntuplizer.cc
create_name EDM/plugins EDMNtuplizer.cc
create_name EDM/plugins Filter.cc

