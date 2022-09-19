#!/bin/sh

export WORKDIR=`dirname $0`
cd ${WORKDIR}
ln -sf ../../../NtuAnalysis
ln -sf ../../../PDAnalysis

if [ ! -s "${NTU_TOOL_DIR}" ]; then
  echo "NTU_TOOL_DIR not defined, only light version available"
fi

ln -sf ${NTU_TOOL_DIR}/NtuTool

export PREFIX="PD"
export PRETMP="TMP"
export LOWFIX=`echo ${PREFIX} | awk '{print tolower($0)}'`
export LOWTMP=`echo ${PRETMP} | awk '{print tolower($0)}'`

create_proof() {
rm ${PREFIX}$1
sed s#TMPNTUANADIR#${WORKDIR}# ${PRETMP}$1 > ${PREFIX}$1
}

create_prlow() {
rm ${LOWFIX}$1
sed s#TMPNTUANADIR#${WORKDIR}# ${LOWTMP}$1 > ${LOWFIX}$1
}

create_prlow Proof.mac
create_proof ProofAnalyzer.cc
create_proof ProofLightAnalyzer.cc


