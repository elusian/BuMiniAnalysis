#!/bin/sh

cd `dirname $0`/..
export PD_PACK_DIR=`/bin/pwd`
cd ../..
export PD_ROOT_DIR=`/bin/pwd`
echo $PD_PACK_DIR
echo $PD_ROOT_DIR

export SA_STRING="==="`printenv | awk -F= '($1=="SCRAM_ARCH") {print $2}'`"==="
if [ ${SA_STRING} = "======" ]; then
  export NTU_TOOL_LIB="-L "${NTU_TOOL_DIR}"/lib -lNtupleTool"
else
  export NTU_TOOL_LIB="-L "${NTU_TOOL_DIR}"/lib/"${SCRAM_ARCH}" -lNtupleTool_"${CMSSW_VERSION}
fi

export PD_EXE_FILE=${PD_PACK_DIR}/bin/pdTreeAnalyze
rm -f ${PD_EXE_FILE}

c++ -Wall `root-config --cflags`              \
-I ${PD_ROOT_DIR} -I ${NTU_TOOL_DIR}/include \
-o ${PD_EXE_FILE}                            \
`ls -1 $PD_PACK_DIR/src/*cc |\
       awk '{printf($0" ")}'`                 \
${PD_PACK_DIR}/bin/treeAnalyze.cc            \
${PD_ROOT_DIR}/NtuAnalysis/*/src/*cc         \
${NTU_TOOL_LIB}                               \
`root-config --glibs`

echo ${PD_EXE_FILE}
