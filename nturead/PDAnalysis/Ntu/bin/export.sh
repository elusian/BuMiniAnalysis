#!/bin/sh

cd `dirname $0`/../../..
export PD_PATH=PDAnalysis/Ntu
export PD_PACK=`pwd`/pdAnalysis.tgz
rm -f ${PD_PACK}
tar -czf ${PD_PACK}                    \
NtuAnalysis/Common                      \
NtuAnalysis/Read                        \
${PD_PATH}/interface/*.h               \
${PD_PATH}/interface/*.hpp             \
${PD_PATH}/src/*.cc                    \
${PD_PATH}/bin/*.h                     \
${PD_PATH}/bin/*.cc                    \
${PD_PATH}/bin/compile.sh              \
${PD_PATH}/bin/*Macro.C                \
${PD_PATH}/bin/*Proof.mac              \
${PD_PATH}/bin/setupMacroAnalysis.sh   \
${PD_PATH}/bin/installProofAnalysis.sh \
${PD_PATH}/bin/*_test.list             \
${PD_PATH}/bin/export.sh
echo ${PD_PACK}
