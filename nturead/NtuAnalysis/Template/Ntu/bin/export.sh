#!/bin/sh

cd `dirname $0`/../../..
export TMP_PATH=TMPAnalysis/Ntu
export TMP_PACK=`pwd`/tmpAnalysis.tgz
rm -f ${TMP_PACK}
tar -czf ${TMP_PACK}                    \
NtuAnalysis/Common                      \
NtuAnalysis/Read                        \
${TMP_PATH}/interface/*.h               \
${TMP_PATH}/interface/*.hpp             \
${TMP_PATH}/src/*.cc                    \
${TMP_PATH}/bin/*.h                     \
${TMP_PATH}/bin/*.cc                    \
${TMP_PATH}/bin/compile.sh              \
${TMP_PATH}/bin/*Macro.C                \
${TMP_PATH}/bin/*Proof.mac              \
${TMP_PATH}/bin/setupMacroAnalysis.sh   \
${TMP_PATH}/bin/installProofAnalysis.sh \
${TMP_PATH}/bin/*_test.list             \
${TMP_PATH}/bin/export.sh
echo ${TMP_PACK}
