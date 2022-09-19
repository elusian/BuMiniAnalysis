#!/bin/sh

export DIR=`dirname $0`
export TOP=`dirname ${DIR}`

${DIR}/changeReturnType.sh bool void ${TOP}/Write/interface/NtuFilter.h
${DIR}/changeReturnType.sh bool void ${TOP}/Write/interface/NtuFilter.hpp

${DIR}/changeArgumentType.sh "      edm::Run" "const edm::Run" ${TOP}/Write/interface/NtuFilter.h
${DIR}/changeArgumentType.sh "      edm::Run" "const edm::Run" ${TOP}/Write/interface/NtuFilter.hpp

eval `find ${TOP} | grep BuildFile | grep 'xml$' | awk -v DIR=${DIR} '{print DIR"/cleanupBuildFile.sh "$0" ;"}'`
