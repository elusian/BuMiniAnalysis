#!/bin/sh

export DIR=`dirname $0`
export TOP=`dirname ${DIR}`

${DIR}/changeReturnType.sh void bool ${TOP}/Write/interface/NtuFilter.h
${DIR}/changeReturnType.sh void bool ${TOP}/Write/interface/NtuFilter.hpp

${DIR}/changeArgumentType.sh "const edm::Run" "      edm::Run" ${TOP}/Write/interface/NtuFilter.h
${DIR}/changeArgumentType.sh "const edm::Run" "      edm::Run" ${TOP}/Write/interface/NtuFilter.hpp

