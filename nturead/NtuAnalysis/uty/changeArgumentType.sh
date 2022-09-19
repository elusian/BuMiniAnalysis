#!/bin/sh

export TYPE_I="$1"
export TYPE_F="$2"
export FILE=$3

export FTMP=/tmp/`basename ${FILE}``date +%s`

awk -v TYPE_I="${TYPE_I}" -v TYPE_F="${TYPE_F}" '(index($0,"beginRun")==0){print $0} (index($0,"beginRun")>0){line=$0; sub(TYPE_I,TYPE_F,line); print line;}' ${FILE} | awk -v TYPE_I="${TYPE_I}" -v TYPE_F="${TYPE_F}" '(index($0,"endRun")==0){print $0} (index($0,"endRun")>0){line=$0; sub(TYPE_I,TYPE_F,line); print line;}' > ${FTMP}

cp ${FTMP} ${FILE}
rm ${FTMP}
