#!/bin/sh

export TYPE_I=$1
export TYPE_F=$2
export FILE=$3

export FTMP=/tmp/`basename ${FILE}``date +%s`

awk -v TYPE_I=${TYPE_I} -v TYPE_F=${TYPE_F} '(index($0,"beginRun")==0){print $0} (index($0,"beginRun")>0){line=$0; sub(TYPE_I,TYPE_F,line); print line}' ${FILE} | awk -v TYPE_I=${TYPE_I} -v TYPE_F=${TYPE_F} '(index($0,"endRun")==0){print $0} (index($0,"endRun")>0){line=$0; sub(TYPE_I,TYPE_F,line); print line}' | awk -v TYPE_F=${TYPE_F} 'BEGIN{ff=0} ((ff==0)||(index($0,"return")==0)){print $0} (index($0,"beginRun")>0){ff=1} ((ff==1)&&(index($0,"return")>0)){printf(substr($0,1,index($0,"return")+5)); if(TYPE_F=="bool") printf(" true"); print ";"; ff=0}' | awk -v TYPE_F=${TYPE_F} 'BEGIN{ff=0} ((ff==0)||(index($0,"return")==0)){print $0} (index($0,"endRun")>0){ff=1} ((ff==1)&&(index($0,"return")>0)){printf(substr($0,1,index($0,"return")+5)); if(TYPE_F=="bool") printf(" true"); print ";"; ff=0}' > ${FTMP}
cp ${FTMP} ${FILE}
rm ${FTMP}
