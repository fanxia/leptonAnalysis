#!/bin/bash

WORK_DIR=`pwd`

CHANGED=false

while read file
do

    if [[ ../$file -nt ../src.tgz ]];
    then
	CHANGED=true
	break
    fi

done < doc/listOfImportantFiles.txt

if [[ "$CHANGED" = true ]];
then
    
    export VO_CMS_SW_DIR=/uscmst1/prod/sw/cms
    source $VO_CMS_SW_DIR/cmsset_default.sh
    export SCRAM_ARCH=slc5_amd64_gcc462
    
    cd ../
    eval `scramv1 runtime -sh`
    root -b -q -l compileAnalyzer.C

    tar -czf src.tgz libSusyEvent.so SusyEventAnalyzer_cc.so jan3_pileup.root Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt rootlogon.C
fi