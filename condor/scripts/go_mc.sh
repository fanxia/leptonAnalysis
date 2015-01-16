#!/bin/bash

JOB_NUMBER=$1
STAGING=$2
WORK_DIR=`pwd`

export CONDOR_SECTION=$1

tar -xzf fileLists.tgz

#export VO_CMS_SW_DIR=/uscmst1/prod/sw/cms
#source $VO_CMS_SW_DIR/cmsset_default.sh
source /cvmfs/cms.cern.ch/cmsset_default.csh
export SCRAM_ARCH=slc5_amd64_gcc462
scramv1 project CMSSW CMSSW_5_3_8_patch3

cd CMSSW_5_3_8_patch3/src/
mv $WORK_DIR/src.tgz .
tar -xzf src.tgz
eval `scramv1 runtime -sh`

mv $WORK_DIR/ANALYZER .
mv $WORK_DIR/filelist_$JOB_NUMBER .

if [ "$STAGING" == "acceptance" ];
then
	mv $WORK_DIR/pileupReweighting_DATASETNAME.root .
	mv $WORK_DIR/btagEfficiency_DATASETNAME.root .
fi

while read file
do

        if [[ $file == /pnfs/* ]];
        then
            	sed -i '13i chain.Add("dcap://'$file'");' ANALYZER
        else
            	sed -i '13i chain.Add("'$file'");' ANALYZER
        fi

done < filelist_$JOB_NUMBER

sed -i "s/STAGING/$STAGING/g" ANALYZER

root -b -q -l ANALYZER

if [ "$STAGING" == "pileup" ];
then
	mv pileupReweighting_DATASETNAME.root $WORK_DIR/pileupReweighting_DATASETNAME_$JOB_NUMBER.root
elif [ "$STAGING" == "btag" ];
then
	mv btagEfficiency_DATASETNAME.root $WORK_DIR/btagEfficiency_DATASETNAME_$JOB_NUMBER.root
elif [ "$STAGING" == "acceptance" ];
then
	mv signal_contamination_DATASETNAME.root $WORK_DIR/signal_contamination_DATASETNAME_$JOB_NUMBER.root
else
	mv *.root $WORK_DIR
fi

rm pileup_22Jan2013_69400.root
rm pileup_22Jan2013_72870.root
rm pileup_22Jan2013_65930.root

cd $WORK_DIR
rm -rf CMSSW_5_3_8_patch3/
rm fileLists.tgz
rm filelist_*
