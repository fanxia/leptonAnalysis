#!/bin/bash

JOB_NUMBER=$1

WORK_DIR=`pwd`

export CONDOR_SECTION=$1

export VO_CMS_SW_DIR=/uscmst1/prod/sw/cms
source $VO_CMS_SW_DIR/cmsset_default.sh
export SCRAM_ARCH=slc5_amd64_gcc462
scramv1 project CMSSW CMSSW_5_3_8_patch3

cd CMSSW_5_3_8_patch3/src/
mv $WORK_DIR/src.tgz .
tar -xzf src.tgz
eval `scramv1 runtime -sh`

mv $WORK_DIR/ana_signal.C .

root -b -q -l ana_signal.C

mv *.root $WORK_DIR

rm pileup_22Jan2013_69400.root
rm pileup_22Jan2013_72870.root
rm pileup_22Jan2013_65930.root

cd $WORK_DIR
rm -rf CMSSW_5_3_8_patch3/
