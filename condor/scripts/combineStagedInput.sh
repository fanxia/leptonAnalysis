#!/bin/bash

WORK_DIR=`pwd`

#export VO_CMS_SW_DIR=/uscmst1/prod/sw/cms
#source $VO_CMS_SW_DIR/cmsset_default.sh
source /cvmfs/cms.cern.ch/cmsset_default.csh
export SCRAM_ARCH=slc5_amd64_gcc462
eval `scramv1 runtime -sh`

hadd pileupReweighting_temp_DATASETNAME.root pileupReweighting_DATASETNAME_*.root
rm pileupReweighting_DATASETNAME_*.root 

hadd btagEfficiency_temp_DATASETNAME.root btagEfficiency_DATASETNAME_*.root
rm btagEfficiency_DATASETNAME_*.root 

root -b -q -l combineHistograms.C\(\"DATASETNAME\"\)

rm pileupReweighting_temp_DATASETNAME.root
rm btagEfficiency_temp_DATASETNAME.root
