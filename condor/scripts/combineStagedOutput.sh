#!/bin/bash

WORK_DIR=`pwd`

export VO_CMS_SW_DIR=/uscmst1/prod/sw/cms
source $VO_CMS_SW_DIR/cmsset_default.sh
export SCRAM_ARCH=slc5_amd64_gcc462
eval `scramv1 runtime -sh`

hadd signal_contamination_DATASETNAME.root signal_contamination_DATASETNAME_*.root
rm signal_contamination_DATASETNAME_*.root
