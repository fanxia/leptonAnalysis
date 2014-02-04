#!/bin/sh

JOB_NUMBER=$1
let FILE_NUMBER=$JOB_NUMBER+1

WORK_DIR=`pwd`

export VO_CMS_SW_DIR=/uscmst1/prod/sw/cms
source $VO_CMS_SW_DIR/cmsset_default.sh
export SCRAM_ARCH=slc5_amd64_gcc462
scramv1 project CMSSW CMSSW_5_3_8_patch3

cd CMSSW_5_3_8_patch3
mv $WORK_DIR/src.tgz .
rm -r src/
tar -xzf src.tgz
cd src/
eval `scramv1 runtime -sh`
scramv1 b -j 4

cd SusyAnalysis/SusyNtuplizer/leptonAnalysis/
make
cd ../
mv $WORK_DIR/runOverAOD.py .

INFILE=`sed -n "${FILE_NUMBER}p" $WORK_DIR/fileList.txt`
rm $WORK_DIR/fileList.txt
sed -i "s/FILENAME/${INFILE}/g" runOverAOD.py 

cmsRun runOverAOD.py

OUTFILE="tree_"$INFILE  

iscopy=1
numfail=0
while [[ $iscopy -ne 0 && $numfail -lt 10 ]]
do
	cp susyEvents.root /eos/uscms/store/user/lpcpjm/PrivateMC/FastSim/533p3_full/naturalHiggsinoNLSP/SusyNtuple/cms538v1/$OUTFILE
        iscopy=$?
        numfail=$(( $numfail + 1 ))
done

echo
echo "Tried to copy $numfail times, ended with status $iscopy"
echo

cd $WORK_DIR
rm -rf CMSSW_5_3_8_patch3/
