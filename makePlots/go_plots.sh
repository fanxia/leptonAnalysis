#!/bin/bash

if [ $# -ne 1 ]; then
	echo
	echo "Usage: ./go_plots.sh nPhotons"
	echo
	exit 0
fi

#export VO_CMS_SW_DIR=/uscmst1/prod/sw/cms
#source $VO_CMS_SW_DIR/cmsset_default.sh
source /cvmfs/cms.cern.ch/cmsset_default.csh
export SCRAM_ARCH=slc5_amd64_gcc462
eval `scramv1 runtime -sh`

NUM_PHOTONS_REQUIRED=$1
ELE_FILE_TO_RUN=/eos/uscms/store/user/bfrancis/inputs/SingleElectron.root
MUON_FILE_TO_RUN=/eos/uscms/store/user/bfrancis/inputs/SingleMu.root

rm *.root

cat makePlots_template.C | sed s:ELE_FILE_TO_RUN:$ELE_FILE_TO_RUN: | sed s:MUON_FILE_TO_RUN:$MUON_FILE_TO_RUN: | sed s:NUM_PHOTONS_REQUIRED:$NUM_PHOTONS_REQUIRED: > makePlots.C
root -b -q -l makePlots.C 2>&1 | sed '/does not exist/d' | sed '/has been created/d'
rm makePlots.C

file=`date +"%b%d"`

for x in ele_bjj muon_bjj
do
    
    cp template_errorTable errorTable_$x.tex
    while read line
    do
	code=`echo $line | cut -d : -f 1`
	value=`echo $line | cut -d : -f 2`
	
	sed -i "s/${code}/${value}/g" errorTable_$x.tex
    done < errorTable_$x.temp
    rm errorTable_$x.temp

    latex errorTable_$x.tex
    dvips -Ppdf -t landscape errorTable_$x.dvi
    ps2pdf errorTable_$x.ps
    
    rm errorTable_$x.log errorTable_$x.dvi errorTable_$x.aux errorTable_$x.ps
    
done

mkdir $file
mv *.pdf *.txt *.tex $file
cp *.root $file/
cd $file

tar -czf $file.tgz *
scp $file.tgz $HEP

mv $file.tgz ..
cd ..
rm -r $file
