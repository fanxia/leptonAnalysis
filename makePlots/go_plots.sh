#!/bin/bash

if [ $# -ne 3 ]; then
	echo
	echo "Usage: ./go_plots.sh nPhotons electrons.root muons.root"
	echo
	exit 0
fi

export VO_CMS_SW_DIR=/uscmst1/prod/sw/cms
source $VO_CMS_SW_DIR/cmsset_default.sh
export SCRAM_ARCH=slc5_amd64_gcc462
eval `scramv1 runtime -sh`

NUM_PHOTONS_REQUIRED=$1
ELE_FILE_TO_RUN=$2
MUON_FILE_TO_RUN=$3

cat makePlots_template.C | sed s:ELE_FILE_TO_RUN:$ELE_FILE_TO_RUN: | sed s:MUON_FILE_TO_RUN:$MUON_FILE_TO_RUN: | sed s:NUM_PHOTONS_REQUIRED:$NUM_PHOTONS_REQUIRED: > makePlots.C
root -b -q -l makePlots.C 2>&1 | sed '/does not exist/d' | sed '/has been created/d'
rm makePlots.C

file=`date +"%b%d"`

rm signal_contamination_stop.root
rm contamination_stop.root
rm plots_*.root

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

cp template_datacard_ele datacard_ele_bjj.dat
while read line
do
    code=`echo $line | cut -d : -f 1`
    value=`echo $line | cut -d : -f 2`
    
    sed -i "s/${code}/${value}/g" datacard_ele_bjj.dat
done < datacard_ele_bjj.temp

cp template_datacard_muon datacard_muon_bjj.dat
while read line
do
    code=`echo $line | cut -d : -f 1`
    value=`echo $line | cut -d : -f 2`
    
    sed -i "s/${code}/${value}/g" datacard_muon_bjj.dat
done < datacard_muon_bjj.temp

hadd contamination_stop.root signal_*.root
rm signal_*.root
mv contamination.root signal_contamination_stop.root
	
mkdir $file
mv *.pdf *.txt *.tex $file
cp *.root $file/
cd $file

tar -czf $file.tgz *
scp $file.tgz $HEP

mv $file.tgz ..
cd ..
rm -r $file
