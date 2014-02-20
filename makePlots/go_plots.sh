#!/bin/bash

if [ $# -ne 1 ]; then
	echo
	echo "Need a root file input"
	echo
	exit 0
fi

export VO_CMS_SW_DIR=/uscmst1/prod/sw/cms
source $VO_CMS_SW_DIR/cmsset_default.sh
export SCRAM_ARCH=slc5_amd64_gcc462
eval `scramv1 runtime -sh`

FILE_TO_RUN=$1
cat makePlots_template.C | sed s:FILE_TO_RUN:$FILE_TO_RUN: > makePlots.C
root -b -q -l makePlots.C | sed '/.root does not exist/d'
rm makePlots.C

file=`date +"%b%d"`

rm signal_contamination_stop.root
rm contamination_stop.root
rm plots_*.root

if [ 1 -eq 0 ]; then
    for x in ele
    do
	
	cp template_errorTable.tex errorTable_$x.tex
	
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
fi

hadd contamination_stop.root signal_*.root
rm signal_*.root
mv contamination.root signal_contamination_stop.root
	
mkdir $file
mv *.pdf *.txt *.tex $file
cp *.root $file/
cd $file

tar -czf $file.tgz *
scp $file.tgz $DUMP

mv $file.tgz ..
cd ..
rm -r $file
