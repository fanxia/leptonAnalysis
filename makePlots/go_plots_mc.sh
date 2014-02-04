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
cat makePlots_mc_template.C | sed s:FILE_TO_RUN:$FILE_TO_RUN: > makePlots.C
root -b -q -l makePlots.C | sed '/.root does not exist/d'
rm makePlots.C

file=`date +"%b%d"`

rm signal_contamination_stop.root
rm contamination_stop.root
rm plots_*.root

#for x in nojet j b jj bj bb jjj bjj bbj 4j 4j1b 4j2b 5j 5j1b 5j2b 6j 6j1b 6j2b l lj lb ljj lbj lbb ljjj lbjj lbbj ljjjj lbjjj lbbjj ll llj llb lljj llbj llbb
for x in nojet j b bj muJets eleJets hadronic
do

	cp template_errorTable_mc.tex errorTable_$x.tex

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
