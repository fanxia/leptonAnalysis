#!/bin/bash

if [ $# -ne 1 ]; then
    echo
    echo "Need a target!"
    echo
    exit 0
fi

FILE=$1
NAME=${FILE%.tex}

latex $NAME.tex
dvips -Ppdf -t landscape $NAME.dvi
ps2pdf $NAME.ps
pdfcrop $NAME.pdf
mv $NAME-crop.pdf $NAME.pdf

#convert -density 300x300 stop-bino_diagram-crop.pdf -transparent white durp.png
#dvipng -D 1000 -bg Transparent -pp 1 $NAME.dvi

rm $NAME.dvi
rm $NAME.ps
rm $NAME.aux
rm $NAME.log

echo
echo $NAME.pdf created!
echo
