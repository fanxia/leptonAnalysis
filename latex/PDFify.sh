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

rm $NAME.dvi
rm $NAME.ps
rm $NAME.aux
rm $NAME.log

echo
echo $NAME.pdf created!
echo
