#!/bin/bash

WORK_DIR=`pwd`

cd ../../../../
tar -czf src.tgz `cat $WORK_DIR/doc/tarballFileList.txt`
mv src.tgz $WORK_DIR/../

cd $WORK_DIR