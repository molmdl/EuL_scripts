#!/bin/bash
CWD=`pwd` ; for i in `cat list_2cv.txt`; do cd $i ; gnuplot ${CWD}/chk_fm7b2cv.gpi ; cd $CWD ; done
