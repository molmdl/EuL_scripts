#!/bin/bash
CWD=`pwd` ; for i in `cat list.txt`; do cd $i ; gnuplot ${CWD}/chk_funnel2d-time.gpi ; cd $CWD ; done
