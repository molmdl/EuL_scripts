#!/bin/bash
date | tee -a chk_sel.log
for i in `cat list.txt` ; do
	echo `grep Oct ${i} | tail -n1 | awk '{print $6, $7, $8, $9, $10, $2, 0.002*$4}'` $i | tee -a chk_sel.log
done
echo | tee -a chk_sel.log
