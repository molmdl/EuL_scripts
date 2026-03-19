date | tee -a last_colvar.log
for i in `cat list.txt`; do
	echo $i | tee -a last_colvar.log
	tail -n1 ${i}/COLVAR | tee -a last_colvar.log
done
