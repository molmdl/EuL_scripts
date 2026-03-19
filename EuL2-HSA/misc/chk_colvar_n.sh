date | tee -a last_colvar_n.log
for i in `cat list_n.txt`; do
	echo $i | tee -a last_colvar_n.log
	tail -n1 ${i}/COLVAR | tee -a last_colvar_n.log
done
