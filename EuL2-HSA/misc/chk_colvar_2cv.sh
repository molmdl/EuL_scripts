date | tee -a last_colvar_2cv.log
for i in `cat list_2cv.txt`; do
	echo $i | tee -a last_colvar_2cv.log
	tail -n1 ${i}/COLVAR | tee -a last_colvar_2cv.log
done
