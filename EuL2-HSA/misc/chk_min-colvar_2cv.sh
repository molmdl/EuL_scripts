 for i in `cat list_2cv.txt`; do echo $i ; grep -v '#' ${i}/COLVAR | sort -gk 3 | head -n1 ; done
