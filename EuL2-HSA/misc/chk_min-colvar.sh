 for i in `cat list.txt`; do echo $i ; grep -v '#' ${i}/COLVAR | sort -gk 2 | head -n1 ; done
