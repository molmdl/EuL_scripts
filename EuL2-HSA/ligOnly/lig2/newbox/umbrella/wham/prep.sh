for i in {0..20} ; do grep -v '#' ../window${i}/COLVAR | awk '{print $1, $3}' > window${i}.dat ; done
