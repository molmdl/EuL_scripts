for i in * ; do grep -A6 Energy ${i}/mmpbsa/averaged_energies.txt |grep -v '\-\-\-\-\-\-' > ${i}_avg_mmpbsa.dat ; done
