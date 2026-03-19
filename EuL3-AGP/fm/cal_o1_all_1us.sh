#!/bin/bash

CWD=`pwd`

#for s in aagp1_EL3p aagp2_EL3p aagp1_V92E_EL3p aagp2_E92V_EL3p ; do
#for s in aagp1_EL3p aagp2_E92V_EL3p ; do
#for s in aagp2_EL3p ; do
#for s in aagp2_EL3p aagp1_V92E_EL3p ; do
for s in aagp2_EL3p aagp1_V92E_EL3p aagp2_E92V_EL3p ; do
#for s in aagp1_V92E_EL3p aagp2_E92V_EL3p ; do
#for s in aagp1_V92E_EL3p ; do
#for s in aagp2_E92V_EL3p ; do
#for s in aagp2_EL3p ; do
	cd $s
#	for i in mgpu_t0 mgpu_t0_1 mgpu_t1 mgpu_t2 mgpu_t3 ; do
#		cd $i 
#		bash ../../trjconv.sh 
#		cd .. 
#	done
	cd ana_trj_1500ns/
#	for i in mgpu_t0 mgpu_t0_1 mgpu_t1 mgpu_t2 mgpu_t3 ; do
	for i in mgpu_t0 mgpu_t0_1 mgpu_t1 mgpu_t2 mgpu_t3  mgpu_t1_1 mgpu_t2_1 mgpu_t3_1 ; do
#	for i in mgpu_t1_1 mgpu_t2_1 mgpu_t3_1 ; do
		if [[ -e "../${i}/v1.xtc" ]] ; then
			mkdir ${i}_o1_d_800ns
			cd ${i}_o1_d_800ns
			echo calculating $s $i ...
			sbatch << EOF
#!/bin/bash
#SBATCH -J fes_o1_${s}_$i
#SBATCH -p workq
#SBATCH -c 1
#SBATCH -n 1 

##rm o1.dat
##mv o1.dat o1.dat.bak
#mv o1.dat o1.dat.bak.1
#
#ln -s ../../${i}/v1.??? . 
#python ../../../cal_o1_abs.py
##paste -d " " ../${i}/trj.COLVAR.m o1.dat > trj.COLVAR.o
#paste -d " " ../${i}/trj.COLVAR o1.dat > trj.COLVAR.o
#sed -i s/fps.ld/fd.ld/ trj.COLVAR.o
#sed -i s/fps.lp/fp.lp/ trj.COLVAR.o

#head -50002 ../${i}_o1_d/trj.COLVAR.o > trj.COLVAR.o
head -40002 ../${i}_o1_d/trj.COLVAR.o > trj.COLVAR.o

plumed driver --plumed ${CWD}/plumed_reweight_o1.dat --kt 2.494339 --noatoms

EOF
			echo
		else
			echo ${s}/${i} does not exist!
		fi
		cd ..
	done
	cd $CWD
done

#find -name fes_d1o1.dat | wc -l
