#!/bin/bash

CWD=`pwd`

#for s in aagp1 aagp2 aagp1_V92E aagp2_E92V ; do
for s in aagp1 ; do
#for s in aagp2 ; do
#for s in aagp1_V92E ; do
#for s in aagp2_E92V ; do
	cd ${s}_EL3p
	#for i in 0 1 2 3 ; do 
	for i in 1 2 3 ; do 
	#	cd mgpu_t$i ; bash ../../trjconv.sh $i
		cd mgpu_t${i}_1 ; bash ../../trjconv.sh $i
		cd .. 
	done
	#cd mgpu_t0_1
	#bash ../../trjconv.sh 0
	#
	#cd ../ana_trj_1500ns/
	cd ./ana_trj_1500ns/
#	cp ../../${s}_reweight.dat plumed_reweight.dat
	#mkdir mgpu_t0 mgpu_t1 mgpu_t2 mgpu_t3 mgpu_t0_1
	#for i in mgpu_t0 mgpu_t1 mgpu_t2 mgpu_t3 mgpu_t0_1 ; do 
	for i in mgpu_t1_1 mgpu_t2_1 mgpu_t3_1 ; do 
#	for i in mgpu_t2_1 mgpu_t3_1 ; do 
		mkdir $i
		cd $i 
		cp ../../${i}/HILLS . 
		ln -s ../../${i}/fm_1500ns_?.xtc . 
		#ln -s ../../${i}/v1.xtc .
		sbatch ../plumed_reweight.sh
		cd ..
	done
	cd $CWD
done

#for s in aagp1 aagp2 aagp1_V92E aagp2_E92V ; do
#	rm ${s}_EL3p/mgpu_*/fm_1500ns_?.xtc
#	rm ${s}_EL3p/ana_trj_1500ns/mgpu_*/fm_1500ns_?.xtc
#done

#for s in aagp1 aagp2 aagp1_V92E aagp2_E92V ; do
#for s in aagp2 ; do
#	cd ${s}_EL3p
#	cd ./ana_trj_1500ns/
#	cp ../../${s}_merged_reweight.dat plumed_reweight_merged.dat
#	sbatch ${CWD}/plumed_reweight_merged.sh
#	cd $CWD
#done


