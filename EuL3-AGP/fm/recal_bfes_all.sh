#!/bin/bash

###bash cal_o1_all.sh 
##for i in aagp1 aagp2 aagp1_V92E aagp2_E92V ; do 
##	#mkdir  ${i}_EL3p/ana_trj_1500ns/bfes_o1_5t
##	cat ${i}_EL3p/ana_trj_1500ns/mgpu_t[0123]_o1/trj.COLVAR.o ${i}_EL3p/ana_trj_1500ns/mgpu_t0_1_o1/trj.COLVAR.o > ${i}_EL3p/ana_trj_1500ns/bfes_o1_5t/COLVAR_merged
##	wc -l ${i}_EL3p/ana_trj_1500ns/bfes_o1_5t/COLVAR_merged 
##done
##
#CWD=`pwd`
#
#for i in aagp1 aagp2 aagp1_V92E aagp2_E92V ; do 
#	cd ${i}_EL3p/ana_trj_1500ns/bfes_o1_5t
#	sbatch << EOF
##!/bin/bash
##SBATCH -J fes_5t_$i
##SBATCH -p workq
##SBATCH -c 1
##SBATCH -n 1
#
#plumed driver --plumed ${CWD}/plumed_reweight_o1_m.dat --kt 2.494339 --noatoms
#
#EOF
#	cd $CWD
#done
#CWD=`pwd`
###for i in aagp1 aagp2 aagp1_V92E aagp2_E92V ; do 
##for i in aagp2 ; do 
##	mkdir  ${i}_EL3p/ana_trj_1500ns/bfes_o1_8t
##	cat ${i}_EL3p/ana_trj_1500ns/mgpu_t[0123]_o1/trj.COLVAR.o ${i}_EL3p/ana_trj_1500ns/mgpu_t[0123]_1_o1/trj.COLVAR.o > ${i}_EL3p/ana_trj_1500ns/bfes_o1_8t/COLVAR_merged
##	wc -l ${i}_EL3p/ana_trj_1500ns/bfes_o1_8t/COLVAR_merged 
##done
#
##for i in aagp1 aagp2 aagp1_V92E aagp2_E92V ; do 
#for i in aagp2 ; do 
#	cd ${i}_EL3p/ana_trj_1500ns/bfes_o1_8t
#	sbatch << EOF
##!/bin/bash
##SBATCH -J fes_8t_$i
##SBATCH -p workq
##SBATCH -c 1
##SBATCH -n 1
#
#plumed driver --plumed ${CWD}/plumed_reweight_o1_m.dat --kt 2.494339 --noatoms
#
#EOF
#	cd $CWD
#done


CWD=`pwd` 
#for i in aagp1 aagp2 aagp1_V92E aagp2_E92V ; do
#	cd ${i}_EL3p/ana_trj_1500ns/bfes_o1_5t
#	rm delta*txt bfes_cal.log
#	bash ${CWD}/ffs_cal.sh
#	cd $CWD
#done
#
for i in aagp1 aagp2 aagp1_V92E aagp2_E92V ; do
#for i in aagp2 aagp1_V92E ; do
#for i in aagp2_E92V ; do
#for i in aagp1_E92V ; do
	cd ${i}_EL3p/ana_trj_1500ns/bfes_o1_8t_d
	rm delta*txt bfes_cal2.log
	bash ${CWD}/ffs_cal.sh
	cd $CWD
done

grep DeltaG */ana_trj_1500ns/bfes_o1_?t_d/bfes_cal2.log

#grep DeltaG */ana_trj_1500ns/bfes_o1_?t_d/bfes_cal1.log
#for i in aagp2 ; do
#	cd ${i}_EL3p/ana_trj_1500ns/bfes_o1_8t_d
##	rm delta*txt bfes_cal.log
#	bash ${CWD}/ffs_cal.sh
#	cd $CWD
#done
#
#grep DeltaG */ana_trj_1500ns/bfes_o1_?t_d/bfes_cal1.log
#
