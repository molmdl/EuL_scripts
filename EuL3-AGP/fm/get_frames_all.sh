#!/bin/bash

CWD=`pwd`

#for s in aagp1_EL3p aagp2_EL3p aagp1_V92E_EL3p aagp2_E92V_EL3p ; do
#for s in aagp1_EL3p aagp2_E92V_EL3p ; do
#for s in aagp2_EL3p ; do
#for s in aagp2_EL3p aagp1_V92E_EL3p ; do
#for s in aagp2_EL3p aagp1_V92E_EL3p ; do
#for s in aagp2_EL3p aagp1_V92E_EL3p aagp2_E92V_EL3p ; do
#for s in aagp2_E92V_EL3p ; do
#for s in aagp1_V92E_EL3p aagp2_E92V_EL3p ; do
#for s in aagp2_EL3p ; do
for s in aagp1_EL3p ; do
	cd $s
#	for i in mgpu_t0 mgpu_t0_1 mgpu_t1 mgpu_t2 mgpu_t3 ; do
#		cd $i 
#		bash ../../trjconv.sh 
#		cd .. 
#	done
	cd ana_trj_1500ns/
#	for i in mgpu_t0 mgpu_t0_1 mgpu_t1 mgpu_t2 mgpu_t3 ; do
	for i in mgpu_t0 mgpu_t0_1 mgpu_t1 mgpu_t2 mgpu_t3  mgpu_t1_1 mgpu_t2_1 mgpu_t3_1 ; do
		mkdir ${i}_o1_d_ext_fes8t
		cd ${i}_o1_d_ext_fes8t
		ln -s ../../${i}/v1.??? . 
		#cp ../bfes_o1_8t/fes_d1o1.dat fes.dat
		cp ../bfes_o1_8t_d/merged_fes.dat fes.dat
		cp ../${i}_o1_d/trj.COLVAR.o .
		sbatch << EOF
#!/bin/bash
#SBATCH -J getframes_${s}_$i
#SBATCH -p workq
##SBATCH -p cpu
#SBATCH -n 1
#SBATCH -c 1

CWD="/share/home/nglokwan/dparker/aagp/fm"
#python ${CWD}/get_sel_frames.py --colvar trj.COLVAR.o --cv_cols d1 o1 --filter_cvs dm0.min dm1.min dm2.min
python ${CWD}/get_sel_frames.py --colvar trj.COLVAR.o --cv_cols d1 o1 --filter_cvs dm0.min dm2.min

EOF
		cd ..
	done
	cd $CWD
done
