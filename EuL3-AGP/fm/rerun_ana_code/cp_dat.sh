#!/bin/bash
CWD=`pwd` 
#for i in  aagp1 aagp2 ; do 
#for i in aagp1_V92E aagp2_E92V ; do 
for i in  aagp1 aagp2 aagp1_V92E aagp2_E92V ; do 
	for j in 0 1 2 3 0_1 1_1 2_1 3_1 ; do 
		mkdir -p ${i}/mgpu_t${j}/rerun 
		cp ../${i}_EL3p/ana_trj_1500ns/mgpu_t${j}_o1_d_ext_fes8t/trj.COLVAR.o.xz ${i}/mgpu_t${j}/rerun 
		cp ../${i}_EL3p/mgpu_t${j}/rerun/*xvg.xz ${i}/mgpu_t${j}/rerun 
		cd ${i}/mgpu_t${j}/rerun 
		xz -d trj.COLVAR.o.xz *xvg.xz
		cd $CWD 
	done
	cd $i
	bash ${CWD}/plot.sh
	cd $CWD
	rm ${i}/mgpu_t*/rerun/trj.COLVAR.o ${i}/mgpu_t*/rerun/*.xvg
done
