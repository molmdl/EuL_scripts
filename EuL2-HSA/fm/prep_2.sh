source /share/scripts/gromacs-2019.env 
for i in 0 1 2 3 ; do 
	mkdir fm7b_${i}
	cp plumed.dat.7b fm7b_${i}/plumed.dat
	cp protein.pdb fm7b_${i}
	gmx_mpi grompp -f local_fm_full_300ns.mdp -c pr_${i}.gro -p sys.top -n index.ndx -o fm7b_${i}/fm_300ns_${i}.tpr
	cd fm7b_$i
	bash ~/test/gen_pbs_plumed.sh -q gpu2 -n fm7b_lig2_${i} -gpu $i -tpr fm_300ns_${i}
	cd ..
done
#for i in 0 1 2 ; do cd funnel7_${i} ; qsub fm7_lig2_${i}.pbs ; cd .. ; done

