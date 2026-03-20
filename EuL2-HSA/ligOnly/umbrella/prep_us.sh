for l in lig2m lig2 ; do
	cd $l
	#for i in `seq 0 20` ; do gmx_mpi grompp -f ../em0.mdp -p sys.top -r window${i}.gro -c window${i}.gro -o window${i}/em0.tpr ; done
	#for i in `seq 0 20` ; do gmx_mpi grompp -f ../em.mdp -p sys.top -r window${i}/em0.gro -c window${i}/em0.gro -o window${i}/em.tpr ; done
	#for i in `seq 0 20` ; do gmx_mpi grompp -f ../pr_pos.mdp -p sys.top -r window${i}/em.gro -c window${i}/em.gro -o window${i}/pr_pos.tpr ; done
	#for i in `seq 0 20` ; do gmx_mpi grompp -f ../pr0.mdp -p sys.top -c window${i}/pr_pos.gro -o window${i}/pr_0.tpr ; done
	for i in `seq 0 20` ; do gmx_mpi grompp -f ../umbrella_30ns.mdp -p sys.top -c window${i}/pr_0.gro -o window${i}/umbrella${i}.tpr ; done
	cd ..
done

