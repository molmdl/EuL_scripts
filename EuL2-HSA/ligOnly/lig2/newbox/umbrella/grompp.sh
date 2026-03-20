##for i in {0..20} ; do gmx_mpi grompp -f em0.mdp -c window${i}/ion.gro -p sys.top -o window${i}/em0.tpr -r window${i}/ion.gro -maxwarn 1 ; done
#for i in {0..20} ; do gmx_mpi grompp -f em0.mdp -c window${i}.gro -p sys.top -o window${i}/em0.tpr -r window${i}.gro -maxwarn 2 ; done
#for i in {0..20} ; do gmx_mpi grompp -f pr_pos.mdp -c window${i}/em0.gro -p sys.top -o window${i}/pr_pos.tpr -r window${i}.gro ; done
#for i in {0..20} ; do gmx_mpi grompp -f pr0.mdp -c window${i}/pr_pos.gro -p sys.top -o window${i}/pr0.tpr ; done
for i in {0..20} ; do gmx_mpi grompp -f umbrella_10ns.mdp -c window${i}/pr0.gro -p sys.top -o window${i}/window${i}.tpr ; done
