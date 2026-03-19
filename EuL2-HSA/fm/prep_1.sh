#source /share/scripts/gromacs-2019.env 
#for i in 0 1 2 ; do gmx_mpi grompp -f local_fm_full_100ns.mdp -c pr_${i}.gro -p sys.top -n index.ndx -o funnel6b_${i}/fm_100ns_${i}.tpr -maxwarn 1 ; done
#for i in 0 1 2 ; do cd funnel6b_$i ; bash ~/test/gen_pbs_plumed.sh -q gpu2 -n fm6b_lig2_${i} -gpu $i -tpr fm_100ns_${i} ; cd .. ; done
#for i in 0 1 2 ; do cp plumed.dat protein.pdb funnel6b_${i} ; done
#for i in 0 1 2 ; do mv funnel6b_${i} funnel7_${i} ; done
#for i in 0 1 2 ; do cd funnel7_${i} ; bash ../clean_err.sh ; cd .. ; done
#for i in 0 1 2 ; do cp plumed.dat.7 funnel7_${i}/plumed.dat ; done
#for i in 0 1 2 ; do cd funnel7_${i} ; mv fm6b_lig2_${i}.pbs fm7_lig2_${i}.pbs ; sed -i s/fm6b/fm7/ fm7_lig2_${i}.pbs ; cd .. ; done
#for i in 0 1 2 ; do cd funnel7_${i} ; qsub fm7_lig2_${i}.pbs ; cd .. ; done
for i in 0 1 2 3 ; do mkdir fm7b2cv_${i} ; cp plumed.dat.7b.2cv fm7b2cv_${i}/plumed.dat ; cp protein.pdb fm7b2cv_${i} ; done
#for i in 0 1 2 3 ; do gmx_mpi grompp -f local_fm_full_100ns.mdp -c pr_${i}.gro -p sys.top -n index.ndx -o fm7b2cv_${i}/fm_100ns_${i}.tpr -maxwarn 2 ; done
for i in 0 1 2 3 ; do gmx_mpi grompp -f local_fm_full_100ns.mdp -c pr_${i}.gro -p sys.top -n index.ndx -o fm7b2cv_${i}/fm_100ns_${i}.tpr ; done
