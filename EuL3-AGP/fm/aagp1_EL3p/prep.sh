source ~/data/scripts/gmx2023.5_plumed2.9.3_gpu_shared.env
gmx_mpi solvate -cp sys0.gro -cs spc216 -p sys.top -o solv.gro
gmx_mpi grompp -f em.mdp -c solv.gro -p sys.top -o ion.tpr -maxwarn 1
echo SOL | gmx_mpi genion -s ion.tpr  -p sys.top  -neutral -conc 0.15 -nname CL -pname NA -o ion.gro
