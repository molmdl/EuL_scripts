source ~/scripts/gmx2023.5_plumed_gpu_py3.env 
gmx_mpi pdb2gmx -f 3APX_E92V.pdb -ignh -o prot.gro
gmx_mpi editconf -f prot.gro -o box.gro -d 1 -c -bt dodecahedron #-bt cubic
gmx_mpi solvate -cp box.gro -cs spc216 -p topol.top -o solv.gro
gmx_mpi grompp -f em.mdp -c solv.gro -p topol.top -o ion.tpr -maxwarn 1
echo SOL | gmx_mpi genion -s ion.tpr  -p topol.top  -neutral -conc 0.15 -nname CL -pname NA -o ion.gro
