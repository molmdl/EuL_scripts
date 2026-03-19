#source ~/scripts/gmx2023.5_plumed_gpu_py3.env 
gmx pdb2gmx -f 3APX_E64V.pdb -ignh -o prot.gro
gmx editconf -f prot.gro -o box.gro -d 1 -c -bt dodecahedron #-bt cubic
gmx solvate -cp box.gro -cs spc216 -p topol.top -o solv.gro
gmx grompp -f em.mdp -c solv.gro -p topol.top -o ion.tpr -maxwarn 1
echo SOL | gmx genion -s ion.tpr  -p topol.top  -neutral -conc 0.15 -nname CL -pname NA -o ion.gro
