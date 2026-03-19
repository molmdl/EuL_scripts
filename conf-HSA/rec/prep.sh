#source ~/scripts/gmx2023.5_plumed_gpu_py3.env 
source ~/data/scripts/gmx2023.5_plumed2.9.3_gpu_shared.env

#gmx pdb2gmx -f 2bxo_processed.pdb -ignh -o 2bxo.gro
#gmx editconf -f 2bxo.gro -o box.gro -d 1 -c -bt dodecahedron #-bt cubic
#gmx solvate -cp box.gro -cs spc216 -p topol.top -o solv.gro
gmx grompp -f em.mdp -c solv.gro -p topol.top -o ion.tpr -maxwarn 1
echo SOL | gmx genion -s ion.tpr  -p topol.top  -neutral -conc 0.15 -nname CL -pname NA -o ion.gro
