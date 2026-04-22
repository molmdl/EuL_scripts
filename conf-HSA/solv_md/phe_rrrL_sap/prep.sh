source ~/data/scripts/gmx2023.5_plumed2.9.3_gpu_shared.env
gmx editconf -f lig.gro -o box.gro -d 2 -bt dodecahedron -c
gmx solvate -p lig.top -cp box.gro -cs spc216 -o solv.gro
gmx grompp -f em.mdp -c solv.gro -p lig.top -o ion.tpr -maxwarn 2
gmx genion -s ion.tpr -p lig.top -neutral -nname CL -pname  NA -conc 0.15 -o ion.gro
gmx grompp -f em.mdp -c ion.gro -p lig.top -o em.tpr
#gmx make_ndx -f em.tpr 
#gmx make_ndx -f lig.gro 
#gmx genrestr -f lig.gro 
gmx mdrun -deffnm em -ntmpi 1 -ntomp 8
gmx grompp -f pr.mdp -c em.gro -p lig.top  -r em.gro -o pr.tpr
#for i in 0 1 2 3 ; do 
#	gmx grompp -f pr0.mdp -c pr.gro -p lig.top  -o pr_${i}.tpr
#done
#for i in 0 1 2 3 ; do gmx grompp -f md.mdp -c pr_${i}.gro -p lig.top   -o prod_${i}.tpr ; done
