source ~/data/scripts/gmx2023.5_plumed2.9.3_gpu_shared.env
gmx editconf -f ${1}.gro -o box.gro -d 2 -bt dodecahedron -c
gmx solvate -p ${1}.top -cp box.gro -cs spc216 -o solv.gro
gmx grompp -f em.mdp -c solv.gro -p ${1}.top -o ion.tpr -maxwarn 2
echo SOL | gmx genion -s ion.tpr -p ${1}.top -neutral -nname CL -pname  NA -conc 0.15 -o ion.gro
gmx grompp -f em.mdp -c ion.gro -p ${1}.top -o em.tpr
#gmx make_ndx -f em.tpr 
#gmx make_ndx -f ${1}.gro 
#gmx genrestr -f ${1}.gro 
gmx mdrun -deffnm em -ntmpi 1 -ntomp 8
gmx grompp -f pr.mdp -c em.gro -p ${1}.top  -r em.gro -o pr.tpr
#for i in 0 1 2 3 ; do 
#	gmx grompp -f pr0.mdp -c pr.gro -p ${1}.top  -o pr_${i}.tpr
#done
#for i in 0 1 2 3 ; do gmx grompp -f md.mdp -c pr_${i}.gro -p ${1}.top   -o prod_${i}.tpr ; done
