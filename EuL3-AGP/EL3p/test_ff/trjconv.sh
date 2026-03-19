source ~/scripts/gmx2023.5_plumed_gpu_shared.env 
gmx_mpi trjconv -f lig_eq_0.xtc -s lig_eq_0.tpr -pbc mol -center -o pbcmol.xtc
gmx_mpi trjconv  -o v1.xtc -s lig_eq_0.tpr -fit rot+trans -f pbcmol.xtc
rm pbcmol.xtc 
gmx_mpi trjconv  -f v1.xtc -s lig_eq_0.tpr -o v1.pdb -dump 0
