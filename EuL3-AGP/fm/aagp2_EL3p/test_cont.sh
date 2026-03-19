#source ~/data/scripts/gmx2023.5_plumed2.9.3_gpu_shared.env
#gmx trjconv -s ../em.tpr -f fm_1200ns_0.cpt -fit progressive -o fm_1200ns_0.fit.cpt -n ../index.ndx 
gmx grompp -f cont.mdp -p sys.top -c fm_1200ns_0.fit.gro -t fm_1200ns_0.cpt -o fm_1500ns_0.tpr -v -n ../index.ndx -maxwarn 1
