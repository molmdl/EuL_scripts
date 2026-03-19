source ~/scripts/gmx2023.5_plumed_gpu_py3.env 
gmx_mpi trjconv -f apo_100ns.xtc -s apo_100ns.tpr -pbc mol -center -o pbcmol.xtc
gmx_mpi trjconv -f pbcmol.xtc -s apo_100ns.tpr -fit rot+trans -o apo_100ns_pbc.xtc
rm pbdmol.xtc 
gmx_mpi trjconv -s apo_100ns.tpr -dump 0 -f apo_100ns_pbc.xtc -o apo_100ns_pbc_0.pdb
gmx_mpi rms -f apo_100ns_pbc.xtc -s apo_100ns.tpr -o apo_100ns_bb.rmsd
gmx_mpi rmsf -s apo_100ns.tpr -f apo_100ns_pbc.xtc -oq apo100ns_rmsf_bfac.pdb -o apo100ns_bb_rmsf.xvg -res -dir apo100ns_bb_rmsf.log -oc apo100ns_bb_rmsf.xvg
gmx_mpi rmsf -s apo_100ns.tpr -f apo_100ns_pbc.xtc -oq apo100ns_rmsf_noh_bfac.pdb -o apo100ns_noh_rmsf.xvg -res -dir apo100ns_noh_rmsf.log -oc apo100ns_noh_rmsf.xvg
 
