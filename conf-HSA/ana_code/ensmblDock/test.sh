python dock2com_1.py -i dzp.itp -s hsa*-dzp.sdf -t dzp.mol2 -r '#topol.top.1#' --ff-path ../../amber19SB_OL21_OL3_lipid17.ff/forcefield.itp --water-itp ../../amber19SB_OL21_OL3_lipid17.ff/opc3.itp --ions-itp ../../amber19SB_OL21_OL3_lipid17.ff/ions.itp --lig-gro best.gro --com-gro com.gro --rec-itp rec.itp --sys-top sys.top --metric CNNscore

gmx grompp -p sys.top -c com.gro -f em.mdp -maxwarn 2 -o test.tpr
