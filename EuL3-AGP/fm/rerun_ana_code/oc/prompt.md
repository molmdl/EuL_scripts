Write a python script to process data from multiple replicas of funnel metadynamics, and create publication-quality plots of the data.
Use only the python libraries in the current conda environment, if prefered also can use the gmx/plumed in the current env. The script should do the following:
1. For each xvg files, for each frame, sum the LJ-SR term and Coul-SR term, calling this the IE value.
  a. from the bias value, get the weight of the frame as in the reweight_metad option in plumed, or just print the weight to a plumed data file  with a plumed run (https://www.plumed.org/doc-v2.9/user-doc/html/_r_e_w_e_i_g_h_t__m_e_t_a_d.html)
  b. If necessary, normalize the weights
  c. note the IE, d1, and w value of each frame
  d. assign these IE values along N bins of the CV of d1
  e. calculate the weighted average of IE of each bin, and standard deviation as error (+/-)
2. Generate a heatmap to plot res-ie.png:
  a. y-axis is residue ID, major tick label every 20 residues, minor ticks without labels every residue, axis label "Residue ID"
  b. x-axis as d1 values range 0-4, major tick label every 0.5, minor ticks every 0.1, axis label "D1 (nm)"
  c. heatmap colored with the weighted average of IE with blue-white-red color scale, where blue is negative value, white is zero, red is positive values
  d. colorbar with colorbar label "Average IE (kJ/mol)"
3. If the user point to multiple replica directory, after calculating for one replica, take average of each replica and provide the plot with the 'merged' prefix.

sample plumed input in the current directory for your reference
* `plumed_reweight_o1.dat`: a working plumed input that plot the reweighted FES of selected CV from trj.COLVAR.o
* plumed command for the run: `plumed driver --plumed ./plumed_reweight_o1.dat --kt 2.494339 --noatoms`

sample files to be analyzed in ./rerun/ of a single replica, for your reference:
* `rerun_*.xvg`: Interaction energy (LJ-SR an Coul-SR) of the residue with ligand from `gmx energy` runs after `gmx mdrun -rerun`
* `trj.COLVAR.o`: plumed data file with weights and other CVs
* `plumed_reweight_o1.dat`: a working plumed input that plot the reweighted FES of selected CV from trj.COLVAR.o

Ask me for any uncertainties or decisions that multiple options available, in that case also provide me your suggestions.
sample files to be analyzed in ./t1_1/ and ./t1/ of two replicas, for your reference:
* `rerun_*.xvg`: Interaction energy (LJ-SR an Coul-SR) of the residue with ligand from `gmx energy` runs after `gmx mdrun -rerun`
* `trj.COLVAR.o`: plumed data file with weights and other CVs

Ask me for any uncertainties or decisions that multiple options available, in that case also provide me your suggestions.
