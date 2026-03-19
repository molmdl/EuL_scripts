Write a python script to analyze the geometry of metal coordination in MD trajectory. Use the currently loaded python environment only. Refer to the example files in me_sssL_sap and phe_sssL_sap folders and run the test analysis on this data.

Given that:
1. The complex has two major geometries around metal center, square antiprism (SAP) and twisted square antiprism (TSAP).
2. The metal center has 9-coordination.
3. The molecule that complex with metal has a 12-N-4 ring with four N-C-C-N torsions, and a 12-N-4 ring to chromophore N-C-C-N torsion. In an SAP geometry, the torsions of 12-N-4 ring have opposite sign of ring-to-chromophore torsion. In a TSAP geometry, the signs of the torsions are the same.

The analysis script should do the following:

I. trajectory loading and processing
1. use MDanalysis to load the topology and trajectory. note that the trajectory may have stripped solvent and hence differnt number of atoms than topology
2. handle pbc and align on the heavy atoms of the metal-ligand complex

II. Geometry analysis A based on alignment to ideal geometry
1. locate the metal center and the 9 coordinating atoms
2. align the ideal gemetry of SAP and TSAP onto the metal 
3. calculate the RMSD to ideal SAP and TSAP geometry 
4. calculate the time-series data of the two RMSD as well as the two histogram of RMSD. output the csv of time-series and histogram, SAP (blue line) and TSAP (purple line) in the same png plot.
5. for each of the peak of the histogram, select a representive frame. for the selected frames, place dummy atoms of the aligned ideal gemontry on the metal center. output the pdb of the frame with the dummy atom

III. Geomtry analysis based on torsions
1. calculate the torsion of the four 12-N-4 N-C-C-N torsions, and the N-C-C-N torsion from 12-N-4 ring to chromophore. Torsional angles should be in degree, from -180 to 180.
2. classify frames of opposite sign as SAP, same sign as TSAP
3. as in II, output the time-series and histogram csv and png plots (for all torsional angles and torsion-based classification)
