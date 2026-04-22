## Aim
* consider the test files in rec_align. 
* create a standalone python script to align any of the specified pdb and gro structures onto a specified reference pdb. in this example directory, hsa*gro and hsa*pdb are specified structures to be aligned, 2bxf_A.pdb is the reference.
* note that the pdb are directly converted from gro via `gmx editconf`, thus they are identical
* expected outputs from manual alignment of pdbs in e.g. UCSF chimera are provided in rec_align/expected. 
* note that in some cases, target and reference may not have the same number of residues. Resarch for methods th consider such mismatch. For example, UCSF chimera or USalign/TMAlign can do it.

## Workflow
* spawn research agent, planning agent, execute agent, and verifier agent.

## Behavior
* remember do not hardcode paths and file names and any variables, we need a generalized tool
* define all variables at the top of the file
* use only the library in the current environment
* replace for-loops with vectorization whenever possible to speed up
* do not use rm command
* do not alter any existing files
