work on files in the ensmblDock directory

## Aim
* combine mol2_reorder.py and dock2com.py into a new dock2com_1.py script that do the same job as the test_original.sh which call these two scripts with cli options
* include a function to get the correct receptor gro, since we have an ensemble here. currently I get the receptor id  after running mol2_reorder.py which find the best model, which is manual, find a way to make it auto
* no hardcode if it can be variable, we need a general script

## Hints
* consider test_original.sh that run both of the python scripts to get the output in expected/
   * consider the script mol2_reorder.py that process ligand from multiple docking output sdf 
   * consider the script dock2com.py that combine a docked ligand and receptor into complex file for MD simulation.
* consider the gnina_test.sh being used to generate do all the docking, see which file corresponds to which step in this script
   * you may call gnina to check its help text, its in the path.

## Workflow
* spawn research agent, planning agent, execute agent, and verifier agent.

## Behavior
* ask me for any uncertainties
* remember do not hardcode paths and file names and any variables, we need a generalized tool
* define all variables at the top of the file
* use only the library in the current environment
* replace for-loops with vectorization whenever possible to speed up
* do not use rm command
* do not alter any existing files
