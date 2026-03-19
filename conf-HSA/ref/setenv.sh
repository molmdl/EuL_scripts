#!/bin/bash
export PATH="/opt/openmpi/3.1.3/bin:$PATH"
export LD_LIBRARY_PATH="/opt/openmpi/3.1.3/lib:$LD_LIBRARY_PATH"
export PATH=$ORCAPATH:$PATH

export OMP_STACKSIZE=200M
ulimit -s unlimited

export Multiwfnpath="/opt/sob/Multiwfn_3.8_dev_bin_Linux"
export PATH="$PATH:/opt/sob/Multiwfn_3.8_dev_bin_Linux"
