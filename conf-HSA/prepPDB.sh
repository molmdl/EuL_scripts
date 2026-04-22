#!/bin/bash
pdb2pqr ${1}.pdb ${1}_fixed.pdb
echo -e '1 \n 2' | gmx pdb2gmx -f ${i}_fixed.pdb -ignh -o ${i}.gro
