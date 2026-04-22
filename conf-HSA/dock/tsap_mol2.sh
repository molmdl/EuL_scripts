#!/bin/bash

for i in sssD sssL rrrD rrrL ; do
       mkdir phe_${i}_tsap
       python ../scripts/dock/gro_itp_to_mol2.py --gro ../solv_md/phe_${i}_tsap/phe_${i}_tsap.gro --itp ../solv_md/phe_${i}_tsap/phe_${i}_tsap.itp --out phe_${i}_tsap/phe_${i}_tsap.mol2 
done
