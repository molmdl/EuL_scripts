for i in rrrL rrrD sssD sssL ; do echo $i sap `echo 13 0 | gmx energy -f phe_${i}_sap/prod_0.edr 2>/dev/null | grep Potential | awk '{print $2}'` ; echo $i tsap `echo 13 0 | gmx energy -f phe_${i}_tsap/prod_0.edr 2>/dev/null | grep Potential | awk '{print $2}'` ; done
for i in rrrL rrrD sssD sssL ; do echo $i sap `echo 13 0 | gmx energy -f me_${i}_sap/prod_0.edr 2>/dev/null | grep Potential | awk '{print $2}'` ; echo $i tsap `echo 13 0 | gmx energy -f me_${i}_tsap/prod_0.edr 2>/dev/null | grep Potential | awk '{print $2}'` ; done
rm energy.xvg \#energy.xvg*
