#!/bin/bash
#SBATCH -J orca
#SBATCH -w node2
##SBATCH -p cpu
#SBATCH -p workq
##SBATCH -n 36
#SBATCH -n 16
#SBATCH -c 1
##SBATCH --mem=120GB
#SBATCH --mem-per-cpu=8000
##SBATCH --exclusive


source /data/nglokwan/scripts/openmpi_4.1.8.env

#Set actual paths of ORCA and orca_2mkl utility here
ORCA="/data//nglokwan/orca/6.1.0/orca"
orca_2mkl="/data/nglokwan/orca/6.1.0/orca_2mkl"

$ORCA ${1}.inp > ${1}.out

if grep -Fq "ORCA TERMINATED NORMALLY" ${1}.out
then
	echo Done!
else
	echo The task has failed! Please check content of ${1}.out to find reason
	echo The script is terminated
	exit 1
fi

### Convert to .molden file
echo Running orca_2mkl...
$orca_2mkl $1 -molden > /dev/null
cat Nval.txt ${1}.molden.input > ${1}.molden

rm -f ${1}.molden.input Nval.txt

