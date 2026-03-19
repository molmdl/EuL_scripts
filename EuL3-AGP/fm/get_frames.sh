#!/bin/bash
#SBATCH -J getframes
##SBATCH -p workq
#SBATCH -p cpu
#SBATCH -n 1
#SBATCH -c 1

CWD="/share/home/nglokwan/dparker/aagp/fm"
#python ${CWD}/get_min_frames.py --colvar trj.COLVAR.m --cv_cols fp.lp fd.ld
#python ${CWD}/get_sel_frames.py --colvar trj.COLVAR.m --cv_cols fp.lp fd.ld --filter_cvs dm0.min dm1.min dm2.min

python ${CWD}/get_sel_frames.py --colvar trj.COLVAR.o --cv_cols d1 o1 --filter_cvs dm0.min dm1.min dm2.min

#python ${CWD}/get_sel_frames.py --colvar trj.COLVAR.m --cv_cols fp.lp fd.ld --filter_cvs dm0.min dm2.min
