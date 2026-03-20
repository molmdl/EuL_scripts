# find_windows.py - search a distance time series and find the 
# closest value to a supplied target distance, then print the
# time value corresponding to that distance

import numpy
import os 

infile = open('trj.COLVAR.sel', 'r')
lines = infile.readlines()

count = 0
for d in numpy.arange(-1.2,1.21,0.12):
    diff = 100000
    keeptime = 0
    for line in lines:
        if not (line.startswith('#') or (line.startswith('@'))):
            tmp = line.split()
            time = float(tmp[0])
            distance = float(tmp[2])
            tmpdiff = abs(d - distance)
            if (tmpdiff < diff):
                diff = tmpdiff
                keeptime = time
                keepdist = distance
    cmd = "echo 0 | gmx_mpi trjconv -quiet -s ../prod_300ns.tpr -f ../prod_300ns.xtc -dump %f -o window%d.gro" % (keeptime,count)
    os.system(cmd)
    print ("target value: %.2f best value: %.3f at time: %.1f" % (d, keepdist, keeptime))
    count += 1

exit()
