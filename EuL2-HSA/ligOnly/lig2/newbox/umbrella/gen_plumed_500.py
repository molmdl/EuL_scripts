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
    cmd = "sed s/D2VAL/%.2f/ plumed_template_500.dat > window%d/plumed.dat" % (d, count)
    os.system(cmd)
    count += 1

exit()
