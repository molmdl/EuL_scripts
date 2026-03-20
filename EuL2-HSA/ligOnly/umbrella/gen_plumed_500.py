# find_windows.py - search a distance time series and find the 
# closest value to a supplied target distance, then print the
# time value corresponding to that distance

import numpy
import os 

infile = open('trj.COLVAR', 'r')
lines = infile.readlines()

count = 0
for d in numpy.arange(-1.5,1.6,0.15):
    diff = 100000
    keeptime = 0
    for line in lines:
        if not (line.startswith('#') or (line.startswith('@'))):
            tmp = line.split()
            time = float(tmp[0])
            dih = float(tmp[1])
            tmpdiff = abs(d - dih)
            if (tmpdiff < diff):
                diff = tmpdiff
                keeptime = time
                keepdist = dih
    cmd = "sed s/DIHVAL/%.1f/ ../plumed_template_500.dat > window%d/plumed.dat" % (d, count)
    os.system(cmd)
    count += 1

exit()
