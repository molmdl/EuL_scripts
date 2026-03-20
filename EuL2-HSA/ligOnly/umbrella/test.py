# find_windows.py - search a distance time series and find the 
# closest value to a supplied target distance, then print the
# time value corresponding to that distance

import numpy
import os 

infile = open('COLVAR', 'r')
lines = infile.readlines()

count = 0.5
for d in numpy.arange(0.2,2.2,0.2):
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
    print (count)
    count += 1

exit()
