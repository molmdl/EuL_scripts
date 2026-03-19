#python ../../../fes_ana.py calculate --fes ../fes_fps.dat --xmin 0.2 --xmax 1.5 --ymin 0 --ymax 0.6 --wref 6 --rcyl 0.1 | tee bfes_cal.log
#python ../../../fes_ana.py calculate --fes ./fes_fps.dat --xmin 0.2 --xmax 1.5 --ymin 0 --ymax 0.6 --wref 6 --rcyl 0.1 >> bfes_cal.log 2>> bfes_cal.log
#python ../../../fes_ana.py calculate --fes ./fes_fps.dat --xmin 0.2 --xmax 2 --ymin 0 --ymax 0.8 --wref 6 --rcyl 0.1 >> bfes_cal.log 2>> bfes_cal.log
#python ../../../fes_ana.py calculate --fes ./fes_d1o1.dat --xmin 0.2 --xmax 2 --ymin 0 --ymax 90 --wref 6 --rcyl 0.1 >> bfes_cal.log 2>> bfes_cal.log
#python ../../../fes_ana.py calculate --fes ./fes_d1o1.dat --xmin 0.2 --xmax 3 --ymin -180 --ymax 180 --wref 6 --rcyl 0.1 >> bfes_cal.log 2>> bfes_cal.log
#python ../../../fes_ana.py calculate --fes ./merged_fes.dat --xmin 0.2 --xmax 3 --ymin -180 --ymax 180 --wref 6 --rcyl 0.1 >> bfes_cal.log 2>> bfes_cal.log

#python ../../../fes_ana.py calculate --fes ./merged_fes.dat --xmin 0.2 --xmax 2 --ymin -180 --ymax 180 --wref 6 --rcyl 0.1 >> bfes_cal1.log 2>> bfes_cal1.log
python ../../../bfes_code/fes_ana.py calculate --fes ./merged_fes.dat --xmin 0.2 --xmax 1 --ymin -180 --ymax 180 --wref 4 --rcyl 0.1 >> bfes_cal2.log 2>> bfes_cal2.log
#python ../../../fes_ana.py calculate --fes ./merged_fes.dat --xmin 0.2 --xmax 1 --ymin -180 --ymax 180 --wref 4.1 --rcyl 0.1 >> bfes_cal1.log 2>> bfes_cal1.log

