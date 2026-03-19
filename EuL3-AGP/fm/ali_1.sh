echo 1 | gmx_mpi editconf -f aagp10.gro -d 1 -bt cubic -c -o aagp1_box.gro -princ
echo 1 | gmx_mpi editconf -f aagp20_ali.gro -d 1 -bt cubic -c -o aagp2_box.gro -princ
