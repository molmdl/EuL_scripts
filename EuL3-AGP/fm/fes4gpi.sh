#!/bin/bash
sort -k1,1g -k2,2g $1 | gawk '
/^[[:blank:]]*#/ {next}
NF < 3 {next}
$1 != prev {printf "\n"; prev=$1}
{print}
' > output.dat
