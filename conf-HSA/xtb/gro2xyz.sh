for i in `cat list.txt` ; do mkdir -p $i ; cp ../solv_md/${i}/box.gro $i ; cd $i ; echo -e "xyz \n input.xyz \n q" | Multiwfn box.gro ; cd .. ; done
