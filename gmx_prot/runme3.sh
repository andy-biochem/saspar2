 gmx trjconv -s moldyn.tpr -f moldyn.trr -o gmx_`cat protein.name`.pdb
 cat gmx_`cat protein.name`.pdb |grep -v HW1 |grep -v HW2 |grep -v MW |gzip  >gmx_`cat protein.name`.pdb.gz
