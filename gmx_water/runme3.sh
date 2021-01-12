 gmx trjconv -s moldyn.tpr -f moldyn.trr -o gmx_water.pdb
 cat gmx_water.pdb |grep -v HW1 |grep -v HW2 |grep -v MW |gzip  >gmx_water.pdb.gz
