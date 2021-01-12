 gmx mdrun -v -deffnm em
 gmx grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr -r em.gro
 gmx mdrun -v -deffnm nvt
 gmx grompp -f npt.mdp -c nvt.gro -p topol.top -o npt.tpr -r em.gro -t nvt.cpt
 gmx mdrun -v -deffnm npt
 gmx grompp -f moldyn.mdp -c npt.gro -p topol.top -o moldyn.tpr -t npt.cpt
 gmx mdrun -v -deffnm moldyn
