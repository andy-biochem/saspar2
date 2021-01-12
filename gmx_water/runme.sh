gmx pdb2gmx -f water_pep.pdb -o water.gro -ignh
gmx editconf -f water.gro -o water_box.gro -c -d 10.0 -bt cubic
gmx solvate -cp water_box.gro -cs tip4p.gro  -o water_solv.gro -p topol.top
gmx grompp -f ions.mdp -c water_solv.gro -p topol.top -o ions.tpr -maxwarn 1
gmx genion -s ions.tpr -o water_ions.gro -p topol.top -pname NA -nname CL -neutral
gmx grompp -f minim.mdp -c 6lyz_ions.gro -p topol.top -o em.tpr 
