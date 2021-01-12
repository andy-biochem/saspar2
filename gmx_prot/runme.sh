gmx pdb2gmx -f `cat protein.name`.pdb -o `cat protein.name`_gmx.gro -ignh
gmx editconf -f `cat protein.name`_gmx.gro -o `cat protein.name`_box.gro -c -d 2.7 -bt cubic
gmx solvate -cp `cat protein.name`_box.gro -cs tip4p.gro  -o `cat protein.name`_solv.gro -p topol.top
gmx grompp -f ions.mdp -c `cat protein.name`_solv.gro -p topol.top -o ions.tpr -maxwarn 1
gmx genion -s ions.tpr -o `cat protein.name`_ions.gro -p topol.top -pname NA -nname CL -neutral
gmx grompp -f minim.mdp -c `cat protein.name`_ions.gro -p topol.top -o em.tpr 
