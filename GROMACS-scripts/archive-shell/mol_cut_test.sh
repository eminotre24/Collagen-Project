# Collagen Molecule Cutted Test
gmx pdb2gmx -f colmol_restrim.pdb -o colmol.gro -water tip3p -ff amber14sb
gmx editconf -f colmol.gro -o colmol_boxd.gro -c -d 2.0
gmx solvate -cp colmol_boxd.gro -cs spc216.gro -o colmol_solvtd.gro -p topol.top
gmx grompp -f ions.mdp -c colmol_solvtd.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o colmol_solvtd_neut.gro -p topol.top -pname NA -nname CL -neutral
# Seems as the new collagen molecule was negatively charged.. NA ions where added.
gmx grompp -f minimize.mdp -c colmol_solvtd_neut.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em
gmx energy -f em.edr -o datafiles/potential.xvg
gmx grompp -f nvt_heating.mdp -c em.gro -r em.gro -p topol.top -o nvt_heating.tpr
gmx mdrun -v -deffnm nvt_heating

gmx grompp -f nvt_heating.mdp -c em.gro -r em.gro -p topol.top -o nvt_heating.tpr
gmx mdrun -v -deffnm nvt_heating

# Prepare and Run NPT Phase
gmx grompp -f npt.mdp -c nvt_equilibration.gro -r nvt_equilibration.gro -t nvt_equilibration.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt

gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr
gmx mdrun -deffnm md_0_1


