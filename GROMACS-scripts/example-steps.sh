gmx pdb2gmx -f ticf.pdb -o ticf.gro -water tip3p -ff amber14sb

# Organize Chain Structure files .itp
mkdir itp_chains
mv *.itp ./itp_chains/
# Update the direction in the topology file
sed -i '' '/^; Include chain topologies$/,/^; Include water topology$/ {
  s#^\(\#include "\)\(.*\.itp\)"#\1itp_chains/\2#
}' topol.top

gmx editconf -f ticf.gro -o ticf-boxd.gro -c -box 16.0 16.0 72.0 -angles 90.0 90.0 90.0

gmx solvate -cp ticf-boxd.gro -cs scp216.gro -o ticf-solvtd.gro -p topol.top

gmx grompp -f ions.mdp -c ticf-solvtd.gro -p topol.top -o ions.tpr

gmx genion -s ions.tpr -o ticf-solvtd-ions.gro -p topol.top -pname NA -nname CL -neutral

# We pass the next steps up to ACENET
# Minimization, NVT Heating Phase, NVT Equilibration Phase, NPT Phase, MD Production
