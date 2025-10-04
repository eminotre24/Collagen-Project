#!/bin/bash
#SBATCH --job-name=addcg-0110
#SBATCH --output=/home/aenovt/out_err/addcg-0110/slurm-%j.out
#SBATCH --error=/home/aenovt/out_err/addcg-0110/slurm-%j.err
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=a00836485@tec.mx

module purge
module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 cuda/12.2 gromacs/2024.4
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

# NON MD run
cd /home/aenovt/scratch/emlndscp-2309/twist0.005-files/
srun gmx grompp -f minimize2.mdp -c em.gro -p topol.top -o em2.tpr
srun gmx mdrun -deffnm em2 -ntomp $OMP_NUM_THREADS -ntmpi $SLURM_NTASKS -nb gpu

cd /home/aenovt/scratch/emlndscp-2309/twist0.01-files/
srun gmx grompp -f minimize2.mdp -c em.gro -p topol.top -o em2.tpr
srun gmx mdrun -deffnm em2 -ntomp $OMP_NUM_THREADS -ntmpi $SLURM_NTASKS -nb gpu

cd /home/aenovt/scratch/emlndscp-2309/twist0.015-files/
srun gmx grompp -f minimize2.mdp -c em.gro -p topol.top -o em2.tpr
srun gmx mdrun -deffnm em2 -ntomp $OMP_NUM_THREADS -ntmpi $SLURM_NTASKS -nb gpu

cd /home/aenovt/scratch/emlndscp-2309/twist0.025-files/
srun gmx grompp -f minimize2.mdp -c em.gro -p topol.top -o em2.tpr
srun gmx mdrun -deffnm em2 -ntomp $OMP_NUM_THREADS -ntmpi $SLURM_NTASKS -nb gpu

cd /home/aenovt/scratch/emlndscp-2309/twist0.03-files/
srun gmx grompp -f minimize2.mdp -c em.gro -p topol.top -o em2.tpr
srun gmx mdrun -deffnm em2 -ntomp $OMP_NUM_THREADS -ntmpi $SLURM_NTASKS -nb gpu

cd /home/aenovt/scratch/emlndscp-2309/twist0.04-files/
srun gmx grompp -f minimize2.mdp -c em.gro -p topol.top -o em2.tpr
srun gmx mdrun -deffnm em2 -ntomp $OMP_NUM_THREADS -ntmpi $SLURM_NTASKS -nb gpu

cd /home/aenovt/scratch/emlndscp-2309/twist0.05-files/
srun gmx grompp -f minimize2.mdp -c em.gro -p topol.top -o em2.tpr
srun gmx mdrun -deffnm em2 -ntomp $OMP_NUM_THREADS -ntmpi $SLURM_NTASKS -nb gpu

cd /home/aenovt/scratch/emlndscp-2309/twist0.055-files/
srun gmx grompp -f minimize2.mdp -c em.gro -p topol.top -o em2.tpr
srun gmx mdrun -deffnm em2 -ntomp $OMP_NUM_THREADS -ntmpi $SLURM_NTASKS -nb gpu

cd /home/aenovt/scratch/emlndscp-2309/twist0.06-files/
srun gmx grompp -f minimize2.mdp -c em.gro -p topol.top -o em2.tpr
srun gmx mdrun -deffnm em2 -ntomp $OMP_NUM_THREADS -ntmpi $SLURM_NTASKS -nb gpu

cd /home/aenovt/scratch/emlndscp-2309/twist0.065-files/
srun gmx grompp -f minimize2.mdp -c em.gro -p topol.top -o em2.tpr
srun gmx mdrun -deffnm em2 -ntomp $OMP_NUM_THREADS -ntmpi $SLURM_NTASKS -nb gpu

cd /home/aenovt/scratch/emlndscp-2309/twist0.07-files/
srun gmx grompp -f minimize2.mdp -c em.gro -p topol.top -o em2.tpr
srun gmx mdrun -deffnm em2 -ntomp $OMP_NUM_THREADS -ntmpi $SLURM_NTASKS -nb gpu

