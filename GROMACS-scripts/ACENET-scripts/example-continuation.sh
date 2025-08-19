#!/bin/bash
#SBATCH --job-name=fullcolmdcont2
#SBATCH --output=/home/aenovt/out_err/fullcolmd-1508/out.out
#SBATCH --error=/home/aenovt/out_err/fullcolmd-1508/err.err
#SBATCH --time=5-00:00:00             # Max Time
#SBATCH --nodes=1					  # Resources based on Benchmark ID=205
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12			  #	OPENMP
#SBATCH --gres=gpu:1				  # GPU
#SBATCH --mem-per-cpu=2G              # RAM
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=a00836485@tec.mx

# Load GROMACS (In parallel for optimization of the system)
module purge
module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 cuda/12.2 gromacs/2024.4          # Load gromacs and other

# Move to output dir to save files there
cd /home/aenovt/scratch/fullcolmd100-files

# Get the threads from the variable cpus per task - should be same number
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

# Continue Production Phase
srun gmx mdrun -s md_1.tpr -cpi md_1.cpt -deffnm md_1 -append -ntomp $OMP_NUM_THREADS -ntmpi $SLURM_NTASKS -nb gpu -pme gpu -update gpu -bonded cpu
