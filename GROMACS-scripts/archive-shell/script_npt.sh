#!/bin/bash
#SBATCH --job-name=npt
#SBATCH --output=npt.out
#SBATCH --error=npt.err
#SBATCH --time=1-00:00:00             # Max Time
#SBATCH --partition=gpubase_bygpu_b3  # Look up using the command: sinfo -o "%P %D %C %m %l" (-o for clear output)
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=10			  #	OPENMP
#SBATCH --gres=gpu:1				  # GPU
#SBATCH --mem=64G                     # RAM
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=a00836485@tec.mx

# Load GROMACS (CPU version)
module purge
module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 cuda/12.2 gromacs/2024.4          # Load gromacs and other

# Move to output dir to save files there
cd /home/aenovt/scratch/minnvt

# Get the threads from the variable cpus per task - should be same number
export OMP_NUM_THREADS=10

# Run the simulation - EM
srun gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
srun gmx mdrun -deffnm npt -ntomp $OMP_NUM_THREADS -ntmpi $SLURM_NTASKS