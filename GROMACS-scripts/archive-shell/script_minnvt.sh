#!/bin/bash
#SBATCH --job-name=minnvt
#SBATCH --output=minnvt.out
#SBATCH --error=minnvt.err
#SBATCH --time=1-00:00:00               # Max Time
#SBATCH --partition=gpubase_bygpu_b3  # Look up using the command: sinfo -o "%P %D %C %m %l" (-o for clear output)
#SBATCH --nodes=1
#SBATCH --ntasks=1                    # 1 MPI rank
#SBATCH --cpus-per-task=8             # Use 8 threads
#SBATCH --mem=32G                     # RAM
#SBATCH --gres=gpu:1 				  #	Request 1 GPU
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=a00836485@tec.mx

# Load GROMACS (CPU version)
module purge
module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 cuda/12.2 gromacs/2024.4          # Load gromacs and other

# Move to output dir to save files there
cd /home/aenovt/scratch/minnvt

# Get the threads from the variable cpus per task
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

# Run the simulation - EM
srun gmx mdrun -deffnm em -ntomp $OMP_NUM_THREADS -ntmpi 1
srun gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
srun gmx mdrun -deffnm nvt -ntomp $OMP_NUM_THREADS -ntmpi 1

