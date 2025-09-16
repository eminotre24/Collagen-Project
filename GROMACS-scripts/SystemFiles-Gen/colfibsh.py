# Libraries
import os
import datetime
import shutil

def generate_files_prep(name, pdb_file, ff, watertype, solvent_config, pion, nion, dimensions = None, angles = None, distance = None, ccon = None, chain_org = True):
    # Generate preparation files, the structure of the folder of the "project" and the steps scripts ready to jsut execute

    # Folders
    date_today = datetime.datetime.today().strftime('-%d%m')
    name_folder = name + date_today
    files_folder = name_folder + "/" + name + "-files"
    os.makedirs(name_folder, exist_ok=True)
    os.makedirs(files_folder, exist_ok=True)

    # Files
    if os.path.isdir("templates/" + ff + ".ff"):
        shutil.copytree("templates/" + ff + ".ff", files_folder + "/" + ff + ".ff")
    shutil.copy("templates/ions.mdp", files_folder)
    shutil.copy(pdb_file, files_folder)

    with open(files_folder + "/steps.sh", "w") as file:
        # Generate script
        file.write(f"gmx pdb2gmx -f {pdb_file} -o colfib.gro -ff {ff} -water {watertype}\n")

        if chain_org:
            file.write("mkdir itp_chains\n")
            file.write("mv *.itp ./itp_chains/\n")
            file.write(r"""sed -i '' -E '/^;[[:space:]]*Include chain topologies[[:space:]]*$/, /^;[[:space:]]*Include water topology[[:space:]]*$/ s@^([[:space:]]*#include ")@\1itp_chains/@' topol.top""")
            file.write("\n")
        
        if distance is None:
            file.write(f"gmx editconf -f colfib.gro -o colfib-boxd.gro -c -box {dimensions} -angles {angles}\n")
        else:
            file.write(f"gmx editconf -f colfib.gro -o colfib-boxd.gro -c -d {distance}\n")

        file.write(f"gmx solvate -cp colfib-boxd.gro -cs {solvent_config} -o colfib-solvtd.gro -p topol.top\n")
        file.write("gmx grompp -f ions.mdp -c colfib-solvtd.gro -p topol.top -o ions.tpr\n")

        if ccon is None:
            file.write(f"echo 13 | gmx genion -s ions.tpr -o colfib-solvion.gro -p topol.top -pname {pion} -nname {nion} -neutral\n")
        else:
            file.write(f"echo 13 | gmx genion -s ions.tpr -o colfib-solvion.gro -p topol.top -pname {pion} -nname {nion} -neutral -conc {ccon}\n")

    print("Files generated in folder: " + name_folder)

def generate_script(name, runtime, user):
    # Generate ACENET script
    date_today = datetime.datetime.today().strftime('-%d%m')
    name_dated = name + date_today
    name_folder = name + "-files"
    script_path = f"{name_dated}/{name}-script.sh"

    # Script content
    script_content = f"""#!/bin/bash
                        #SBATCH --job-name={name_dated}
                        #SBATCH --output=/home/aenovt/out_err/{name_dated}/out.out
                        #SBATCH --error=/home/aenovt/out_err/{name_dated}/err.err
                        #SBATCH --time={runtime}
                        #SBATCH --nodes=1
                        #SBATCH --ntasks=1
                        #SBATCH --cpus-per-task=12
                        #SBATCH --gres=gpu:1
                        #SBATCH --mem-per-cpu=2G
                        #SBATCH --mail-type=BEGIN,END,FAIL
                        #SBATCH --mail-user={user}

                        module purge
                        module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 cuda/12.2 gromacs/2024.4

                        cd /home/aenovt/scratch/{name_folder}
                        export OMP_NUM_THREADS="${{SLURM_CPUS_PER_TASK:-1}}"

                        # Energy Minimization Phase
                        srun gmx grompp -f minimize.mdp -c colfib-solvion.gro -p topol.top -o em.tpr
                        srun gmx mdrun -deffnm em -ntomp $OMP_NUM_THREADS -ntmpi $SLURM_NTASKS -nb gpu

                        # Heating Phase
                        srun gmx grompp -f nvt_heating.mdp -c em.gro -r em.gro -p topol.top -o nvt_heating.tpr
                        srun gmx mdrun -deffnm nvt_heating -ntomp $OMP_NUM_THREADS -ntmpi $SLURM_NTASKS -nb gpu -pme gpu -update gpu -bonded cpu

                        # NVT Equilibration
                        srun gmx grompp -f nvt_equilibration.mdp -c nvt_heating.gro -r nvt_heating.gro -t nvt_heating.cpt -p topol.top -o nvt_equilibration.tpr
                        srun gmx mdrun -deffnm nvt_equilibration -ntomp $OMP_NUM_THREADS -ntmpi $SLURM_NTASKS -nb gpu -pme gpu -update gpu -bonded cpu

                        # NPT Phase
                        srun gmx grompp -f npt.mdp -c nvt_equilibration.gro -r nvt_equilibration.gro -t nvt_equilibration.cpt -p topol.top -o npt.tpr
                        srun gmx mdrun -deffnm npt -ntomp $OMP_NUM_THREADS -ntmpi $SLURM_NTASKS -nb gpu -pme gpu -update gpu -bonded cpu

                        # Production Run
                        srun gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_1.tpr -maxwarn 2
                        srun gmx mdrun -deffnm md_1 -ntomp $OMP_NUM_THREADS -ntmpi $SLURM_NTASKS -nb gpu -pme gpu -update gpu -bonded cpu
                        """

    # Write the script to file
    with open(script_path, "w") as f:
        f.write(script_content)


    print("SLURM script generated with name: " + script_path)

def mdp_parms(temperature , pressure, taup, nvt_time, npt_time, md_time, dt_eq, dt_md):
    # Edit configuration of mdp files
    

    print("MDP files generated")