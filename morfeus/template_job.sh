#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:01:00

module purge
module load anaconda3/2021.5
conda activate qml

script=/home/zuranski/auto-qchem/morfeus/morfeus_descriptors.py
output_path=/tigress/zuranski/jobs

python $script CCCC butane $output_path --n_confs=5 --descriptors --solvent THF