import os
import pandas as pd

data_dir = "/home/zuranski"
input_a = pd.read_csv(f"{data_dir}/smiles.csv")
input_asf = pd.read_csv(f"{data_dir}/smiles_sulfonate_esters.csv")

# input_a = input_a[input_a['subtype']=='alc'][['can', 'name']]
input_asf['name'] = input_asf['alcohol'] + "_" + input_asf['sf']
# input_asf = input_asf[['can', 'name']]

input_tot = input_asf

for _, row in input_tot.iterrows():
    name, smiles = row['name'], row['can']

    job_file = f"/tigress/zuranski/jobs/{name}.sh"

    file_str = f"#!/bin/bash\n" \
               f"#SBATCH --nodes=1\n" \
               f"#SBATCH --ntasks=1\n" \
               f"#SBATCH --cpus-per-task=4\n" \
               f"#SBATCH --mem-per-cpu=2G\n" \
               f"#SBATCH --time=23:59:00\n" \
               f"\n" \
               f"module purge\n" \
               f"module load anaconda3/2021.5\n" \
               f"conda activate qml\n" \
               f"\n" \
               f"script=/home/zuranski/auto-qchem/morfeus/morfeus_descriptors.py\n" \
               f"output_path=/tigress/zuranski/jobs\n" \
               f"\n"
    file_str += f"python $script '{smiles}' {name} $output_path --n_confs 200 --solvent THF --descriptors " \
                f"--substructure COS --substructure_labels C 0 S --sterimol_pairs 1,0"

    with open(job_file, "w") as f:
        f.write(file_str)

    os.system(f"sbatch {job_file}")
