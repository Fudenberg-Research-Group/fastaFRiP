import pandas as pd
import yaml
import glob
import os
from utils import create_frip_table_from_bed

with open("config/create_frip_table_config.yml", "r") as f:
    config = yaml.safe_load(f)

GENOME_SIZE = {"hg38": 3.1 * 10**9, "mm10": 2.7 * 10**9}

######## Read parameters from config file##############################################
species = config["parameters"]["species"]
condition = config["parameters"]["condition"]
peak_protein = config["parameters"]["peak_protein"]
nproc = config["parameters"]["nproc"]

path_to_metadata = config["input"]["path_to_metadata"]
path_to_data = config["input"]["path_to_data"]
path_to_bed = config["input"]["path_to_bed"]

output_dir = config["output"]["output_directory"]
######## Create table ################################################################
try:
    genome_size = GENOME_SIZE[species]
except KeyError:
    genome_size = config["parameters"]["genome_size"]
    if genome_size is None:
        raise Exception('Please use genome_size in create_frip_table_config.yml to specify the genome size. Currently, we only include genome size of hg38 and mm10')

df = pd.read_table(path_to_metadata)
if condition == 'all':
    conditions = df.Condition.unique().tolist()
else:
    conditions = [condition]

for condition in conditions:
    samples_metadata = df[df["Condition"] == condition].reset_index(
        drop=True
    )

    if path_to_bed != "":
        beds = [path_to_bed]
        peak_protein_sruns = [""]
    else:
        peak_proteins = df[
            df["Condition"] == condition
            & (df["Antibody"].str.contains(peak_protein, case=False))
        ].reset_index(drop=True)

        if len(peak_proteins) == 0:
            raise IndexError(f"There is no {peak_protein} sample under this condition ({condition}). Please try to use <path_to_bed> parameter in config file to specify the bed file you want to use for this condition")
        peak_protein_sruns = peak_proteins["SRUN"].to_list()
        beds = [
            glob.glob(f"{path_to_data}/{p_srun}/*.narrowPeak")[0]
            for p_srun in peak_protein_sruns
        ]

    frip_dfs = []
    for i, bed in enumerate(beds):
        frip_df = create_frip_table_from_bed(
            samples_metadata,
            bed,
            path_to_data,
            genome_size,
            species,
            nproc,
            peak_protein_srun=peak_protein_sruns[i],
        )
        frip_dfs.append(frip_df)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    frip_table = pd.concat(frip_dfs, axis=0)
    frip_table.to_csv(
        output_dir + f"/{condition}_{peak_protein}_frips.txt", sep="\t", index=False
    )
