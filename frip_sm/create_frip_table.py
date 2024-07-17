import pandas as pd
import yaml
import glob
from utils import create_frip_table_from_bed

with open("config/create_frip_table_config.yml", "r") as f:
    config = yaml.safe_load(f)

GENOME_SIZE = {"hg38": 3.1 * 10**9, "mm10": 2.7 * 10**9}

######## Read parameters from config file##############################################
species = config["parameters"]["species"]
condition = config["parameters"]["condition"]
peak_protein = config["parameters"]["peak_protein"]
nproc = config["parameters"]["nproc"]
path_to_bed = config["parameters"]["path_to_bed"]
path_to_metadata = config["input"]["path_to_metadata"]
path_to_data = config["input"]["path_to_data"]
output_dir = config["output"]["output_directory"]

######## Create table ################################################################
genome_size = GENOME_SIZE[species]
df = pd.read_table(path_to_metadata)
samples_metadata = df[df["Condition"].str.contains(condition, case=False)].reset_index(
    drop=True
)

if path_to_bed != "":
    beds = [path_to_bed]
    peak_protein_sruns = [""]
else:
    peak_proteins = df[
        df["Condition"].str.contains(condition, case=False)
        & (df["Antibody"].str.contains(peak_protein, case=False))
    ].reset_index(drop=True)
    peak_protein_sruns = peak_proteins["SRUN"].to_list()
    beds = [
        glob.glob(f"{path_to_data}/{p_srun}/*.narrowPeak")[0]
        for p_srun in peak_protein_sruns
    ]

frip_dfs = []
for i, bed in enumerate(beds):
    frip_df = create_frip_table_from_bed(
        samples_metadata,
        path_to_bed,
        path_to_data,
        genome_size,
        species,
        nproc,
        peak_protein_srun=peak_protein_sruns[i],
    )
    frip_dfs.append(frip_df)

frip_table = pd.concat(frip_dfs, axis=0)
frip_table.to_csv(
    output_dir + f"/{condition}_{peak_protein}_frips.txt", sep="\t", index=False
)
