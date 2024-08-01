import numpy as np
import pandas as pd
import os
import yaml
import argparse
from utils import fetch_metadata

parser = argparse.ArgumentParser(description='Fetch a metadata table')
parser.add_argument('config_path', type=str, help='Path to the configuration file.')
args = parser.parse_args()

with open(args.config_path, "r") as f:
    config = yaml.safe_load(f)

######## Parameters ####################################################
dataset = config["parameters"]["dataset"]
geo_accession = config["parameters"]["geo_accession"]

path_to_accessions = config["input"]["path_to_accessions"]

output_dir = config["output"]["output_directory"]
########################################################################

accessions = np.loadtxt(path_to_accessions, dtype="str")

### fetch title from SRUN accession
experi_infos = []
gsm = []
for i, a in enumerate(accessions):
    while True:
        try:
            data = fetch_metadata(a)
            title = data[a]["title"].split(";")
        except TypeError:
            continue
        break
    experi_info = title[1].split(":")
    gsm.append(experi_info[0])
    experi_infos.append(experi_info[1])
    print(i, a, "title has been fetched")

df = pd.DataFrame(
    {"SRUN": accessions, "Experiment": experi_infos, "GSM_accession": gsm}
)

### fetch attributes from GSM accession
condition = []
antibody = []
celltype = []
organism = []

KEY_WORDS = ["gen", "cell", "organism"]
attribute_keys = []
for i, row in df.iterrows():
    a = row["GSM_accession"].replace(" ", "")
    data = fetch_metadata(a)
    attributes = data[a]["samples"][list(data[a]["samples"].keys())[0]]["attributes"]

    if len(attribute_keys) == 0:
        key_df = pd.DataFrame({"Keys": attributes.keys()})
        for kw in KEY_WORDS:
            try:
                attribute_keys.append(
                    key_df[key_df["Keys"].str.contains(kw, case=False)].iloc[0, 0]
                )
            except IndexError:
                print(f"There is no attributes match with the substring {kw}")

    condition.append(attributes[attribute_keys[0]])
    celltype.append(attributes[attribute_keys[1]])
    organism.append(attributes[attribute_keys[2]])

    try:
        antibody_key = key_df[key_df["Keys"].str.contains("anti", case=False)].iloc[
            0, 0
        ]
        antibody.append(attributes[antibody_key])
    except IndexError:
        antibody.append(row["Experiment"].split("_")[-1])
    print(i, accessions[i], "attributes have been fetched")

df["Condition"] = condition
df["Antibody"] = antibody
df["Celltype"] = celltype
df["Organism"] = organism

df["Peak_ChIP"] = ["N/A"] * len(df)
df["author_year"] = [dataset] * len(df)
df["GEO"] = [geo_accession] * len(df)

cols = [
    "Organism",
    "Celltype",
    "Condition",
    "Antibody",
    "Peak_ChIP",
    "author_year",
    "SRUN",
    "GSM_accession",
    "GEO",
    "Experiment",
]

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

df[cols].to_csv(output_dir + '/metadata.txt', sep="\t", index=False)
