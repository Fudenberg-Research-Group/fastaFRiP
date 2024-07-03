import subprocess
import json
import numpy as np
import pandas as pd
import yaml

with open('config/fetch_metadata_config.yml', 'r') as f:
    config = yaml.safe_load(f)

######## Parameters ####################################################
dataset = "Haarhuis_2017"
geo_accession = "GSE90994"
path_to_accessions = (
    f"/home1/yxiao977/sc1/frip_sm_data/download_fastq/{dataset}/accessions.txt"
)

# you need to check if the author upload antibody attribute
contain_antibody_attribute = True

# you need to check how the author names these attributes
condition_key = "genotype/variation"
antibody_key = "antibody"
celltype_key = "cell line"
organism_key = "organism"

output_path = f"/home1/yxiao977/sc1/frip_sm_data/frip_result/{dataset}/metadata.txt"
########################################################################

accessions = np.loadtxt(path_to_accessions, dtype="str")


def fetch_metadata(accession):
    # Use subprocess to run ffq and capture the output
    result = subprocess.run(["ffq", accession], capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error fetching data for {accession}: {result.stderr}")
        return None
    return json.loads(result.stdout)


### fetch title from SRUN accession
experi_infos = []
gsm = []
for i, a in enumerate(accessions):
    # while True:
    #     try:
    #         data = fetch_metadata(a)
    #     except SomeSpecificException:
    #         continue
    #     break
    data = fetch_metadata(a)
    title = data[a]["title"].split(";")
    experi_info = title[1].split(":")
    gsm.append(experi_info[0])
    experi_infos.append(experi_info[1])
    print(i, a, "done")

df = pd.DataFrame(
    {"SRUN": accessions, "Experiment": experi_infos, "GSM_accession": gsm}
)

### fetch attributes from GSM accession
condition = []
antibody = []
celltype = []
organism = []

# for data has complete attributes
for i, row in df.iterrows():
    a = row["GSM_accession"].replace(" ", "")
    data = fetch_metadata(a)
    attributes = data[a]["samples"][list(data[a]["samples"].keys())[0]]["attributes"]
    condition.append(attributes[condition_key])
    if contain_antibody_attribute:
        antibody.append(attributes[antibody_key])
    else:
        antibody.append(row["Experiment"].split("_")[-1])
    celltype.append(attributes[celltype_key])
    organism.append(attributes[organism_key])
    print(a, "done")

df["Condition"] = condition
df["Antibody"] = antibody
df["Celltype"] = celltype
df["Organism"] = organism

df["Peak BED"] = ["CTCF"] * len(df)
df["author_year"] = [dataset] * len(df)
df["GEO"] = [geo_accession] * len(df)

cols = [
    "Organism",
    "Celltype",
    "Condition",
    "Antibody",
    "Peak BED",
    "author_year",
    "SRUN",
    "GSM_accession",
    "GEO",
    "Experiment",
]

df[cols].to_csv(output_path, sep="\t", index=False)
