import pandas as pd
import yaml
from utils import calculate_frip

with open('config/create_frip_table_config.yml', 'r') as f:
    config = yaml.safe_load(f)

GENOME_SIZE = {'human':3.1 * 10**9, 'mouse':2.7 * 10**9}

######## Read parameters from config file##############################################
dataset = "Haarhuis_2017"
species = 'human'
path_to_metadata = (
    f"/home1/yxiao977/sc1/frip_sm_data/frip_result/{dataset}/metadata.txt"
)
path_to_data = f"/scratch1/yxiao977/frip_sm_data/ChIP_fastqs_maps/{dataset}"

nproc = 49
condition = "SCC4KO"
peak_protein = "CTCF"
path_to_bed = ""

output_path = f"/scratch1/yxiao977/frip_sm_data/frip_result/{dataset}/{condition}_{peak_protein}_frips.txt"

######## Create table ################################################################
genome_size = GENOME_SIZE[species]
df = pd.read_table(path_to_metadata)
samples = df[
        df["Condition"].str.contains(condition, case=False)
    ].reset_index(drop=True)

if path_to_bed != "":
    beds = [path_to_bed]
    peak_protein_sruns = [""]
else:
    peak_proteins = df[
        df["Condition"].str.contains(condition, case=False)
        & (df["Antibody"].str.contains(peak_protein, case=False))
    ].reset_index(drop=True)
    peak_protein_sruns = peak_proteins["SRUN"].to_list()
    beds = [f"{path_to_data}/{p_srun}/{p_srun}.q30.dedup_peaks.narrowPeak" for p_srun in peak_protein_sruns]
    
frip_dfs = []
for i, bed in enumerate(beds):
    sruns = samples["SRUN"].to_list()
    peaks = pd.read_table(bed, header=None).iloc[:, :3]
    peaks.columns = ["chrom", "start", "end"]
    num_peaks = [peaks.shape[0]] * len(sruns)
    peaks_width = peaks["end"] - peaks["start"]
    total_bp_in_peaks = [peaks_width.sum()] * len(sruns)

    total_reads = []
    frip_enrich = []
    samples_frips = []
    for sample in sruns:
        bam = f"{path_to_data}/{sample}/{sample}.q30.dedup.bam"
        result = calculate_frip(bam, bed, nproc=nproc)
        samples_frips.append(result[0])
        total_reads.append(result[2])
        frip_enrich.append(result[0] / (total_bp_in_peaks[0] / genome_size))
        print(sample, "done")

    frip_df = pd.DataFrame({"FRiP": samples_frips})
    extra_df = pd.DataFrame(
        {
            "FRiP enrichment": frip_enrich,
            "#Peaks": num_peaks,
            "Total #basepairs in peaks": total_bp_in_peaks,
            "Total #reads": total_reads,
        }
    )
    frip_df = pd.concat([frip_df, samples, extra_df], axis=1)
    frip_df["GSM_accession"] = [peak_protein_sruns[i]] * len(frip_df)
    frip_df = frip_df.rename(columns={"GSM_accession": "peaks-SRA"})
    frip_dfs.append(frip_df)

frip_table = pd.concat(frip_dfs, axis=0)
frip_table.to_csv(output_path, sep="\t", index=False)