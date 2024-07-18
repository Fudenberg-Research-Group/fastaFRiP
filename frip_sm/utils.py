import pysam
import pandas as pd
import os
import subprocess
import json
import bioframe as bf
import deeptools.countReadsPerBin as crpb

pysam.set_verbosity(0)


# Helper function
def fetch_metadata(accession):
    # Use subprocess to run ffq and capture the output
    result = subprocess.run(["ffq", accession], capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error fetching data for {accession}: {result.stderr}")
        return None
    return json.loads(result.stdout)


def calculate_frip(bam, bed, nproc=40):
    alignment = pysam.AlignmentFile(bam)
    reads_counter = crpb.CountReadsPerBin([bam], bedFile=bed, numberOfProcessors=nproc)
    reads_at_peaks = reads_counter.run()
    total_reads = alignment.mapped
    total_reads_at_peaks = reads_at_peaks.sum(axis=0)
    frip = float(total_reads_at_peaks[0]) / total_reads

    return frip, reads_at_peaks, total_reads


def create_frip_table_from_bed(
    samples_metadata,
    path_to_bed,
    path_to_data,
    genome_size,
    species,
    nproc,
    peak_protein_srun="",
):
    sruns = samples_metadata["SRUN"].to_list()

    peaks = bf.read_table(path_to_bed, schema="bed").iloc[:, :3]
    num_peaks = [peaks.shape[0]] * len(sruns)
    peaks_width = peaks["end"] - peaks["start"]
    total_bp_in_peaks = [peaks_width.sum()] * len(sruns)

    if os.path.exists(f"{path_to_data}/{sruns[0]}/{sruns[0]}.q30.{species}.sort.bam"):
        suffix = f"{species}.sort"
    else:
        suffix = "dedup"

    total_reads = []
    frip_enrich = []
    samples_frips = []
    for i, sample in enumerate(sruns):
        bam = f"{path_to_data}/{sample}/{sample}.q30.{suffix}.bam"
        result = calculate_frip(bam, path_to_bed, nproc=nproc)
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
    frip_df = pd.concat([frip_df, samples_metadata, extra_df], axis=1)
    frip_df["GSM_accession"] = [peak_protein_srun] * len(frip_df)
    frip_df = frip_df.rename(columns={"GSM_accession": "peaks-SRA"})

    return frip_df