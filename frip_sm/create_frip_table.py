import pandas as pd
import deeptools.countReadsPerBin as crpb
import pysam
pysam.set_verbosity(0)

import plotly.io as pio
import plotly.graph_objects as go
pio.templates["mytheme"] = go.layout.Template(
    layout=go.Layout(
        # title_font=dict(family="Arial, sans-serif", size=24, color="RebeccaPurple"),
        font=dict(size=18, color="Black"),
        # legend=dict(bgcolor="LightSteelBlue", title_font=dict(size=20), font=dict(size=18)),
        # plot_bgcolor='ivory',
        # paper_bgcolor='ivory',
    )
)
pio.templates.default = "simple_white+mytheme"

######## Parameters ####################################################
dataset = 'Haarhuis_2017' 
path_to_metadata = f'/home1/yxiao977/sc1/frip_sm_data/frip_result/{dataset}/metadata.txt' 
ctcf = 'SRR5266528' # SRUN of ctcf
condition = 'SCC4KO'
peak_protein = 'CTCF'
genome_size = 3.1*10**9

df = pd.read_table(path_to_metadata)
samples = df[df['Experiment'].str.contains(condition) & ~(df['Antibody'].str.contains(peak_protein))].reset_index(drop=True)
########################################################################

# Helper function 
def calculate_frip(bam, bed, nproc=40):
    alignment = pysam.AlignmentFile(bam)
    reads_counter = crpb.CountReadsPerBin([bam], bedFile=bed, numberOfProcessors=nproc)
    reads_at_peaks = reads_counter.run()
    total_reads = alignment.mapped
    total_reads_at_peaks = reads_at_peaks.sum(axis=0)
    frip = float(total_reads_at_peaks[0])/total_reads

    return frip, reads_at_peaks, total_reads 

samples = pd.concat([df[df['SRUN'] == ctcf], samples]).reset_index(drop=True)
sruns = samples['SRUN'].to_list()

bed = f'/scratch1/yxiao977/frip_sm_data/ChIP_fastqs_maps/{dataset}/{ctcf}/{ctcf}.q30.dedup_peaks.narrowPeak'
peaks = pd.read_table(bed, header=None).iloc[:,:3]
peaks.columns = ['chrom', 'start', 'end']
num_peaks = [peaks.shape[0]] * len(sruns)
peaks_width = peaks['end'] - peaks['start']
total_bp_in_peaks = [peaks_width.sum()] * len(sruns)

total_reads = []
frip_enrich = []
samples_frips = []
for i, sample in enumerate(sruns):
    bam = f'/scratch1/yxiao977/frip_sm_data/ChIP_fastqs_maps/{dataset}/{sample}/{sample}.q30.dedup.bam'
    result = calculate_frip(bam, bed, nproc=49)
    samples_frips.append(result[0])
    total_reads.append(result[2])
    frip_enrich.append(result[0] / (total_bp_in_peaks[0] / genome_size))
    print(sample, 'done')

frip_df = pd.DataFrame({'FRiP': samples_frips})
extra_df = pd.DataFrame({'FRiP enrichment': frip_enrich, '#Peaks': num_peaks, 'Total #basepairs in peaks': total_bp_in_peaks, 'Total #reads': total_reads})
samples = pd.concat([frip_df, samples, extra_df], axis=1)
samples['GSM_accession'] = [ctcf] * len(samples)
samples = samples.rename(columns={'GSM_accession': 'peaks-SRA'})
samples.to_csv(f'/scratch1/yxiao977/frip_sm_data/frip_result/{dataset}/{condition}_ctcf({ctcf})_frips.txt', sep='\t', index=False)