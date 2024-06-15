# fastaFRiP

## Description
This is a pipeline for calculating FRiP ("fraction of reads in peaks") from fasta (or bed + fasta).

## Prerequisites
- python
- Snakemake

## Installation
```
git clone https://github.com/Fudenberg-Research-Group/fastaFRiP.git
cd fastaFRiP/frip_sm
```
## Getting Started

You will need to specifiy where your fastq and bowtie index files located, and choose whether to include spike-in normalization procedure, and specify your output directory path in the configuration file `config.yml` (the explanation of each parameter is included in config.yml)

* Tips: To rescale bigwig file and call peaks based on input(control) sample, the pipeline will require a metadata table (stored in .txt file) have the two columns as the below example (column names must be the same as the example), and the sample name should be fetched from the fastq files (e.g. SRR5085155.fastq), and if a sample does not have input, you don't include it in the table:

<center>

|        ChIP|       Input|
|-----------:|-----------:|
| SRR5085155 | SRR5085156 |
| SRR5085157 | SRR5085158 |
</center>

Once you set up configuration file, you can run the below command line in the terminal.

```
snakemake --use-conda --cores $Ncores --configfile config/config.yml
```
* Tips: [Number of cores] * [number of process in config.yml] <= the total number of cpus you have

### Create FRiP table
After you generate bam files and bed files with the above command line, you can first use `fetch_metadata.py` to create a metadata table. Example metadata table:
<center>

| SRUN       | Experiment                | GSM_accession | Treatment  | Antibody | Celltype | Organism    | Peak BED | author_year   | GEO       |
|------------|---------------------------|---------------|------------|----------|----------|-------------|----------|---------------|-----------|
| SRR5266522 | Hap1 IgG ChIPseq          | GSM2493874    | WT         | IgG      | Hap1     | Homo sapiens| CTCF     | Haarhuis_2017 | GSE90994  |
| SRR5266523 | WaplKO_3.3 IgG ChIPseq    | GSM2493875    | WaplKO_3.3 | IgG      | Hap1     | Homo sapiens| CTCF     | Haarhuis_2017 | GSE90994  |
</center>

And then you can use `create_frip_table.py` to create a FRiP table. Example FRiP table:
<center>

| FRiP              | Organism      | Celltype | Treatment | Antibody | Peak BED | author_year   | SRUN       | peaks-SRA   | GEO       | Experiment              | FRiP enrichment | #Peaks | Total #basepairs in peaks | Total #reads |
|-------------------|---------------|----------|-----------|----------|----------|---------------|------------|-------------|-----------|-------------------------|-----------------|--------|----------------------------|--------------|
| 0.11113138926362592 | Homo sapiens | Hap1     | SCC4KO    | CTCF     | CTCF     | Haarhuis_2017 | SRR5266528 | SRR5266528  | GSE90994  | SCC4KO CTCF ChIPseq     | 27.17470160067354 | 37415  | 12677501                   | 19977713     |
| 0.00698026098917694 | Homo sapiens | Hap1     | SCC4KO    | IgG      | CTCF     | Haarhuis_2017 | SRR5266524 | SRR5266528  | GSE90994  | SCC4KO IgG ChIPseq      | 1.7068670762832925 | 37415  | 12677501                   | 14485275     |
</center>

### Only calculating FRiP value
you can use `calculate_frip.py` to calculate FRiP value.
```
python calculate_frip.py --nproc [number of cpus] [pathway to the metadata table]
```
For example:
```
python calculate_frip.py --nproc 45 /home1/yxiao977/sc1/frip_sm_data/frip_result/Hansen2017/metadata.txt
```
* The metadata.txt should looks like this: each row indicates one comparison pair, and the FRiP value will be concatenated to this metadata table as a column after run calculate_frip.py

<center>

|                BAM|                BED|
|------------------:|------------------:|
| path to bam file1 | path to bed file1 |
| path to bam file2 | path to bed file2 |

</center>

* If you use this snakemake pipeline to process fastq files, then the BAM files you want to use are with the suffix "q30.dedup.bam", and the BED files you want to use are with the suffix ".q30.dedup_peaks.narrowPeak"

## Log
What: Goal would be to have a Pipeline for calculating FRiP ("fraction of reads in peaks") from fasta (or bed + fasta). 

Why: This will be helpful since many experiments share fasta & bed files on GEO, but not many share their bam files (which are need for FRiP).
We would like to compare FRiPs from simulations with a variety of experiments. We can make the rough steps used by the Nora group more reproducible.

Rough steps are:
- bowtie2 alignment --> filtering --> bam generation
- macs2 --> bed peaks --> consensus peaks (+ additional filtering)
- deeptools.countReadsPerBin ( bam, bed) --> FRiP 

Some design questions:
- [ ] what is good test data for developing/testing a pipeline? (e.g. synthetic, or a very small real dataset).
- [ ] how should functions be organized? should snakemake be used or something else? Is the distiller-sm a good template for where to start?

Some todos (list can expand):
- [ ] files to specify requirements & environment needed
- [ ] options for pipeline to run on regular ChIP or spike-in ChIP data
       