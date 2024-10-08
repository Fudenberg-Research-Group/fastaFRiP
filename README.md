# fastqFRiP

## Description
This is a pipeline for calculating FRiP ("fraction of reads in peaks") from fastq (or bed + fastq).

## Prerequisites
- python

## Installation
```
git clone https://github.com/Fudenberg-Research-Group/fastaFRiP.git
cd fastaFRiP/frip_sm
conda env create -f env/fastq_frip_env.yml -n fastq_frip_env
conda activate fastq_frip_env
```
## Getting Started

*All dependecies mentioed below are include in our conda environment, so you don't need to worry about any further installation :)*

### Downloading fastq data from GEO
FASTQ data is available in the Gene Expression Omnibus (GEO). FASTQ is a common format for storing raw sequencing data generated by next-generation sequencing technologies. In GEO, such raw sequencing data are often included as part of the supplementary files associated with a GEO Series (GSE) record. 

By clicking on the 'SRA Run Selector', users can select and download specific data (e.g., based on organism, gene, condition, or experiment type) from the Sequence Read Archive (SRA) page.

To quickly access the accession codes of ChIP-seq experiments, you can click "Metadata" button on the page to get "SraRunTable.txt" and use the following command:
```
grep "ChIP" SraRunTable.txt | awk -F, '{print $1}' > accessions.txt
```
This command extracts the codes for each file, which can later be used to download the necessary data.

We have provided a script, batch_download.sh, to facilitate the data download (In this script, we use `fasterq-dump`, which comes as part of `sra-tools`) \
You can run the script by having the 'accession.txt' in the same folder:
```
./batch_download.sh
```

### Create Bowtie2 index files
We got our index files from NCBI or UCSC genome browser. From NCBI, you can choose to use bowtie2 index files directly, or download reference genome for alignment to make your own bowtie2 index files.

<center>

|        Species/File Type|       URL Link|
|-----------:|-----------:|
| hg38/reference genome |  [Link](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz)
|
| mm39/reference genome |  [Link](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001635.9_GRCm39_full_analysis_set.fna.gz)|
|
| hg38/bowtie2 index |  [Link](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz)|
|
| mm39/bowtie2 index |  [Link](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001635.9_GRCm39_full_analysis_set.fna.bowtie_index.tar.gz)|
</center>

Most of time, you can use bowtie2 index files directly by running the following command:
```
tar -xvzf GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz
```

However, sometimes you might encounter spike-in ChIP-seq. Then you can use the following way to create a bowtie index files that include two species, here we use hg38 and mm39 as an example:
```
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001635.9_GRCm39_full_analysis_set.fna.gz

sed -i '1s/^>/>hg38_/' GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
sed -i '1s/^>/>mm39_/' GCA_000001635.9_GRCm39_full_analysis_set.fna

cat GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GCA_000001635.9_GRCm39_full_analysis_set.fna > hg38_mm39.fna
mkdir hg38_mm39
cd hg38_mm39
bowtie2-build ../hg38_mm39.fna hg38_mm39.bowtie_index
```

## Running the pipeline

### files specification
Specify the locations of your input files (FASTQ and Bowtie index files) and output files, and choose whether to include spike-in normalization in the configuration file config.yml. Detailed explanations for each parameter are included in config.yml. 
If the experiment includes spike-in, set include_spikein to true, and set index_primary and index_spikein according to the experiment.

### Metadata Table for Peak Calling:
To rescale the bigwig file and call peaks based on an input (control) sample, the pipeline requires a metadata table with two columns: ChIP and Input. The sample names should match those in the FASTQ files (e.g., SRR5085155.fastq). If a sample does not have input, exclude it from the table.
A python script for generating such a table is provided here, create_frip_table.py.

<center>

|        ChIP|       Input|
|-----------:|-----------:|
| SRR5085155 | SRR5085156 |
| SRR5085157 | SRR5085158 |
</center>

### generating BAM/BED files

Once the configuration file is set up, run the following command in the terminal to generate the required BAM/BED files:

```
snakemake --use-conda --cores $Ncores --configfile config/config.yml
```
Ensure that your computing resources are available.\
  Tips: [Number of cores] = [number of jobs] * [number of process in config.yml]. And, [Number of cores] <= the total number of cpus you have

### creating metadata 
to create metadata file, run
```
python fetch_metadata.py config/fetch_metadata_config.yml
```
after modifying `config/fetch_metadata_config.yml`. 
Example metadata table:

<center>

| SRUN       | Experiment                | GSM_accession | Treatment  | Antibody | Celltype | Organism    | Peak BED | author_year   | GEO       |
|------------|---------------------------|---------------|------------|----------|----------|-------------|----------|---------------|-----------|
| SRR5266522 | Hap1 IgG ChIPseq          | GSM2493874    | WT         | IgG      | Hap1     | Homo sapiens| CTCF     | Haarhuis_2017 | GSE90994  |
| SRR5266523 | WaplKO_3.3 IgG ChIPseq    | GSM2493875    | WaplKO_3.3 | IgG      | Hap1     | Homo sapiens| CTCF     | Haarhuis_2017 | GSE90994  |
</center>

### Create FRiP table
After you generate bam files and bed files with the above command line, you can specify path to the bed file, input data, and output data in the config file `config/create_frip_table_config.yml` , and use `create_frip_table.py` to calculate FRiP value, 
```
python create_frip_table.py config/create_frip_table_config.yml
```
 Example FRiP table:
 
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

       
