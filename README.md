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

* Tips: To rescale bigwig file based on input(control) sample, the pipeline assume input sample file has the same file name with the corresponding ChIP sample (We will consider to use metadata table to improve this modality)

Once you set up configuration file, you can run the below command line in the terminal.

```
snakemake --use-conda --cores $Ncores --configfile config/config.yml
```
* Tips: [Number of cores] * [number of process in config.yml] <= the total number of cpus you have

After you generate bam files and bed files with the above command line, you can use `calculate_frip.py` to calculate FRiP value.
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