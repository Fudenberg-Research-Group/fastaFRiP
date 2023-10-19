# fastaFRiP

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
