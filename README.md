# BaitCapture

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

**BaitCapture** is a bioinformatics workflow for processing data obtained from targeted resistome bait-capture sequencing, built using [Nextflow](https://www.nextflow.io/).

## Overview

BaitCapture is an end-to-end bioinformatics workflow used to obtain coverage information for paired-end Illumina sequence reads against a database of target sequences.
The output is a tab-separated table of alignment coverages against each gene target for each sample.
For example:

```
target	seqlen	depth	len_cov	prop_cov	fold_cov
1__dfrA32__NG_047729.1__1	674	0	0	0	0
2__aac6__NG_052380.1__1	752	0	0	0	0
2__aac6__NG_052221.1__2	755	1860	157	0.20794701986755	2.4635761589404
2__aac6__NG_052259.1__3	755	3550	414	0.548344370860927	4.70198675496689
```

The workflow also outputs a table of read depths for each base pair within each gene target for a more granular view of coverage information.

The steps of the workflow are:

1. Report the quality of the raw sequence data using [FastQC](https://github.com/s-andrews/FastQC).
2. (Optional) Trim the raw sequence reads using [Trimmomatic](https://github.com/usadellab/Trimmomatic).
3. (Optional) Decontaminate the trimmed sequence reads using a host reference genome with [Bowtie2](https://github.com/BenLangmead/bowtie2).
4. Align trimmed and/or decontaminated reads against the database of gene targets using [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2).
5. Obtain sequence coverage depth statistics from the alignments and save tables in TSV format using [AlignCov](https://github.com/pcrxn/aligncov).

## Usage

If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
with `-profile test` before running the workflow on actual data.

Please provide pipeline parameters via the CLI or Nextflow `-params-file` option.
Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

BaitCapture can be run using two different input types:

- A **samplesheet**, including sample names and paths to paired-end gzipped FASTQ files, or
- A **folder** containing paired-end gzipped FASTQ files.

### Input type: Samplesheet

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a gzipped FASTQ file.

Now, you can run the pipeline using:

```bash
nextflow run OLC-LOC-Bioinformatics/BaitCapture \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --targets targets.fa \
   --outdir <OUTDIR>
```

### Input type: Folder

Instead of a samplesheet, the user can instead provide a path to a directory containing gzipped FASTQ files.
In this case, the sample name will be the name of the file up until the first period (`.`).

For example, for a folder `data/` that looks as follows:

```bash
data
├── ERR9958133_R1.fastq.gz
├── ERR9958133_R2.fastq.gz
├── ERR9958134_R1.fastq.gz
└── ERR9958134_R2.fastq.gz
```

The pipeline can be run using:

```bash
nextflow run OLC-LOC-Bioinformatics/BaitCapture \
   -profile <docker/singularity/.../institute> \
   --input_folder data/ \
   --targets targets.fa \
   --outdir <OUTDIR>
```

If the names of the gzipped FASTQ files do not end with `.fastq.gz`, an alternate extension can be specified using `--extension`.

## Contributions and support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) initative, including [nf-core/mag](https://github.com/nf-core/mag) and [nf-core/ampliseq](https://github.com/nf-core/ampliseq), and reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).


> Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x. PubMed PMID: 32055031.

In addition, references of tools and data used in this pipeline are as follows:

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

  > Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data [Online].

- [MultiQC](https://pubmed.ncbi.nlm.nih.gov/27312411/)

  > Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016;32(19):3047-8. doi: 10.1093/bioinformatics/btw354. Epub 2016 Jun 16. PubMed PMID: 27312411; PubMed Central PMCID: PMC5039924.

- [Bowtie2](https://www.nature.com/articles/nmeth.1923)

  > Langmead, B, Salzberg, S. Fast gapped-read alignment with Bowtie 2. Nat Methods. 2012;9:357-359. doi: 10.1038/nmeth.1923

- [BWA-MEM2](https://ieeexplore.ieee.org/document/8820962)

  > Vasimuddin M, Misra S, Li H, Aluru S. Efficient architecture-aware acceleration of BWA-MEM for multicore systems. IEEE Parallel and Distributed Processing Symposium (IPDPS). 2019;19028010. doi: 10.1109/IPDPS.2019.00041

- [Trimmomatic](https://academic.oup.com/bioinformatics/article/30/15/2114/2390096)

  > Bolger AM, Lohse M, Usadel B. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics. 2014;30(15):2114-2120. doi: 10.1093/bioinformatics/btu170

- [AlignCov](https://github.com/pcrxn/aligncov)

  > Brown LP. (2023) AlignCov [Online].