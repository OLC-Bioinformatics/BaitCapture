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
4. Report the quality of the pre-processed sequence data using [FastQC](https://github.com/s-andrews/FastQC).
5. Align trimmed and/or decontaminated reads against the database of gene targets using [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2), [BWA](https://github.com/lh3/bwa), or [KMA](https://bitbucket.org/genomicepidemiology/kma).
6. Obtain sequence coverage and depth statistics from the alignments and save tables in TSV format using [AlignCov](https://github.com/pcrxn/aligncov).
7. Obtain further alignment statistics using [SAMtools](https://github.com/samtools/samtools).
8. Create a summary report with [MultiQC](https://github.com/ewels/MultiQC).

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

## Testing the workflow

To check if BaitCapture, Nextflow, and your container manager have been configured properly, a test run of the workflow can be performed by first cloning the GitHub repository and then running the test workflow:

```bash
git clone https://github.com/OLC-LOC-Bioinformatics/BaitCapture
cd BaitCapture/
nextflow run . \
  -profile test,<docker/singularity/.../institute> \
  --outdir <OUTDIR>
```

If your `<OUTDIR>` was `test-results`, you could then run the following command to inspect the SAMtools alignment summary statistics for the test sample:

```bash
cat test-results/samtools_stats/bwamem2/SRR14739083.stats | grep ^SN | cut -f 2-
```

The expected output for this is saved under `assets/SRR14739083.stats`.

## Running the workflow on high-performance compute clusters

Nextflow is capable of running several jobs in parallel using job submission managers (e.g. SLURM) that have been configured on high-performance compute (HPC) clusters.
For your convenience, profiles have been added to simplify running the workflow on commonly used clusters.

### Waffles

To run BaitCapture on National Microbiology Laboratory's HPC cluster Waffles using Singularity, use the following command:

```bash
nextflow run OLC-LOC-Bioinformatics/BaitCapture \
  -profile waffles \
  --outdir <OUTDIR>
```

## Advanced usage

More usage information can be obtained at any time by running `nextflow run OLC-LOC-Bioinformatics/BaitCapture --help`:

```
$ nextflow run . --help
N E X T F L O W  ~  version 23.04.2
Launching `./main.nf` [modest_picasso] DSL2 - revision: f078cc2e0b


------------------------------------------------------
  olc/baitcapture v1.0dev
------------------------------------------------------
Typical pipeline command:

  nextflow run olc/baitcapture --input samplesheet.csv --targets targets.fa -profile docker

Input/output options
  --targets                          [string]  Path to FASTA file of gene targets for alignment.
  --outdir                           [string]  The output directory where the results will be saved. You have to use absolute paths to storage on Cloud 
                                               infrastructure. 
  --input                            [string]  Path to comma-separated file containing information about the samples in the experiment.
  --input_folder                     [string]  Path to folder containing paired-end gzipped FASTQ files.
  --email                            [string]  Email address for completion summary.
  --multiqc_title                    [string]  MultiQC report title. Printed as page header, used for filename if not otherwise specified.

Workflow execution options
  --aligner                          [string]  Alignment tool to use for aligning (preprocessed) reads to the provided database of gene targets). (accepted: 
                                               bwamem2, kma, bwa) [default: bwamem2] 
  --extension                        [string]  Naming of sequencing files. [default: /*.fastq.gz]
  --host                             [string]  Path to FASTA file of host genome to use for host DNA removal (decontamination).
  --skip_trimmomatic                 [boolean] Indicate whether to skip trimming of raw reads.
  --trimmomatic                      [string]  Trimmomatic parameters. [default: ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 
                                               MINLEN:36] 

Generic options
  --multiqc_methods_description      [string]  Custom MultiQC yaml file containing HTML including a methods description.

 !! Hiding 23 params, use --validationShowHiddenParams to show them !!
------------------------------------------------------
If you use olc/baitcapture for your analysis please cite:

* The nf-core framework
  https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
  https://github.com/olc/baitcapture/blob/master/CITATIONS.md
------------------------------------------------------
```

## Contributions and support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) initative, including [nf-core/mag](https://github.com/nf-core/mag) and [nf-core/ampliseq](https://github.com/nf-core/ampliseq), and reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x. PubMed PMID: 32055031.