# BaitCapture

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

A bioinformatics workflow for processing data obtained from targeted resistome bait-capture sequencing, built using [Nextflow](https://www.nextflow.io/).

## Future updates

- [ ] Add host DNA decontamination step.
- [ ] Make Trimmomatic step optional.

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
2. Trim the raw sequence reads using [Trimmomatic](https://github.com/usadellab/Trimmomatic).
3. Summarize the FastQC reports and Trimmomatic logs using [MultiQC](https://multiqc.info/).
4. Index the database of gene targets and align trimmed sequence reads against the database using [BWA MEM](https://github.com/lh3/bwa) and [KMA](https://bitbucket.org/genomicepidemiology/kma).
5. Obtain sequence coverage depth statistics from the alignments and save tables in TSV format using [AlignCov](https://github.com/pcrxn/aligncov).

## Installation

1. Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) (22.10.1+).
2. Install any of [Docker](https://docs.docker.com/engine/install/), [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/), or [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) package manager.


3. Clone the BaitCapture GitHub repo:

    ```bash
    git clone https://github.com/pcrxn/BaitCapture.git
    ```

4. Run your analysis with your preferred package manager (`-profile <docker,singularity,conda>`).

## Usage

### Quick start

After you've installed all of the workflow dependencies, BaitCapture can be run with the following commands:

```bash
cd BaitCapture/
nextflow run main.nf -profile singularity --reads "*_R{1,2}_001.fastq.gz" --targets targets.fa
```

This command will:

  - Use Singularity to manage dependencies.
  - Use all paired-end sequence reads with the file name pattern of `*_R{1,2}_001.fastq.gz` as input (e.g. CL02392_R1_001.fastq.gz, CL02392_R2_001.fastq.gz).
    - Important: **You must use double-quotes** for the file naming pattern to be recognized.
  - Align paired-end sequence reads to the FASTA file `targets.fa`.

### Advanced usage

More usage information can be obtained at any time by running `nextflow run main.nf --help`:

```
Typical pipeline command:

  nextflow run main.nf --reads "*_R{1,2}_001.fastq.gz" --targets targets.fa

Input/output options
  --reads       [string]  A naming pattern for the .fastq.gz files which will be aligned to the targets FASTA file.
  --outdir      [string]  The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. [default: 
                          results/] 

Other parameters
  --targets     [string]  The FASTA file containing the DNA sequences used for bait-capture sequencing.
  --trimmomatic [string]  Command-line arguments for custom Trimmomatic parameters.
```
