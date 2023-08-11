# BaitCapture

A bioinformatics workflow for processing data obtained from targeted resistome bait-capture sequencing, built using [Nextflow](https://www.nextflow.io/).

## Future updates

- [ ] Make Trimmomatic step optional.
- [ ] Add Docker image for `aligncov.py` to allow for the use of `-profile` when executing the workflow.

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
4. Index the database of gene targets and align trimmed sequence reads against the database using [BWA MEM](https://github.com/lh3/bwa).
5. Obtain sequence coverage depth statistics from the alignments using a custom Python script and save tables in TSV format.

## Installation

Currently, the workflow must be executed from within a local environment, such as a Conda environment with the following dependencies:

```
fastqc=0.12.1
trimmomatic=0.39
multiqc=1.15
bwa=0.7.17
samtools=1.17
pandas=2.0.3
nextflow=23.04.1
```

Other package versions may work but are untested.
To create a new Conda environment for running the workflow:

1. [Install Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).
2. [Add Bioconda](https://bioconda.github.io/#usage) to your list of Conda channels and change the channel priorities.
3. Run the following command:

    ```bash
    conda create -n baitcapture \
        fastqc=0.12.1 \
        trimmomatic=0.39 \
        multiqc=1.15 \
        bwa=0.7.17 \
        samtools=1.17 \
        pandas=2.0.3 \
        nextflow=23.04.1
    ```

   - Note: For quicker package installation, I recommend the use of [`mamba`](https://mamba.readthedocs.io/en/latest/mamba-installation.html#mamba-install) within Conda.

4. Clone the BaitCapture GitHub repo:

    ```bash
    git clone https://github.com/pcrxn/BaitCapture.git
    ```

## Usage

### Quick start

After you've installed all of the workflow dependencies, BaitCapture can be run with the following commands:

```bash
cd BaitCapture/
conda activate baitcapture # or whatever you named your Conda environment
nextflow run main.nf --reads "*_R{1,2}_001.fastq.gz" --targets targets.fa
```

This command will use all paired-end sequence reads with the file name pattern of `*_R{1,2}_001.fastq.gz` as input (e.g. CL02392_R1_001.fastq.gz, CL02392_R2_001.fastq.gz).
**You must use double-quotes** for the file naming pattern to be recognized.

The sequence reads will be aligned against the FASTA file provided to `--targets`.

### Advanced usage

More usage information can be obtained at any time by running `nextflow main.nf --help`:

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