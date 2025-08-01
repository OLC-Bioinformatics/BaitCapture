<!-- omit in toc -->
# BaitCapture

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13737733.svg)](https://doi.org/10.5281/zenodo.13737733)

<p align='center'><img src='assets/baitcapture-banner_v02.png' alt="BaitCapture banner" width="75%"></p>

**BaitCapture** is a bioinformatics workflow designed for processing sequencing data obtained from targeted resistome bait-capture sequencing, built using [Nextflow](https://www.nextflow.io/).

<!-- omit in toc -->
## Contents

- [Introduction](#introduction)
- [Quick start](#quick-start)
- [Pipeline summary](#pipeline-summary)
- [Usage](#usage)
  - [Input type: Samplesheet](#input-type-samplesheet)
  - [Input type: Folder](#input-type-folder)
    - [Case #1: Default file name pattern](#case-1-default-file-name-pattern)
    - [Case #2: Alternate file name pattern](#case-2-alternate-file-name-pattern)
  - [(Optional) Input type: Target metadata](#optional-input-type-target-metadata)
- [Output](#output)
- [Testing the workflow](#testing-the-workflow)
- [Running the workflow on high-performance compute clusters](#running-the-workflow-on-high-performance-compute-clusters)
  - [Waffles](#waffles)
- [Advanced usage](#advanced-usage)
- [Contributions and support](#contributions-and-support)
- [Citations](#citations)

## Introduction

BaitCapture is based upon a Bash script workflow originally created by [Shay et al. (2023)](http://doi.org/10.1186/s40793-023-00482-0).
Though it was designed in particular consideration of bait-capture sequencing data, BaitCapture can be used for any paired-end Illumina sequencing dataset where the user needs to align many sequence reads to a reference database of gene targets.

BaitCapture offers the following features:

- **Quality control**: Assess the quality of raw and pre-processed sequence data.
- **Pre-processing**:
  - Read decontamination using a host reference genome
  - Quality-based trimming
  - Adapter removal
- **Read alignment**: Align reads against a reference database of gene targets using [KMA](https://bitbucket.org/genomicepidemiology/kma), [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2), or [BWA](https://github.com/lh3/bwa).
- **[Alignment reports](docs/alignment_reports.md)**:
  - `mapstats.tsv`: A table of read alignment statistics against each gene target for each sample, including KMA-specific alignment statistics.
  - `sumstats.tsv`: A table of summary statistics for each sample, including the on-target alignment rate, and the number of reads lost from host decontamination and filtering.
  - `presence_absence.tsv`: A table of presence-absence calls for each gene target in each sample, based upon user-defined thresholds.
  - `presence_absence_clusters.tsv`: A table of presence-absence calls for each gene target cluster in each sample, with clusters defined by a target metadata file (e.g. resistance mechanism).



## Quick start

<img alt="BaitCapture terminal demo" src=assets/baitcapture-demo.gif width="100%" />

## Pipeline summary

The steps of the workflow are:

1. Report the quality of the raw sequence data using [FastQC](https://github.com/s-andrews/FastQC).
2. (Optional) Trim the raw sequence reads using [fastp](https://github.com/OpenGene/fastp).
3. (Optional) Decontaminate the trimmed sequence reads using a host reference genome with [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2).
4. Report the quality of the pre-processed sequence data using [FastQC](https://github.com/s-andrews/FastQC).
5. Obtain total read and bp counts from the raw and pre-processed sequence data using [fastq-scan](https://github.com/rpetit3/fastq-scan).
6. Align trimmed and/or decontaminated reads against the database of gene targets using [KMA](https://bitbucket.org/genomicepidemiology/kma), [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2), or [BWA](https://github.com/lh3/bwa).
7. Obtain sequence coverage and depth statistics from the alignments using [AlignCov](https://github.com/pcrxn/aligncov), [Mosdepth](https://github.com/brentp/mosdepth), and [SAMtools](https://github.com/samtools/samtools).
8. Create a [MultiQC](https://github.com/ewels/MultiQC) report and other summary reports with custom scripts.

## Usage

If you are new to Nextflow and the nf-core framework, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow.

Please provide pipeline parameters via the CLI or Nextflow `-params-file` option.
Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

BaitCapture can be run using two different input types:

- A **samplesheet**, including sample names and paths to paired-end gzipped FASTQ files, or
- A **folder** containing paired-end gzipped FASTQ files.

Additionally, a target metadata file can optionally be provided to generate a `presence_absence_clusters.tsv` file.

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
nextflow run OLC-Bioinformatics/BaitCapture \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --targets targets.fa \
   --outdir <OUTDIR>
```

### Input type: Folder

Instead of a samplesheet, the user can instead provide a path to a directory containing gzipped FASTQ files.
In this case, the sample name will be the name of the file up until the first underscore (`_`).

#### Case \#1: Default file name pattern

The default file naming pattern that `--input_folder` searches for is `"/*_R{1,2}.fastq.gz"`.
For example, for a folder `data/` containing sequencing files that looks as follows:

```bash
data
├── ERR9958133_R1_001.fastq.gz
├── ERR9958133_R2_001.fastq.gz
├── ERR9958134_R1_001.fastq.gz
└── ERR9958134_R2_001.fastq.gz
```

The workflow can be run simply by using:

```bash
nextflow run OLC-Bioinformatics/BaitCapture \
   -profile <docker/singularity/.../institute> \
   --input_folder data/ \
   --targets targets.fa \
   --outdir <OUTDIR>
```

And the sample names will be:
- ERR9958133
- ERR9958134

#### Case \#2: Alternate file name pattern

If the names of the gzipped FASTQ files **do not** end with `_R{1,2}_001.fastq.gz`, an alternate sequencing file pattern must be specified using `--pattern`.
For example, for a folder `more-data/` that looks as follows:

```bash
more-data
├── SAMN000214_R1.fastq.gz
├── SAMN000214_R2.fastq.gz
├── SAMN000215_R1.fastq.gz
└── SAMN000215_R2.fastq.gz
```

The workflow can be run using:

```bash
nextflow run OLC-Bioinformatics/BaitCapture \
   -profile <docker/singularity/.../institute> \
   --input_folder more-data/ \
   --pattern "/*_R{1,2}.fastq.gz" \
   --targets targets.fa \
   --outdir <OUTDIR>
```

And the sample names will be:
- SAMN000214
- SAMN000215

> [!NOTE]
> When providing an argument to `--pattern`, the string must be enclosed in double quotes (`""`) and must be prepended with a forward slash (`/`).

### (Optional) Input type: Target metadata

The path to a target metadata file can optionally be provided with `--target_metadata` to generate a `presence_absence_clusters.tsv` file, which describes the presence or absence of clusters of gene targets (e.g. AMR resistance mechanism, target drug class, replicon type).
The requirements for this file are:

- A comma-separated value (CSV) file with a header row, with one column header named `target`.
- At least one column in addition to the `target` column, which describes a categorical variable to cluster the targets by.
- Each value in the target column must exactly match the first word of a sequence header in the `--targets` FASTA database.

For example, if we provide a `--targets` database containing two sequences:

```fasta
>gb|X91840.1|+|1-1146|CMY-2 [Klebsiella pneumoniae]
ATGATGAAAAAATCGTTATGCTGCGCTCTGCTGCTGACAGCCTCTTTCTCCACATTTGCTGCCGCAAAAACAGAACAACAGATTGCCGAT
ATCGTTAATCGCACCATCACCCCGTTGATGCAGGAGCAGGCTATTCCGGGTATGGCCGTTGCCGTTATCTACCAGGGAAAACCCTATTAT
TTCACCTGGGGTAAAGCCGATATCGCCAATAACCACCCAGTCACGCAGCAAACGCTGTTTGAGCTAGGATCGGTTAGTAAGACGTTTAAC
GGCGTGTTGGGCGGCGATGCTATCGCCCGCGGCGAAATTAAGCTCAGCGATCCGGTCACGAAATACTGGCCAGAACTGACAGGCAAACAG
TGGCAGGGTATCCGCCTGCTGCACTTAGCCACCTATACGGCAGGCGGCCTACCGCTGCAGATCCCCGATGACGTTAGGGATAAAGCCGCA
TTACTGCATTTTTATCAAAACTGGCAGCCGCAATGGACTCCGGGCGCTAAGCGACTTTACGCTAACTCCAGCATTGGTCTGTTTGGCGCG
CTGGCGGTGAAACCCTCAGGAATGAGTTACGAAGAGGCAATGACCAGACGCGTCCTGCAACCATTAAAACTGGCGCATACCTGGATTACG
GTTCCGCAGAACGAACAAAAAGATTATGCCTGGGGCTATCGCGAAGGGAAGCCCGTACACGTTTCTCCGGGACAACTTGACGCCGAAGCC
TATGGCGTGAAATCCAGCGTTATTGATATGGCCCGCTGGGTTCAGGCCAACATGGATGCCAGCCACGTTCAGGAGAAAACGCTCCAGCAG
GGCATTGCGCTTGCGCAGTCTCGCTACTGGCGTATTGGCGATATGTACCAGGGATTAGGCTGGGAGATGCTGAACTGGCCGCTGAAAGCT
GATTCGATCATCAACGGCAGCGACAGCAAAGTGGCATTGGCAGCGCTTCCCGCCGTTGAGGTAAACCCGCCCGCCCCCGCAGTGAAAGCC
TCATGGGTGCATAAAACGGGCTCCACTGGTGGATTTGGCAGCTACGTAGCCTTCGTTCCAGAAAAAAACCTTGGCATCGTGATGCTGGCA
AACAAAAGCTATCCTAACCCTGTCCGTGTCGAGGCGGCCTGGCGCATTCTTGAAAAGCTGCAATAA
>gb|AF231986.2|+|3308-4522|floR [Escherichia coli]
ATGACCACCACACGCCCCGCGTGGGCCTATACGCTGCCGGCAGCACTGCTGCTGATGGCTCCTTTCGACATCCTCGCTTCACTGGCGATG
GATATTTATCTCCCTGTCGTTCCAGCGATGCCCGGCATCCTGAACACGACGCCCGCTATGATCCAACTCACGTTGAGCCTCTATATGGTG
ATGCTCGGCGTGGGCCAGGTGATTTTTGGTCCGCTCTCAGACAGAATCGGGCGACGGCCAATTCTACTTGCGGGCGCAACGGCTTTCGTC
ATTGCGTCTCTGGGAGCAGCTTGGTCTTCAACTGCACCGGCCTTTGTCGCTTTCCGTCTACTTCAAGCAGTGGGCGCGTCGGCCATGCTG
GTGGCGACGTTCGCGACGGTTCGCGACGTTTATGCCAACCGTCCTGAGGGTGTCGTCATCTACGGCCTTTTCAGTTCGGTGCTGGCGTTC
GTGCCTGCGCTCGGCCCTATCGCCGGAGCATTGATCGGCGAGTTCTTGGGATGGCAGGCGATATTCATTACTTTGGCTATACTGGCGATG
CTCGCACTCCTAAATGCGGGTTTCAGGTGGCACGAAACCCGCCCTCTGGATCAAGTCAAGACGCGCCGATCTGTCTTGCCGATCTTCGCG
AGTCCGGCTTTTTGGGTTTACACTGTCGGCTTTAGCGCCGGTATGGGCACCTTCTTCGTCTTCTTCTCGACGGCTCCCCGTGTGCTCATA
GGCCAAGCGGAATATTCCGAGATCGGATTCAGCTTTGCCTTCGCCACTGTCGCGCTTGTAATGATCGTGACAACCCGTTTCGCGAAGTCC
TTTGTCGCCAGATGGGGCATCGCAGGATGCGTGGCGCGTGGGATGGCGTTGCTTGTTTGCGGAGCGGTCCTGTTGGGGATCGGCGAACTT
TACGGCTCGCCGTCATTCCTCACCTTCATCCTACCGATGTGGGTTGTCGCGGTCGGTATTGTCTTCACGGTGTCCGTTACCGCGAACGGC
GCTTTGGCAGAGTTCGACGACATCGCGGGATCAGCGGTCGCGTTCTACTTCTGCGTTCAAAGCCTGATAGTCAGCATTGTCGGGACATTG
GCGGTGGCACTTTTAAACGGTGACACAGCGTGGCCCGTGATCTGTTACGCCACGGCGATGGCGGTACTGGTTTCGTTGGGGCTGGTGCTC
CTTCGGCTCCGTGGGGCTGCCACCGAGAAGTCGCCAGTCGTCTAA
```

We might optionally provide a `--target_metadata` CSV file which looks as follows:

| target                             | drug_class  | resistance_mechanism    |
|------------------------------------|-------------|-------------------------|
| gb\|X91840.1\|+\|1-1146\|CMY-2     | beta_lactam | antibiotic_inactivation |
| gb\|AF231986.2\|+\|3308-4522\|floR | phenicol    | antibiotic_efflux       |

So, BaitCapture would generate the following `presence_absence_clusters.tsv` file if only the CMY-2 target is detected in sampleA:

| sampleid | metavar_name         | metavar_value           | presence_absence |
|----------|----------------------|-------------------------|------------------|
| sampleA  | drug_class           | beta_lactam             | 1                |
| sampleA  | drug_class           | phenicol                | 0                |
| sampleA  | resistance_mechanism | antibiotic_inactivation | 1                |
| sampleA  | resistance_mechanism | antibiotic_efflux       | 0                |

## Output

The pipeline will output the following directories:

```
results/
├── fastp
├── fastqc
│   ├── preprocessed
│   └── raw
├── multiqc
│   ├── multiqc_data
│   └── multiqc_plots
├── pipeline_info
└── summary
```

- `fastp/`: Contains fastp reports for trimmed FASTQ files.
- `fastqc/`: Contains FastQC reports for raw and pre-processed FASTQ files.
- `multiqc/`: Contains MultiQC reports.
- `pipeline_info/`: Contains Nextflow logs and reports.
- `summary/`: Contains [alignment reports](docs/alignment_reports.md) for each sample.

## Testing the workflow

To check if BaitCapture, Nextflow, and your container manager have been configured properly, a test run of the workflow can be performed by first cloning the GitHub repository and then running the test workflow as follows:

```bash
git clone https://github.com/OLC-Bioinformatics/BaitCapture
cd BaitCapture/
nextflow run . \
  -profile test,<docker/singularity/.../institute> \
  --outdir <OUTDIR>
```

If your `<OUTDIR>` was `results/`, you could then run the following command to inspect the summary statistics for the test sample:

```bash
$ cat results/summary/sumstats.tsv 
sampleid        raw_total_reads raw_total_bp    fastp_total_reads       fastp_total_bp  decontam_total_reads        decontam_total_bp   mapped_total_reads      mapped_total_bp percent_reads_lost_fastp        percent_reads_lost_decontam    percent_reads_on_target
SRR14739083     553624  83043600        464776  69619648        455764  68278988        98485   14753774        16.05   1.94    21.61
```

## Running the workflow on high-performance compute clusters

Nextflow is capable of running several jobs in parallel using job submission managers (e.g. SLURM) that have been configured on high-performance compute (HPC) clusters.
For your convenience, profiles have been added to simplify running the workflow on commonly used clusters.

### Waffles

To run BaitCapture on National Microbiology Laboratory's HPC cluster Waffles using Singularity, use `-profile waffles`.
For example:

```bash
nextflow run OLC-Bioinformatics/BaitCapture \
  -profile waffles \
  --input samplesheet.csv \
  --targets targets.fa \
  --outdir <OUTDIR>
```

## Advanced usage

More usage information can be obtained at any time by running `nextflow run OLC-Bioinformatics/BaitCapture --help`:

```
$ nextflow run OLC-Bioinformatics/BaitCapture --help
N E X T F L O W  ~  version 23.10.1
Launching `https://github.com/OLC-Bioinformatics/BaitCapture` [stupefied_bassi] DSL2 - revision: 0fd15d548b [dev]


------------------------------------------------------
  olc/baitcapture v1.0.0-g0fd15d5
------------------------------------------------------
Typical pipeline command:

  nextflow run olc/baitcapture --input samplesheet.csv --targets targets.fa --outdir results/ -profile singularity

Input/output options
  --input                            [string]  Path to comma-separated file containing information about the samples in the experiment.
  --input_folder                     [string]  Path to folder containing paired-end gzipped FASTQ files.
  --pattern                          [string]  Naming of sequencing files for `--input_folder`. Must use double-quotes (`""`) and a prepended slash (`/`). 
                                               [default: "/*_R{1,2}_001.fastq.gz"] 
  --targets                          [string]  Path to FASTA file of gene targets for alignment.
  --outdir                           [string]  The output directory where the results will be saved. You have to use absolute paths to storage on Cloud 
                                               infrastructure. 
  --host                             [string]  Path to FASTA file of host genome to use for host DNA removal (decontamination).
  --adapters                         [string]  Path to FASTA file of adapter sequences to use for adapter removal with FASTP.
  --target_metadata                  [string]  Path to comma-separated file containing information about the metadata for targets used in the 
                                               experiment. 
  --email                            [string]  Email address for completion summary.
  --multiqc_title                    [string]  MultiQC report title. Printed as page header, used for filename if not otherwise specified.

Workflow execution options
  --aligner                          [string]  Alignment tool to use for aligning (preprocessed) reads to the provided database of gene targets). (accepted: 
                                               bwamem2, kma, bwa) [default: kma] 
  --skip_trimming                    [boolean] Indicate whether to skip trimming of raw reads.
  --report_all                       [boolean] Report undetected targets in merged results files.

Target detection thresholds
  --fold_cov_threshold               [number]  The minimum fold-coverage of a target that must be achieved to call a positive detection. [default: 0.9]
  --len_cov_threshold                [integer] The minimum length (in bp) that a target must be covered by to call a positive detection. [default: 0]
  --mapped_reads_threshold           [integer] The minimum number of reads that must be mapped to a target to call a positive detection. [default: 2]
  --prop_cov_threshold               [number]  The minimum percentage of length (in bp) that a target must be covered by to call a positive detection. 
                                               [default: 0.9] 
  --pident_threshold                 [number]  The minimum percentage identity match to a target that must be achieved to call a positive detection (only 
                                               available with `--aligner kma`). 

Generic options
  --multiqc_methods_description      [string]  Custom MultiQC yaml file containing HTML including a methods description.

 !! Hiding 23 params, use the 'validation.showHiddenParams' config value to show them !!
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

If you use BaitCapture in a publication, please consider citing the software using the Zenodo DOI: [10.5281/zenodo.11283946](https://zenodo.org/doi/10.5281/zenodo.11283946).

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) initative, including [nf-core/mag](https://github.com/nf-core/mag) and [nf-core/ampliseq](https://github.com/nf-core/ampliseq), and reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x. PubMed PMID: 32055031.
