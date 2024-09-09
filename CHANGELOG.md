# BaitCapture: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.1.0dev

### `Fixed`

- Fix nf-schema plugin not checking for unique sample IDs when a samplesheet is provided with `--input` (https://github.com/OLC-Bioinformatics/BaitCapture/pull/19)
- Fix null KMA mapping results in `.res` file breaking MERGE_MAPPING_RESULTS process (https://github.com/OLC-Bioinformatics/BaitCapture/pull/18)

## v2.0.0 'Victoria' - [2024-05-24]

### `Added`

- [#e89bca9](https://github.com/OLC-Bioinformatics/BaitCapture/commit/e89bca9) Prepare release `2.0.0`
- [#858569a](https://github.com/OLC-Bioinformatics/BaitCapture/commit/858569a),[#ce7312e](https://github.com/OLC-Bioinformatics/BaitCapture/commit/ce7312e),[#25077c1](https://github.com/OLC-Bioinformatics/BaitCapture/commit/25077c1),[#7115738](https://github.com/OLC-Bioinformatics/BaitCapture/commit/7115738),[#471c1a2](https://github.com/OLC-Bioinformatics/BaitCapture/commit/471c1a2),[#d50805c](https://github.com/OLC-Bioinformatics/BaitCapture/commit/d50805c),[#6505866](https://github.com/OLC-Bioinformatics/BaitCapture/commit/6505866),[#6710a79](https://github.com/OLC-Bioinformatics/BaitCapture/commit/6710a79),[#b0da2cd](https://github.com/OLC-Bioinformatics/BaitCapture/commit/b0da2cd),[#36eea04](https://github.com/OLC-Bioinformatics/BaitCapture/commit/36eea04),[#edc2e24](https://github.com/OLC-Bioinformatics/BaitCapture/commit/edc2e24),[#7a7797f](https://github.com/OLC-Bioinformatics/BaitCapture/commit/7a7797f),[#285111f](https://github.com/OLC-Bioinformatics/BaitCapture/commit/285111f),[#7a3154c](https://github.com/OLC-Bioinformatics/BaitCapture/commit/7a3154c) Added custom alignment and summary statistics report tables which are outputted to `summary/` directory
- [#d914d8b](https://github.com/OLC-Bioinformatics/BaitCapture/commit/d914d8b),[#debc63a](https://github.com/OLC-Bioinformatics/BaitCapture/commit/debc63a) Added options to specify user-defined thresholds for target detection on the command-line
- [#8a9e0f3](https://github.com/OLC-Bioinformatics/BaitCapture/commit/8a9e0f3),[#9e77673](https://github.com/OLC-Bioinformatics/BaitCapture/commit/9e77673),[#7da9fc1](https://github.com/OLC-Bioinformatics/BaitCapture/commit/7da9fc1),[#61dab09](https://github.com/OLC-Bioinformatics/BaitCapture/commit/61dab09) Added option to provide a target metadata file with `--target_metadata` for target cluster detection on the command-line
- [#ddf7d28](https://github.com/OLC-Bioinformatics/BaitCapture/commit/ddf7d28) Added `--report_all` option to report undetected targets in alignment and summary statistics reports
- [#5b76ac9](https://github.com/OLC-Bioinformatics/BaitCapture/commit/5b76ac9),[#5cc1aeb](https://github.com/OLC-Bioinformatics/BaitCapture/commit/5cc1aeb),[#d76186f](https://github.com/OLC-Bioinformatics/BaitCapture/commit/d76186f),[#a1ea22e](https://github.com/OLC-Bioinformatics/BaitCapture/commit/a1ea22e),[#3a9f1ff](https://github.com/OLC-Bioinformatics/BaitCapture/commit/3a9f1ff),[#5a22010](https://github.com/OLC-Bioinformatics/BaitCapture/commit/5a22010),[#4164fdc](https://github.com/OLC-Bioinformatics/BaitCapture/commit/4164fdc),[#bcc349d](https://github.com/OLC-Bioinformatics/BaitCapture/commit/bcc349d) Added [KMA](https://bitbucket.org/genomicepidemiology/kma/) as a read aligner tool, and set as the default aligner
- [#08b0e12](https://github.com/OLC-Bioinformatics/BaitCapture/commit/08b0e12) Added [fastp](https://github.com/OpenGene/fastp) for read trimming and adapter removal
- [#58ace33](https://github.com/OLC-Bioinformatics/BaitCapture/commit/58ace33) Added [fastq-scan](https://github.com/rpetit3/fastq-scan) for obtaining read and base pair counts for raw and preprocessed reads
- [#773995f](https://github.com/OLC-Bioinformatics/BaitCapture/commit/773995f),[#08e5962](https://github.com/OLC-Bioinformatics/BaitCapture/commit/08e5962) Added [Mosdepth](https://github.com/brentp/mosdepth) for read depth summary in MultiQC report
- [#342829a](https://github.com/OLC-Bioinformatics/BaitCapture/commit/342829a) Added configuration file for the OLC cluster

### `Fixed`

- [#b6efdb1](https://github.com/OLC-Bioinformatics/BaitCapture/commit/b6efdb1) Changed `--extension` to `--pattern` to reflect its actual usage on the command-line

### `Dependencies`

- [#663f13f](https://github.com/OLC-Bioinformatics/BaitCapture/commit/663f13f) Upgraded `nf-validation` to `nf-schema=2.0.0`

### `Deprecated`

- [#f1ef4d7](https://github.com/OLC-Bioinformatics/BaitCapture/commit/f1ef4d7),[#de584b3](https://github.com/OLC-Bioinformatics/BaitCapture/commit/de584b3) Removed Trimmomatic as a read aligner tool
- [#24c4cab](https://github.com/OLC-Bioinformatics/BaitCapture/commit/24c4cab) Removed the `PREPROCESS_STATS` module as it's deprecated with the addition of alignment and summary statistics reports
- [#8efda2d](https://github.com/OLC-Bioinformatics/BaitCapture/commit/8efda2d),[#ff449de](https://github.com/OLC-Bioinformatics/BaitCapture/commit/ff449de) Removed unnecessary output dirs

## v1.0.0 - [2023-12-18]

Initial release of BaitCapture, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
