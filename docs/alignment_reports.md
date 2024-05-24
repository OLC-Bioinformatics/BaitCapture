# Alignment reports

BaitCapture provides detailed alignment reports for each sample.
The reports are generated using information gathered from various tools and compiled by custom R scripts in `bin/`.
To report undetected targets and clusters, use the `--report_all` option on the command-line.

## `mapstats.tsv`

A table of read alignment statistics against each gene target for each sample, including KMA-specific alignment statistics.

- `sampleid`: Sample ID
- `target`: Name of the target sequence
- `seqlen`: Target sequence length
- `depth`: Number of bases aligned to the target
- `len_cov`: Length of target covered by aligned bases
- `prop_cov`: Proportion of `seqlen` covered by aligned bases
- `fold_cov`: Fold-coverage of aligned bases to the target (`depth` / `seqlen`)
- `mapped_reads`: Number of reads aligned to the target

### KMA-specific columns

From [KMA](https://bitbucket.org/genomicepidemiology/kma/)'s documentation:

- `kma_score`: "ConClave score (accumulated alignment score), from all reads that were accepted to match" the target.
- `kma_expected`: "\[E\]xpected Score, if all mapping reads were normally distributed over the entire \[reference\] database" of targets.
- `kma_template_identity`: "\[N\]umber of bases in the consensus sequence that are identical to the template sequence divided by the Template_length. In other words, the percentage of identical nucleotides between template and consensus w.r.t. the template."
- `kma_template_coverage`: "\[P\]ercentage of bases in the template that is covered by the consensus sequence. A Template_Coverage above 100% indicates the presence of more insertions than deletions."
- `kma_query_identity`: "\[N\]umber of bases in the template sequence that are identical to the
consensus sequence divided by the length of the consensus. In other words, the percentage of
identical nucleotides between template and consensus w.r.t. the consensus."
- `kma_query_coverage`: "\[R\]eciprocal values of the Template_Coverage. A Query_Coverage above
100% indicates the presence of more deletions than insertions."
- `kma_depth`: "\[D\]epth of coverage of the template. Commonly referred to as X-coverage, coverage, abundance, etc."
- `kma_q_value`: "\[O\]btained quantile in a Χ<sup>2</sup><sub>1</sub>-distribution, when comparing the obtained Score with the Expected, using a McNemar test."
- `kma_p_value`: "\[O\]btained p-value from the quantile Q_value."

## `sumstats.tsv`

A table of summary statistics for each sample, including the on-target alignment rate, and the number of reads lost from host decontamination and filtering.

- `sampleid`: Sample ID
- `raw_total_reads`: Total number of raw reads in the input sample
- `raw_total_bp`: Total number of base pairs in the input sample
- `fastp_total_reads`: Total number of reads after fastp adapter removal and trimming
- `fastp_total_bp`: Total number of base pairs after fastp adapter removal and trimming
- `decontam_total_reads`: Total number of reads after host decontamination
- `decontam_total_bp`: Total number of base pairs after host decontamination
- `mapped_total_reads`: Total number of reads successfully aligned to target sequences
- `mapped_total_bp`: Total number of base pairs successfully aligned to target sequences
- `percent_reads_lost_fastp`: Percentage of reads lost during fastp adapter removal and trimming
- `percent_reads_lost_decontam`: Percentage of reads lost during host decontamination
- `percent_reads_on_target`: Percentage of reads aligned to target sequences (on-target alignment rate)

## `presence_absence.tsv`

A table of presence-absence calls for each gene target in each sample, based upon user-defined thresholds.

- `sampleid`: Sample ID
- `target`: Name of the target sequence
- `presence_absence`: Presence-absence call for the target in the sample (1 = present, 0 = absent)

## `presence_absence_clusters.tsv`

A table of presence-absence calls for each gene target cluster in each sample, with clusters defined by a target metadata file (e.g. resistance mechanism).
This table is only present if the `--target_metadata` option is used on the command-line.

- `sampleid`: Sample ID
- `metavar_name`: Name of the cluster provided in the target metadata file (e.g. drug class family)
- `metavar_value`: Value of the cluster provided in the target metadata file (e.g. aminoglycoside)
- `presence_absence`: Presence-absence call for the cluster in the sample (1 = present, 0 = absent)