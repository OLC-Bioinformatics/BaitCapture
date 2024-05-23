# Alignment reports

BaitCapture provides detailed alignment reports for each sample.
The reports are generated using information gathered from various tools and compiled by custom R scripts in `bin/`.

## `mapstats.tsv`

A table of read alignment statistics against each gene target for each sample, including KMA-specific alignment statistics.

- `sampleid`: Sample ID
- `target`: Target name
- `seqlen`: Target sequence length
- `depth`: Number of bases aligned to the target
- `len_cov`: Length of target covered by aligned bases
- `prop_cov`: Percentage of `seqlen` covered by aligned bases
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
