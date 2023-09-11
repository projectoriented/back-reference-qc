# back-reference-qc

A bioinformatics pipeline that detects unreliable reads with low QV. Quality values are cross-referenced from its own kmerized e.g. Illumina library or your choice of complementary technology for the sample in query. If filtering of low QV are desired, then this can be fine-tuned via the z_filter parameter which is analogous to z score values (peep the z scores projected on standard distribution with a quick google search).

## Getting started
1. Clone the repo
2. Install Snakemake version >= 7.16.0
3. Fill in the config and manifests.
4. Start the analysis!
    * Begin with a dry-run
    ```
    ./runlocal 10 -np
    ```
    * If dry-run looks good, proceed with:
    ```
    ./runlocal 10
    ```

## Output
```commandline
   results/
   ├── plots
   │   └── unique_sample_name
   ├── read_qv
   │   └── unique_sample_name
   └── reads_filtered
       └── unique_sample_name
           ├── fastq
           ├── filtered_out
           │   ├── fasta
           │   └── kraken2
           └── log
```

```commandline
$ cat results/reads_filtered/unique_sample_name/filtered_out/kraken2/summary.tsv.gz
sample  reads_mapped    taxonomy        cell_name
unique_sample_name        54      Mycoplasmataceae (taxid 2092)   m84046_230412_163336_s3.hifi_reads
total   54  N/A     54037280
```
The row total pertains to the sum of all reads mapped to the taxa (e.g. 54), and 54037280 represents the sum of all reads in cells (e.g. m84046_230412_163336_s3.hifi_reads)

```commandline
./results/plots/unique_sample_name
├── kde-after_filter.png
└── kde-before_filter.png
```
The plots showcases the kernel density estimate distribution for qv values before and after filtering. The expectation is that the reads are normally distributed.

```commandline
./results/reads_filtered/unique_sample_name/log/
└── m84046_230412_163336_s3.hifi_reads-extract_undesirable_reads.log

$ cat
comparison_type: self for sample: unique_sample_name
z_filter: -2 for sample: unique_sample_name
QV-99	99.0	m84046_230412_163336_s3.hifi_reads	0.00
QV-0	0.0	m84046_230412_163336_s3.hifi_reads	0.00
QV-low	19.02	m84046_230412_163336_s3.hifi_reads	0.00
```
If comparison_type and z_filter are throwaway headers. We can read the table where lines start with Q. I describe the columns as such:
  - QV thresholds
  - Median value
  - Cell/movie name
  - Proportion of reads that are less than defined z score/filter.

QV thresholds:
  * QV-99 are perfect k-mer matches
  * QV-0 do not match k-mers whatsoever (definitely non-human reads)
  * QV-low are reads beyond the define z_filter score (e.g. -2). 19.02 is the median value of reads that fall beyond a z-score of -2. If we round up to 20, it equates to an error rate of 1 in 100 bases called or 0.01% probability that the base is incorrect.

## To-do
- [ ] Put in CI tests
- [ ] Add conda enviroments
- [ ] Build bare-bone container to get minimal example running

## FAQ
1. What is an example config and manifest?
   1. Answer coming soon.

## Overview
![pipeline vector](https://github.com/projectoriented/back-reference-qc/blob/main/dag.svg)
