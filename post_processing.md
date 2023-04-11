---
layout: page
title: Post-Processing
permalink: /mad4hatter/post-processing/
nav_order: 4
---

# Post-Processing

This module will calculate sample and amplicon coverage from your demultiplexed reads. 

## File Outputs

### alignments{.RDS, .txt}

The sequences from DADA2 are aligned against reference sequences to detect indels and SNPs. This process also outputs a score is used as a filtration step to remove off target sequences. This file contains the original sequence from DADA2, the aligned version (with alignment gaps), the reference sequence that was aligned to (with alignment gaps), and includes the score of the alignment. This file can be useful if you would like to know:

* Where tandem repeat or homopolymer sequencing errors occured (these are masked in the final output). For example, sequences with long tandem repeat regions may have indels that will not be easy to see in the original DADA2 sequences, and are masked to reduce false positives. In this file, long indel gaps in either the hapseq or refseq columns will be much easier to see.
* The distribution of alignment scores. This is useful if you find human reads in your allele table and want to adjust your filtration step. Generally, a score of 60-70 will be a good enough threshold, but you may find that this needs to be adjusted. 

### allele_data{.RDS, .txt}

This file contains the final processed alleles from your demultiplexed sequencing files. The alleles are grouped by locus in each
sample and include relative abundance metrics. You will likely see ‘N’ characters in the allele sequences. These are here to mask low complexity regions that are known to cause sequencing errors (ie.homopolymers and tandem repeats). Without masking, there would be more
unique alleles in the final output that are actually known false positive. The length or number of these regions is a function of the provided reference sequences, and is not impacted by your sequencing data.

[jekyll-organization]: https://github.com/EPPIcenter

