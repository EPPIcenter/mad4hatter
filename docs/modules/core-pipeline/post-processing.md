---
layout: page
title: Post-Processing
nav_order: 3
parent: Core Pipeline Modules
toc: true
---

# Post-Processing

This module rearranges the counts matrix produced from [sequence inference]({% link docs/modules/core-pipeline/seq-inference.md %}) and filters out sequences that do not map to *Plasmodium falciparum*. Users have the option to additionally mask regions of low complexity that are known to cause sequencing error, and can cause spurious output reducing precision. 

### Masking Low Complexity Regions

There are a number of steps in the postprocessing module (`DADA2_POSTPROC`) to reduce false positives in your final allele table. In particular, tandem repeats and hompolymer regions are masked due to known problems that they cause in Illumina instruments. Masked regions will appear as `N`'s in your sequences such as the example below.

Below are a table of parameters to control how much masking of low complexity regions should occur. These can be modified in the `nextflow.config` file.

|Parameter|Description|
|---|---|
|add_mask|Whether to add a mask or not to the final sequences (default `true`)|
|trf_min_score|Used by Tandem Repeat Finder. This will control the alignment score required to call a sequence a tandem repeat (default `25`)|
|trf_max_period|Used by Tandem Repeat Finder. This will limit the range of the pattern size of a tandem repeat (default `3`)|
|homopolymer_threshold|The number of repeating bases to qualify a sequence as a homopolymer region (default `5`)|

### Examples 

Here is a sequence from DADA:

{: .highlight }
TATATATATATATATATATATATATATATATATATATATATATATGTATGTATGTTGATTAATTTGTTTATATATTTATATTTATTTCTTATGACCTTTTTAGGAACGACACCGAAGCTTTAATTTACAATTTTTTGCTATATCCATGTTAGATGCCTGTTCAGTCATTTTGGCCTTCATAGGTCT

And here is it's masked counterpart using the pipeline defaults. Notice how the homopolymers and tandem repeats are masked by `N`s.

{: .highlight }
**NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN**GTATGTATGTTGATTAATTTGTTTATATATTTATATTTATTTCTTATGACCTTTTTAGGAACGACACCGAAGCTTTAATTTACA**NNNNNNNN**CTATATCCATGTTAGATGCCTGTTCAGTCATTTTGGCCTTCATAGGTCT


{: .note }
The masking is accomplished with [Tandem Repeat Finder](https://github.com/Benson-Genomics-Lab/TRF#trf-definitions). Please refer to their documentation for additional information.



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

