---
layout: default
title: Pipeline Outputs
nav_order: 5
has_children: false
---

# Core Pipeline Outputs

Below is a description of the files that you will find in your output directory. These are files that are part of the core pipeline.

## Amplicon and Sample Coverage

There are two coverage files:

1. sample_coverage.txt
2. amplicon_coverage.txt

These provide sample and amplicon level coverage statistics for your sequencing data. The following metrics can be included for each (further broken down by amplicon in the amplicon_coverage.txt file):

|Metric|Description|
|:--:|:--:|
|Input|This is the starting number of reads that were found for the sample. In the amplicon_coverage.txt file, the starting number of reads for the amplicon for the sample is reported|
|No Dimers|Illumina Adapter dimers are first removed from the sequencing data. This number will inform you how many reads remain after filtering.|
|OutputDada2|This is the number of denoised sequences. DADA2 is the denoising algorithm that the pipeline uses. In this step, reads that do not meet quality thresholds will be filtered out, reducing the number of sequences output from the module|
|OutputPostprocessing|This is the number of sequences that remain after filter out reads that did not pass the specified alignment threshold after aligning to the provided reference sequence. At this step, off target sequences will be filtered out of the final dataset.|

## Allele Data

The allele_data.txt file contains all of the amplicon sequencing variants (ASVs) found in your sequencing dataset. There are 6 columns in this file that will be defined below.

### ASV Identification

There are 5 columns that identify the ASV reported.

|Column|Description|
|:--:|:--:|
|sampleID|The reported sample|
|locus|The reported locus for the sample|
|asv|The denoised ASV|
|reads|The number of reads that support this ASV|
|allele|A unique identifier for the allele that is formed using the locus and an incrementing integer|

### ASV Annotations (PseudoCIGAR)

The `PseudoCIGAR` column provides a pseudocigar string that describes the ASV using *reference* coordinates and keys. The string is a succint representation of all:

* Indels and SNPs that were identified
* Locations that were masked by either user provided masking data, or by built in homopolymer and tandem repeat masking

#### **Mutations (Indels and SNPs)**

##### Indels

The following syntax is used to report insertions:

`{position}I=[ATCG]`

Deletions are reported the same way but with `D=`:

`{position}D=[ACTG]`

In both cases:
 * `position` is where the insertion or deletion occured along the *reference* sequence
 * `[ACTG]` is the base that was inserted in the ASV (does not exist in the reference at that position), or the base that was deleted in the ASV (does exist in the reference at that position, but not in the ASV).

##### SNPs

SNPs use a slightly different syntax:

`{position}[ACTG]`

Where:
 * `position` is where the substitution occured along the *reference* sequence.
 * `[ACTG]` is the new base at the position in the ASV sequence.

#### **Masks**

If you are masking low complexity regions, you may see masking annotation in your PseudoCIGAR sequence. The following syntax is used for masks:

`{start_position}+{mask_length}N`

Where:

* `start_position` is where the mask begins
* `mask_length` is the length of the mask

Any mutations that were identified in this masking region will be superseded by the mask. In other words, a substitution or indel will not be reported in the PseudoCIGAR string if the position is within a masked range. 
