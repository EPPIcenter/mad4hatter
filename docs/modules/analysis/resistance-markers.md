---
layout: page
title: Resistance Markers
nav_order: 1
parent: Analysis Modules
---

# Resistance markers

{: .note }
This module can only be run with the `v4` amplicon panel.

This module is designed to identify variants associated with drug resistance in Plasmodium falciparum by utilizing the **PseudoCIGAR** string extracted from `allele_data.txt`. It works by examining mutations within specific genomic coordinates of interest outlined in (`resistance_markers_amplicon_v4.txt`). Any new indels or single nucleotide polymorphisms (SNPs) that occur within these coordinates are reported. 

In the case that there are multiple loci covering the same drug resistance marker being inspected these will be reported on separate lines. You can identify these as they will have the same GeneID, Gene, and CodonID, but will have different values in the locus column. In most cases we would expect the Codon and amino acid (AA) reported to be the same on both loci, however if you find divergence one option would be to take the one with the highest number of reads. Please reach out to the UCSF team if you have any questions about interpreting these files!  

## File Outputs

There are three files that output by this module:

### resmarker_table.txt

This file contains all codons that are found in the ASV as specified by the genomic coordinates in the provided resistance marker table. The file summarizes what the 3-base sequence was in the sample, what was expected, and whether there was a synonomous or non-synomous amino acid change.

|Column|Description|
|:--:|:--:|
|SampleID|The sample being reported|
|Locus|The locus corresponding to the marker|
|GeneID|A numeric identifier for the *P. falciparum* gene and gene position|
|Gene|The name of the gene|
|CodonID|The codon number|
|RefCodon|The codon in the reference|
|Codon|The codon in the ASV|
|CodonStart|The codon start position|
|CodonRefAlt|Can be 'REF' or 'ALT', depending on whether the ASV codon matches the reference codon ('ALT' if they do not match)|
|RefAA|The amino acid coded by the `RefCodon`|
|AA|The amino acid code by the `Codon`|
|AARefAlt|Can be 'REF' or 'ALT', depending on whether the ASV amino acid matches the reference amino acid ('ALT' if they do not match)|
|Reads|The number of reads that contain this `Codon`|


### resmarker_microhap_table.txt

This file provides the same information as the resistance marker table above, but it is more informative as it strings the mutations from each locus together into haplotypes. 

|Column|Description|
|:--:|:--:|
|SampleID|The sample being reported|
|Locus|The locus corresponding to the marker|
|GeneID|A numeric identifier for the *P. falciparum* gene and gene position|
|Gene|The name of the gene|
|MicrohapIndex|A collapsed verison of the 'CodonID' (see 'resmarker_table.txt'). This will contain all codon IDs that are included in the microhaplotype.|
|RefMicrohap|A collapsed version of the 'Ref' column that reports all reference amino acids in order by `CodonID`|
|Microhaplotype|A collapsed version of the 'Alt' column that reports all ASV amino acids in order by `CodonID`|
|MicrohapRefAlt|Can be 'ALT' or 'REF', depending on whether the ASV microhaplotype matches the reference microhaplotype ('ALT' if they do not match).|
|Reads|Number of reads that have the `Microhaplotype`.|

### resmarker_new_mutations.txt

This file contains DNA mutations that were not within the specified genomic coordinates. The indels and SNPs here are listed in this file to allow end users an opportunity to see other mutations found. They should be interpreted with some caution as they are simply reported and not filtered in any way. 

|Column|Description|
|:--:|:--:|
|SampleID|The sample being reported|
|Locus|The locus corresponding to the marker|
|GeneID|A numeric identifier for the *P. falciparum* gene and gene position|
|Gene|The name of the gene|
|CodonID|The codon number|
|Position|The position in the reference sequence where the indel or SNP was found|
|Alt|The DNA base found at the position in the ASV|
|Ref|The DNA base found at the position in the reference sequence|
|Reads|The number of reads that support the `Alt` base|

[jekyll-organization]: https://github.com/EPPIcenter
