---
layout: page
title: Resistance Markers
nav_order: 1
parent: Analysis Modules
---

# Resistance markers

{: .note }
This module can only be run with the `v4` amplicon panel.

This module will identify variants that are known to provide antibiotic resistance to *Plasmodium falciparum* by using the **PseudoCIGAR** string found in `allele_data.txt`. This module identifies mutations found within the genomic coordinates of interest (`resistance_markers_amplicon_v4.txt`), and reports any new indels or SNPs occurred. 

## File Outputs

There are three files that output by this module:

### resmarker_table.txt

This file contains all codons that are found in the ASV as specified by the genomic coordinates in the provided resistance marker table. The file summarizes what the 3-base sequence was in the sample, what was expected, and whether there was a synonomous or non-synomous amino acid change.

|Column|Description|
|:--:|:--:|
|SampleID|The sample being reported|
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

This file provides the same information as the resistance marker table, but in less granular form and joined by haplotype.   

|Column|Description|
|:--:|:--:|
|SampleID|The sample being reported|
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
|GeneID|A numeric identifier for the *P. falciparum* gene and gene position|
|Gene|The name of the gene|
|CodonID|The codon number|
|Position|The position in the reference sequence where the indel or SNP was found|
|Alt|The DNA base found at the position in the ASV|
|Ref|The DNA base found at the position in the reference sequence|
|Reads|The number of reads that support the `Alt` base|

[jekyll-organization]: https://github.com/EPPIcenter
