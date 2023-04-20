---
layout: page
title: Resistance Markers
nav_order: 1
parent: Analysis Modules
---

# Resistance markers

This module will identify variants that are known to provide antibiotic resistance to *Plasmodium falciparum* by using the previously generated `allele_data.txt` and `Mapping/*mpileup.txt` files, along with a codon table (`codontable.txt`) and resistance marker genomic coordinates (`resistance_markers_amplicon_v4.txt`).

## File Outputs

A table with resistance markers as well as a table with natural haplotypes (resistance markers found in the same amplicon), both including read counts inside brackets.

Whenever multiple variants for a given resistance marker are found in the same sample, they are separated by `_` and the reference variant is presented first followed by the alternate variants. This also applies to the read counts. Cells remain blank if no resistance marker was found.

Nomenclature of resistance markers can be broken down into name of the gene and amino acid position of the marker. For instance, for `dhfr_16`, `dhfr` is the gene and `16` is the amino acid position.

### resmarkers_summary.txt

|SampleName|dhfr_16|dhfr_51|dhfr_59|dhfr_108|dhfr_164|...|
|---|---|---|---|---|---|---|
|sample1|A [350]|I [350]|R [350]||I [314]|...|
|sample2|A [1541]|I_N [49_1492]|R [1541]|N [1460]|I [1460]|...|
|sample3|A [226]|I [226]|R [226]|N [249]|I [249]|...|

### resmarkers_haplotype_summary.txt 

|SampleName|dhfr_16/dhfr_51/dhfr_59|dhfr_108/dhfr_164|mdr1_1034/mdr1_1042|crt_72/crt_73/crt_74/crt_75/crt_76|...|
|---|---|---|---|---|---|
|Sample_A|A/I/R [949]|N/I [798]|S/N [1449]|C/V/M/N/K [859]|...|
|Sample_B|A/I/R [348]|N/I [257]|S/N [399]|C/V/M/N/K [390]|...|
|Sample_C||N/I [94]|S/N [269]|C/V/M/N/K [183]|...|

[jekyll-organization]: https://github.com/EPPIcenter
