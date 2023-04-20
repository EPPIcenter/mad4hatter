---
layout: default
title: Core Pipeline Modules
nav_order: 2
has_children: true
permalink: /mad4hatter/modules/core-pipeline/
---

# Core Pipeline Modules

The mad4hatter amplicon sequencing pipeline is composed of multiple core modules to filter and correct demultiplexed reads to accurately and precisely identify variants in sequencing data.

The following provides a brief synopsis of each module. Each module also has a dedicated page that can be accessed by the table of contents or in the side panel.

## Adapter Removal

This module uses cutadapt to remove illumina adapters and primers from demultiplexed reads. Additionally, reads are quality trimmed and filtered by length, and adapter dimers are removed.

## Quality Control

Coverage statistics table and visualizations are output in this module. Statistics include sample and amplicon specific coverage. *Plasmodium falciparum* specific coverage is calculated during postprocessing after aligning to a reference. 

## Sequence Inference

Real biological sequences are inferred from the reads through a trained error model, and their composition determined on a sample-by-sample basis. Therefore, the composition of the sequences output by this module are not biased by the dataset, and expected to all be corrected by sequencing errors that can impact all samples included in the run. The final output is a counts matrix composed of the determined biological sequences for each sample. 

## Post Processing

The final stage is to filter out target sequences and rearrange the matrix output by the sequence inference module into long form. The table will consist of all samples and amplicons from the panel, and the abundance of alleles found at each locus per sample. 


