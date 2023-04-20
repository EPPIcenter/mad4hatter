---
layout: page
title: Sequence Inference
permalink: /mad4hatter/modules/core-pipeline/seq-inference/
nav_order: 2
parent: Core Pipeline Modules
---

# Sequence Inference

The pipeline uses [DADA2](https://benjjneb.github.io/dada2/index.html) to infer real alleles from sample data. After adapater removal and initial quality filtering from the [Quality Control]({{ site.baseurl }}{% link docs/modules/core-pipeline/quality-control.md %}) module, reads are further filtered and trimmed based on the expected number of errors, and to remove any ambiguous bases. Remaining sequences are used to train the error model, which is used by the core denoising algorithm to generate real sequences with sequencing errors removed. The final output is a counts table with the number of times each sequence was seen in the sample. This sample table is later used by the Post-Processing module to optionally mask low complexity regions, and is rearranged for users to quickly inspect detected alleles at each locus per sample. 

## File Outputs

### DADA2.RData and seqtab.nochim.RDS

These files contain the environment from the DADA2 module and the final sequence table, respectively. seqtab.nochim.RDS is also contained within DADA2.RData, and is primarily used as the file that is passed on to the postprocessing module. The contents of these files are useful if you would like to interrogate aspects of the DADA2 analysis, such as how sequences were clustered under the error model.


[jekyll-organization]: https://github.com/EPPIcenter
