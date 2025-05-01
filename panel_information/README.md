# Panel Information

## Contents 
- [Mad4hatter Pool Information](#pool-information)
    - [Diversity pools](#diversity-pools)
        - [D1.1 (D1)](#d11-d1)
    - [Resistance pools](#resistance-pools)
        - [R1.1](#r11)
        - [R1.2 (R1)](#r12-r1)
        - [R2.1 (R2)](#r21-r2)
- [Other Available Panels](#other-available-panels)
- [Adding new pools](#adding-new-pools)

## Mad4hatter Pool Information

Below is a summary of the pools currently supported by the pipeline. The recommended configuration to maximize information retrieval and sensitivity for low parasitemia samples are two mPCR reactions, one with D1 and R1.2 primers, and one with R2 primers. More detailed information on the protocol can be found [here](https://eppicenter.ucsf.edu/resources).

### Diversity pools

#### D1.1 (D1)

A primer pool targeting 165 high diversity targets and 5 loci targetting the ldh gene in P. falciparum and in 4 non-falciparum Plasmodium species (P. vivax, P. malariae, P. ovale, and P. knowlesi).

### Resistance pools

The resistance module is comprised of two complementary and incompatible primer pools (R1 and R2). 

#### R1.1 

82 targets informing drug and diagnostic resistance, Plasmodium species identification, immune-related targets.

#### R1.2 (R1)

R1.2 is a reduced version of primer pool R1.1, containing 47 targets, designed to reduce primer dimers and increase sensitivity.

#### R2.1 (R2)

31 targets to provide information on drug and diagnostic resistance, Plasmodium species identification, immune-related targets. This complements R1, covering missing codons. 

## Other Available Panels
- [AMPLseq]()
- [4cast]()
- [ama1]()
- [v1]()
- [v2]()

## Adding new pools 

* Create a directory under `panel_information` with the new pool or panel name (e.g. R1.3 or PfPHAST).
* Add the panel information. You could use `panel_information/D1.1/D1.1_amplicon_info.tsv` as an example template. 
* Add the new pool to `conf/panel.config`, including the path to the amplicon info file.
* Generate the targeted reference: 
    * Run the pipeline on an example dataset
    * Use the `genome` flag with a full genome reference, making sure it includes all of your loci (e.g. if you include non-Pf species make sure this is included in your reference). 
    * Set the `pools` flag to only the pool name you set in `conf/panel.config`. Do not set any other pools, even if the example data was generated using other pools. 
    * Copy the file `panel_information/reference.fasta` from the output directory to the new pool directory you set up under `panel_information`. Rename the file to include the pool name (e.g. R1.3_reference.fasta).
    * Add the path to the reference to the `panel.config`. This reference will now be used as default when running with this pool. 
* Resistance markers: 
    * In the output you generated before there will be a file called `panel_information/resmarker_info.tsv` if the pool covers any of the resmarkers in `panel_information/principal_resistance_marker_info_table.tsv`.
    * Check this file includes everything you want. If there are markers missing then add these to `panel_information/principal_resistance_marker_info_table.tsv` and re-run the pipeline to double check your resmarkers have been included. 
    * Add `panel_information/resmarker_info.tsv` from the output directory to the pool directory you created in the repo. This is purely informative as the pipeline will automatically regenerate this file everytime the pipeline is run, unless the path to a resmarker file is set explicitly using the flag `--resmarker_info`.
* If this new pool has superceded an old version (e.g. R1.3 replacing R1.2 as the default for R1) of a pool update the default in `panel.config` under the section `Default pools` and update the warnings in the `check_pools` function in `workflows/validate_inputs.nf`.
* Update documentation 
    * Make sure the help message options are up to date. 
    * Make sure this README and documentation is updated.