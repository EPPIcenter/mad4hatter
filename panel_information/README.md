# Panel Information

## Contents 
- [Pool Information](#pool-information)
    - [MADHatTeR Pools](#madhatter-pools)
        - [Diversity pools](#diversity-pools)
            - [D1.1 (D1)](#d11-d1)
        - [Resistance pools](#resistance-pools)
            - [R1.1](#r11)
            - [R1.2 (R1)](#r12-r1)
            - [R2.1 (R2)](#r21-r2)
    - [PfPHAST Pools](#pfphast-pools)
        - [M1.1 (M1)](#m11-m1)
        - [M2.1 (M2)](#m21-m2)
        - [M1.addon](#m1addon)
    - [Other Pools](#other-pools)
        - [ama1](#ama1)
        - [4cast](#4cast)
        - [v1](#v1)
        - [v2](#v2)
- [Adding new pools](#adding-new-pools)

## Pool Information

Below is a summary of the pools currently supported by the pipeline. 

### MADHatTeR Pools 
The recommended configuration to maximize information retrieval and sensitivity for low parasitemia samples are two mPCR reactions, one with D1 and R1.2 primers, and one with R2 primers. More detailed information on the protocol can be found [here](https://eppicenter.ucsf.edu/resources).

#### Diversity pools

##### D1.1 (D1)

A primer pool targeting 165 high diversity targets and 5 loci targetting the ldh gene in P. falciparum and in 4 non-falciparum Plasmodium species (P. vivax, P. malariae, P. ovale, and P. knowlesi). This pool was previously named 1A.

#### Resistance pools

The resistance module is comprised of two complementary and incompatible primer pools (R1 and R2). 

##### R1.1 

82 targets informing drug and diagnostic resistance, Plasmodium species identification, immune-related targets. This pool was previously called 1B.

##### R1.2 (R1)

R1.2 is a reduced version of primer pool R1.1, containing 47 targets, designed to reduce primer dimers and increase sensitivity. This pool was previously called 5.

##### R2.1 (R2)

31 targets to provide information on drug and diagnostic resistance, Plasmodium species identification, immune-related targets. This complements R1, covering missing codons. This pool was previously called 2.

### PfPHAST Pools 

#### M1.1 (M1)

41 targets for drug and diagnostic resistance, Plasmodium spp. identification, csp, and 10 diversity targets. "Minimal" set for prioritary markers and TES classification

#### M2.1 (M2) 

15 targets for drug and diagnostic resistance, Plasmodium spp. identification, csp, and 10 diversity targets. "Minimal" set for prioritary markers and TES classification. Complements M1.1 for missing codons.

#### M1.addon1

4 targets for P. vivax, P. malariae, P. ovale and P. knowlesii mitochondrial cytb targets for increased sensitivity in non-Pf detection (recommended to add to M1.1).

### Other Pools 

#### ama1 
A single-target panel for analyzing diversity in the apical membrane antigen 1 (AMA1) gene, which is important for invasion of red blood cells by Plasmodium parasites. This panel is useful for studying AMA1 genetic diversity and population structure. For more information, see [Miller et al. (2017)](https://pubmed.ncbi.nlm.nih.gov/29246158/).

#### 4cast 
4CAST is a small multiplex of four highly polymorphic antigenic loci: CSP, AMA1, SERA2 and TRAP. For more information, see [LaVerriere et al. (2021)](https://doi.org/10.1111/1755-0998.13622).

#### v1 
Version 1 of the MADHatTeR panel 

#### v2 
Version 2 of the MADHatTeR panel

#### SpotMalaria
A 136 target panel consisting of three pools: spotmalaria_grc1, spotmalaria_grc2, spotmalaria_spec. Targets cover a 100 SNP barcode, 2 markers for species identification, and key drug resistance markers. For more information, see [Jacob CG et al. (2021)](https://doi.org/10.7554/eLife.62997).

## Adding new pools 

* Create a directory under `panel_information` with the new pool name (e.g. R1.3).
* Add the panel information. You could use `panel_information/D1.1/D1.1_amplicon_info.tsv` as an example template. 
* Add the new pool to `conf/panel.config`, including the path to the amplicon info file.
* Generate the targeted reference: 
    * Run the pipeline on an example dataset
    * Use the `genome` flag with a full genome reference, making sure it includes all of your loci (e.g. if you include non-Pf species make sure this is included in your reference). 
    * Set the `pools` flag to only the pool name you set in `conf/panel.config`. Do not set any other pools, even if the example data was generated using other pools. 
    * Copy the file `panel_information/reference.fasta` from the output directory to the new pool directory you set up under `panel_information`. Rename the file to include the pool name (e.g. R1.3_reference.fasta). Double check the number of targets in the reference matches the number in the panel information.
    * Add the path to the reference to the `panel.config`. This reference will now be used as default when running with this pool. 
* Resistance markers: 
    * In the output you generated before there will be a file called `panel_information/resmarker_info.tsv` if the pool covers any of the resmarkers in `panel_information/principal_resistance_marker_info_table.tsv`.
    * Check this file includes everything you want. If there are markers missing then add these to `panel_information/principal_resistance_marker_info_table.tsv` and re-run the pipeline to double check your resmarkers have been included. 
    * Add `panel_information/resmarker_info.tsv` from the output directory to the pool directory you created in the repo. This is purely informative as the pipeline will automatically regenerate this file everytime the pipeline is run, unless the 
* If this new pool has superceded an old version of a pool update the default in `panel.config` and update the warnings in the `check_pools` function in `workflows/validate_inputs.nf`.
* Update documentation 
    * Make sure the help message options are up to date. 
    * Make sure the README and documentation is updated.