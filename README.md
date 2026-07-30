# Defining an evidence base for differential transcription factor binding across rDNA haplotypes in mice

Source codes, configurations and R scripts can be found for the project "Defining an evidence base for differential transcription factor binding across rDNA haplotypes in _Mus musculus_" as a part of the module "7BBG1006 Extended Research Project in Applied Bioinformatics 25-26". This work was done under my supervisor Dr. Michelle Holland at the Department of Medical and Molecular Genetics. Pipeline was run primarily on Kings College London's HPC. R version 4.5.3

Author: Yanjing Zhang  
Date: 30 July 2026 

## Table of content 

* **01_TFBS_analysis.Rmd** — identifying TFBSs in rDNA using data on JASPAR
* **02_RShiny** — ENCODE single-cell cell data to show expression profiles of 1,241 TF genes 
    * **2.1_ENCODE_data_extraction.sh**
    * **2.2_ENCODE_data_cleaning.sh**
    * **2.3_Building_RShiny.Rmd**
* **03_ChIPseq_protocol.Rmd**
    * **3.1_Processing_ChIP.sh**
    * **3.2_Haplotypes_ChIP.r**
