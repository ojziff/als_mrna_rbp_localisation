## Nucleocytoplasmic mRNA redistribution accompanies RNA binding protein mislocalisation in ALS motor neurons and is restored by VCP D2 ATPase inhibition

Oliver Ziff 2023

This repository contains an R Markdown file that includes the scripts necessary to generate each figure in the manuscript. The Rmd file is organised into sections corresponding to the figures in the manuscript. Schematics used in the figures were created using Biorender. Please refer to the instructions provided within the Rmd file to reproduce the figures. We have included a separate script called "deseq_objects.R" that generates the DESeq2 and DEP objects required for the analysis. Additionally, the script entitled "R_workspace.R" loads all the necessary packages and functions required for the figures. Please run these scripts before executing the R Markdown file to ensure all dependencies are properly set up.

#### [Manuscript link](https://www.sciencedirect.com/science/article/pii/S0896627323004786)

## Data availability

RNA sequencing data have been deposited in the GEO under accession number [GSE214017](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214017). 

Mass Spectrometry data have been deposited in the ProteomeXchange Consortium via the PRIDE partner repository under identifier PXD037107.

Post-mortem RNA sequencing data is accessible at [GSE137810](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137810) and https://collaborators.nygenome.org/. 

Raw fastq files were processed with [nf-core/rnaseq](https://nf-co.re/rnaseq) v3.9 utilising alignment with STAR and read quantification with salmon. 

Differential transcript expression was performed using DESeq2. Differential splicing was analysed with [MAJIQ v2.4](https://majiq.biociphers.org/). Differential protein expression was performed with MaxQuant followed by DEP, which is a wrapper around limma.

<img src="[https://github.com/ojziff/als_mrna_rbp_localisation/figures/graphical abstract.png](https://github.com/ojziff/als_mrna_rbp_localisation/blob/main/figures/graphical%20abstract.png)" height="400">

