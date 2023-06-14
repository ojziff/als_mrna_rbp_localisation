## Nucleocytoplasmic mRNA redistribution accompanies RNA binding protein mislocalisation in ALS motor neurons and is restored by VCP D2 ATPase inhibition

Oliver Ziff 2023

This page contains scripts to analyse the data and reproduce the figures from the manuscript.

<img src="https://github.com/ojziff/als_mrna_rbp_localisation/figures/graphical abstract.png" height="400">

#### [Manuscript link](https://www.cell.com/neuron/home)

Rmarkdown scripts to reproduce each figure in the manuscript can be found in this repository.

## Data availability

RNA sequencing data have been deposited in the GEO under accession number [GSE214017](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214017). 

Mass Spectrometry data have been deposited in the ProteomeXchange Consortium via the PRIDE partner repository under identifier PXD037107.

Post-mortem RNA sequencing data is accessible at [GSE137810](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137810) and https://collaborators.nygenome.org/. 

Raw fastq files were processed with [nf-core/rnaseq](https://nf-co.re/rnaseq) v3.9 utilising alignment with STAR and read quantification with salmon. 

Differential transcript expression was performed using DESeq2. Differential splicing was analysed with [MAJIQ v2.4](https://majiq.biociphers.org/). Differential protein expression was performed with MaxQuant followed by DEP, which is a wrapper around limma.
