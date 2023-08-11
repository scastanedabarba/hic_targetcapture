# Detection of rare plasmid hosts using a targeted Hi-C approach

###  **Authors:** Salvador Castañeda-Barba, Benjamin J. Ridenhour, Thibault Stalder, Eva M. Top

**********

**Abstract:** Despite the significant role plasmids play in the spread of antibiotic resistance genes, there is limited knowledge of their in situ transfer in natural microbial communities. Therefore, we developed and implemented a new approach to identify rare plasmid hosts by combining Hi-C, a proximity ligation method, with target enrichment for plasmid-specific DNA. Through an experimental design aimed to replicate scenarios in which the transfer frequency of a plasmid is increasingly rare, we first established the detection limit of Hi-C for linking a plasmid to its host in soil. We then implemented the target capture approach and increased the sensitivity of Hi-C 100-fold. Furthermore, we demonstrate that when sufficient reads are present, Hi-C+ can be used to identify plasmid hosts at the genus level. Hi-C+ implementation will facilitate the exploration and determination of the ecological and evolutionary pathways leading plasmids to spread and transfer within microbiomes. 

**********

## Information about this repository:  

This repository contains the code for the analysis carried out in the paper titled "Detection of rare plasmid hosts using a targeted Hi-C approach". 

### **Analysis Descriptions:**
The scripts are broken down into three main analysis; enrichment comparison, math model, and de novo analysis. The prefix for each script indicates the analysis it was used for. For each analysis, there are slurm, bash, and python scripts. These are indicated by the suffix. The slurm script in each analysis was used to run array jobs, therefore it is not described in more detail. The rest of the scripts are further described below. 

#### **Enrichment comparison:**
These scripts were used for comparing detection of plasmid-associated reads between Hi-C and Hi-C libraries. 

#### **Math Model**
These scripts were used for subsampling paired end reads, analysing count data and running math model.

#### **De Novo**
These scripts contain the code for carrying out the de novo analysis. 

#### **Classified plasmid sources table: [classifications.csv](https://github.com/scastanedabarba/plasmid_review_paper/blob/89ec7281a2420897379e42651eb493d6c0f4bee9/classifications.csv)**

