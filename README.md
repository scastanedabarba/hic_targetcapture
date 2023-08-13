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
These scripts were used for comparing detection of plasmid-associated reads between Hi-C and Hi-C libraries. The [alignment script](https://github.com/scastanedabarba/hic_targetcapture/blob/a2e40c8f7c3e2a94853462ef104d2afccf45450f/enrichment_comparison.alignment.sh) was first used to trim the paired end reads, then align them to our reference genomes. The data analysis [bash](https://github.com/scastanedabarba/hic_targetcapture/blob/a2e40c8f7c3e2a94853462ef104d2afccf45450f/enrichment_comparison.data_analysis.sh) and [python](https://github.com/scastanedabarba/hic_targetcapture/blob/a2e40c8f7c3e2a94853462ef104d2afccf45450f/enrichment_comparison.data_analysis.py) scripts were then used to generate count data and compare results between Hi-C and Hi-C+ libraries.

#### **Math Model**
These scripts were used for subsampling paired end reads, analysing count data and running math model. Paired end reads were first subsampled using the [subsample bash script](https://github.com/scastanedabarba/hic_targetcapture/blob/a2e40c8f7c3e2a94853462ef104d2afccf45450f/math_model.subsample.sh). The subsampled paired end reads were aligned using the [align_subsampled](https://github.com/scastanedabarba/hic_targetcapture/blob/a2e40c8f7c3e2a94853462ef104d2afccf45450f/math_model.align_subsampled.sh) script and count data was subsequently generated using the [data_analysis](https://github.com/scastanedabarba/hic_targetcapture/blob/a2e40c8f7c3e2a94853462ef104d2afccf45450f/math_model.data_analysis.py) script. Lastly, the count data was used as input for the math model script. 

#### **De Novo**
These scripts contain the code for carrying out the de novo analysis.The [alignment script](https://github.com/scastanedabarba/hic_targetcapture/blob/a2e40c8f7c3e2a94853462ef104d2afccf45450f/denovo.alignment.sh) was first used to align the paired end reads to our reference genomes. The data analysis [bash](https://github.com/scastanedabarba/hic_targetcapture/blob/a2e40c8f7c3e2a94853462ef104d2afccf45450f/denovo.data_analysis.sh) and [python](https://github.com/scastanedabarba/hic_targetcapture/blob/a2e40c8f7c3e2a94853462ef104d2afccf45450f/denovo.data_analysis.py) scripts were then used to generate counts data and compare results between Hi-C and Hi-C+ libraries.