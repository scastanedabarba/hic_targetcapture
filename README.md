# Detection of rare plasmid hosts using a targeted Hi-C approach

###  **Authors:** Salvador Castañeda-Barba, Benjamin J. Ridenhour, Thibault Stalder, Eva M. Top

**********

**Abstract:** Despite the significant role plasmids play in the spread of antibiotic resistance genes, there is limited knowledge of their in situ transfer in natural microbial communities. Therefore, we developed and implemented a new approach to identify rare plasmid hosts by combining Hi-C, a proximity ligation method, with target enrichment for plasmid-specific DNA. Through an experimental design aimed to replicate scenarios in which the transfer frequency of a plasmid is increasingly rare, we first established the detection limit of Hi-C for linking a plasmid to its host in soil. We then implemented the target capture approach and increased the sensitivity of Hi-C 100-fold. Furthermore, we demonstrate that when sufficient reads are present, Hi-C+ can be used to identify plasmid hosts at the genus level. Hi-C+ implementation will facilitate the exploration and determination of the ecological and evolutionary pathways leading plasmids to spread and transfer within microbiomes. 

**********

## Information about this repository:  

This repository contains the code for the analysis carried out in the paper titled "Detection of rare plasmid hosts using a targeted Hi-C approach".  

### **File Description:**
#### **Data obtained from PLSDB: [plsdb.tsv](https://github.com/scastanedabarba/plasmid_review_paper/blob/d8bad770b352d069dc2ec795595e12e17d53c6ce/plsdb.tsv) and [plsdb.abr](https://github.com/scastanedabarba/plasmid_review_paper/blob/d8bad770b352d069dc2ec795595e12e17d53c6ce/plsdb.abr)**
The files named plsdb.tsv and plsdb.abr were downloaded from PLSDB, version 2021_06_23_v2. The file plsdb.tsv contains the meta-data for all the plasmids while the file plsdb.abr contains the resistance gene annotations. 
#### **Code for meta-analysis: [plasmid_meta.py](https://github.com/scastanedabarba/plasmid_review_paper/blob/89ec7281a2420897379e42651eb493d6c0f4bee9/plasmid_meta.py)**
The file named plasmid_meta.py contains the script for carrying out the meta-analysis in Figure 3. The columns 'Host_BIOSAMPLE' and 'IsolationSource_BIOSAMPLE' within plsdb.tsv contain metadata related to the source from which the plasmid was isolated. This information was used to determine whether each plasmid originated from Human, Animal, or Environmental habitats. 
Lines 22 to 289 contain the dictionaries and code used for classification of plasmid sources into habitats. In lines 292 to 390, the ARG annotations are merged with the classification information and data is prepared for plotting.
#### **Classified plasmid sources table: [classifications.csv](https://github.com/scastanedabarba/plasmid_review_paper/blob/89ec7281a2420897379e42651eb493d6c0f4bee9/classifications.csv)**
The table classifications.csv contains the habitat to which each plasmid was assigned. The ACC_NUCCORE, Host_BIOSAMPLE, and IsolationSource_BIOSAMPLE columns from plsdb.tsv were also retained.
