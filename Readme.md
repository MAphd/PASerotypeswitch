This repository contains workflows for 
- Serotype switch detection from core genome alignments coupled with in silico serotyping
- Discovering size and boundaries of serotype O12 islands 
	
	
Howto: 

There are three different analysis in this repository.
 
Dataset1workflow.R covers detection of serotype switches in "DatasetA". 
Assemblylist, tree, and typing data found in /metadata/ DatasetAassemblylist.txt, Dataset1parsnp.tree, and Dataset1meta.Rdata, respectively. 

Dataset2workflow.R covers detection of serotype switches in "DatasetB". 
Assemblylist, tree, and typing data found in /metadata/ Dataset2assemblylist.txt, Dataset2parsnp.tree, and Dataset2meta.Rdata, respectively. 

PA7serotypeislandworkflow.R covers detection and characterization of O12 serotype islands in O12 isolates. 
This analysis was run on genomes from both "dataset1" and "dataset2"

To run analysis, the workflows have been set up so you can simply clone the repository, open the corresponding workflow file, and run all code in said workflow file.

All files required for replication of study have been included in the repository. 

Requirements for replication of study: 

R (4.0.3)
Rstudio (1.3.1093)

Requirements for applying workflow to your own data:

- blastn 2.10.0, build dec 3 2019 
- core genome alignment tool such as parsnp
- phylogeny tool such as fasttree 

R libraries required for analysis:

- readr (1.4.0)
- ggtree (2.4.2)
- treeio (1.14.4)
- tidytree (0.3.3)

R libraries required for plotting:

- ggplot2 (3.3.3)
- ggnewscale (0.4.5)
- viridis (0.6.0)

