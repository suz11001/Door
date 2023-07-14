# Door (Domain Organizer)

Door (short for Domain Organizer) is an implementation of a simple proof-of-concept approach for removing conflicting domain architecture within a gene family. Door takes as input a gene family, all domain families present in the given gene family, and a mapping that exist between domain sequences and gene sequences. Sample data is provided as a template. Door operates in multiple steps (outlined below), briefly, it refines the domain families, extracts the genic region, corrects for domain ordering, and creates a global multiple sequence alignment using MUSCLE. There are two outputs Door generates, Door-S which corrects the domain ordering and pushes the entire sequence (gene + domain) into MUSCLE and Door-A which aligns each domain family and the genic region indepedent via MUSCLE and produces a global MSA by accounting for domain loss and introducing gaps as needed. 


## Dependencies
- Python 2.7.15
  - BioPython
  - Scipy
  - Numpy
- MUSCLE v.5.0.142


## Overview

The main script (bash.sh) will run from start to finish as long as the input data is provided as in the sample_data directory. 

1. 03_find_geneseqs.py - Extract only the genic region of the sequences and refind the domain family sequences.

   command: `python 03_find_geneseqs.py  path/to/01_initial_domains/gene_family/gene_family.fa path/to/01_initial_domains/gene_family/domain_fams`

2. 04_findDomainPos_optimal.py - This script orders the domain copies in the order in which they occur in the biological sequences. 
   04_ordered_seqs.py - Generates the Door-S output (fasta file) which is later aligned via Muscle to produce a global MSA.
   04_cat_seqs.py - Generates the Door-A output (a global MSA) which contains the onlyGene.aln (alignment) and each domain family alignment appropriately.

   
   command for how to run each of these scripts can be found in bash.sh script.

## Data Set Availability

Scripts used to generate the simulated data are available in the sim_data directory.

