# Door (Domain Organizer)

Door (short for Domain Organizer) implements two simple proof-of-concept approachs for improving gene family sequence alignments by correcting for out-of-order protein domains. Door takes as input a gene family, all domain families present in the given gene family, and a mapping of the domain sequences to their corresponding gene sequences. It outputs a global multiple sequence alignment for the input gene family using the two algorithms, Door-S and Door-A, described in the paper cited below. Both algorithms reorder the domains and non-domain regions of the gene sequences to allow for better alignment of homologous sequence regions. Door-S reorders the domain and non-domain regions in the input gene sequences and uses MUSCLE to compute the final global alignment, while Door-A organises and aligns the sequences belonging to each domain family independently (using MUSCLE)and then produces a global alignment by concatenating the different domain family and non-domain region alignments.  

Door can be cited as follows:

Reducing the impact of domain rearrangement on sequence alignment and phylogeny reconstruction<br>
Sumaira Zaman, Mukul S. Bansal<br>
Under review


## Dependencies
- Python 2.7.15
  - BioPython
  - Scipy
  - Numpy
- MUSCLE v.5.0.142


## Overview

The main script (bash.sh) will run from start to finish as long as the input data is provided as in the sample_data directory. 

1. 03_find_geneseqs.py - Extract only the genic region of the sequences and refine the domain family sequences.

   command: `python 03_find_geneseqs.py  path/to/01_initial_domains/gene_family/gene_family.fa path/to/01_initial_domains/gene_family/domain_fams`

2. 04_findDomainPos_optimal.py - This script orders the domain copies in the order in which they occur in the biological sequences. 
   04_ordered_seqs.py - Generates the Door-S output (fasta file) which is later aligned via Muscle to produce a global MSA.
   04_cat_seqs.py - Generates the Door-A output (a global MSA) which contains the onlyGene.aln (alignment) and each domain family alignment appropriately.

   
   command for how to run each of these scripts can be found in bash.sh script.

## Data Set Availability

Scripts used to generate the simulated data are available in the sim_data directory.

