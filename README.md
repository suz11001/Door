# Door (Domain Organizer)

Door (short for Domain Organizer) is an implementation of a simple proof-of-concept approach for removing conflicting domain architecture within a gene family. Door takes as input a gene family, all domain families present in the given gene family, and a mapping that exist
between domain sequences and gene sequences. Sample data is provided as a template. Door operates in multiple steps (outlined below), briefly, it refines the domain families, extracts the genic region, corrects for domain ordering, and creates a global multiple sequence alignment using MUSCLE. There are two outputs Door generates, Door-S which corrects the domain ordering and pushes the entire sequence (gene + domain) into MUSCLE and Door-A which aligns each domain family and the genic region indepedent via MUSCLE and produces a global MSA by accounting for domain loss and introducing gaps as needed. 


## Dependencies
- Python 2.7.15
  - BioPython
  - Scipy
  - Numpy
- RAxML 8.2.11
- MUSCLE v3


## Overview

Given a gene family multiple sequence alignment, partition the alignment in three approximately equal sized alignments, create maximum likelihood trees, create bootstrap trees, and conduct a statistical analysis for presence/absence of sub-gene/partial gene transfer event.

1. 01_find_domainfamilies_Aln.py - Refines the domain families and the sequences present within the domain family:  

   command: `python 01_find_domainfamilies_Aln.py  path/to/gene/family/fasta path/to/domain/family/fastas/`  
 
2. 02_check_overlap.py - Report how much overlap exists between domain sequences from the same domain family and different domain families. This step is not necessary and can be skipped. 
   
   command: `python 02_check_overlap.py  path/to/01_initial_domains/`

3. 03_find_geneseqs.py - Extract only the genic region of the sequences:

   command: `python 03_find_geneseqs.py  path/to/01_initial_domains/`

4. 04_align.sh - Executes mutiple scripts to generate Door outputs and run muscles. Generates the final global MSA.
   
   command: `bash 04_align.sh `

## Data Set Availability

Simulated data sets and the scripts used to generate such data are available in the sim_data directory.
Real data is available in the fly_data directory.
