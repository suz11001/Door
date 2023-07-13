#!/bin/bash

#make sure that sample_data_output and 01_initial_domains directories exist
mkdir ./sample_data_output/
mkdir ./sample_data_output/01_initial_domains/
mkdir ./sample_data_output/02_overlap/
mkdir ./sample_data_output/03_final_domains/

### step 1
genetreeNumber="11755"
mkdir -p ./sample_data_output/01_initial_domains/$genetreeNumber
cd ./sample_data_output/01_initial_domains/$genetreeNumber
python scripts/01_find_domainfamilies_Aln.py --genes $genetreeNumber > $genetreeNumber.out 2> $genetreeNumber.err

### step 2
mkdir -p ./sample_data_output/02_overlap/$genetreeNumber
cd ./02_overlap/$genetreeNumber
python scripts/02_check_overlap.py --genes $genetreeNumber > $genetreeNumber.out 2> $genetreeNumber.err

rm ./sample_data_output/01_initial_domains/$genetreeNumber/domain.fa*
rm ./sample_data_output/01_initial_domains/$genetreeNumber/gene.fa 

#### step 3
mkdir -p ./sample_data_output/03_final_domains/$genetreeNumber
cd ./sample_data_output/03_final_domains/$genetreeNumber

python scripts/03_find_geneseqs.py --genes $genetreeNumber > $genetreeNumber.out 2> $genetreeNumber.err

#### step 4 
mkdir -p ./04_alignments/$line
cd ./04_alignments/$line

in1="./sample_data_output/03_final_domains/$genetreeNumber/onlygenes.fa"
out1="./onlygenes.aln"
muscle -in $in1 -maxiters 12 -out $out1

python ../../scripts/04_findDomainPos_optimal.py --genes $line > $line.out 2> $line.err
rm domain.*

mkdir -p domain_fam_alns
mv domainSeqs_fam_*.aln domain_fam_alns/

python ../../scripts/04_ordered_seqs.py --genes $line > $line.ord.out 2> $line.ord.err

in1="./ordered_seqs.fa"
out1="./ordered_seqs.aln"
muscle -in $in1 -maxiters 12 -out $out1 

python ../../scripts/04_cat_seqs.py --genes $line > $line.cat.out 2> $line.cat.err
