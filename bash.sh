#!/bin/bash

#make sure that sample_data_output and 01_initial_domains directories exist



### step 1

#this is the name/number of your gene family 
genetreeNumber=$1
#path to your muscle (v.5.0.142) executable set to /isg/shared/apps/muscle/5.0.1428/muscle
muscle_exec=$2

## set up your input directory
# mkdir -p sample_data/01_initial_domains/
# mkdir -p sample_data/01_initial_domains/$1
# copy your data into these directory as structured in the sample data

# #set up your output directory
# mkdir -p ./sample_data_output/02_overlap/
mkdir -p ./sample_data_output/03_final_domains/


#### step 2
mkdir -p ./sample_data_output/03_final_domains/$genetreeNumber
cd ./sample_data_output/03_final_domains/$genetreeNumber

python ../../../scripts/03_find_geneseqs.py ../../../sample_data/01_initial_domains/11755/11755.fa /home/FCAM/szaman/Door/sample_data/01_initial_domains/11755/domain_fams/
rm domain.*
rm tempgene.fa 
cd ../../../

#### step 4 
mkdir -p ./sample_data_output/04_alignments/$genetreeNumber
cd ./sample_data_output/04_alignments/$genetreeNumber

in1="../../03_final_domains/$genetreeNumber/onlygenes.fa"
out1="./onlygenes.aln"
$muscle_exec/muscle -align $in1 -output $out1

python ../../../scripts/04_findDomainPos_optimal.py ../../../sample_data/01_initial_domains/$genetreeNumber/$genetreeNumber.fa /home/FCAM/szaman/Door/sample_data/01_initial_domains/$genetreeNumber/domain_fams/ $genetreeNumber $muscle_exec

rm domain.*

mkdir -p domain_fam_alns
mv domainSeqs_fam_*.aln domain_fam_alns/

python ../../../scripts/04_ordered_seqs.py --genes $genetreeNumber > $genetreeNumber.doorS.out 2> $genetreeNumber.doorS.err

in1="./doorS_seqs.fa"
out1="./doorS_seqs.aln"
$muscle_exec/muscle -align $in1 -output $out1

python ../../../scripts/04_cat_seqs.py --genes $genetreeNumber > $genetreeNumber.doorA.out 2> $genetreeNumber.doorA.err

echo "done"

cd ../../../
