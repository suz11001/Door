#!/bin/bash                                                                                                
#SBATCH --job-name=seqgen                                                                                 
#SBATCH --mail-user=sumaira.zaman@uconn.edu                                                                
#SBATCH --mail-type=ALL                                                                                    
#SBATCH -n 1                                                                                               
#SBATCH -N 1                                                                                               
#SBATCH -c 1                                                                                               
#SBATCH --mem=2G                                                                                           
#SBATCH -o seqgen_%j.out                                                                                  
#SBATCH -e seqgen_%j.err                                                                                  
#SBATCH --partition=general                                                                                
#SBATCH --qos=general    

export PATH=/home/FCAM/szaman/researchMukul/pgtr/Seq-Gen-1.3.4/source/:$PATH


cwd="/home/FCAM/szaman/researchMukul/pgtr/biological_data/rearrangement_simulated_multidomain_prot/"
mkdir sequences
for i in {1..100}; do
    mkdir -p sequences/$i
    cd ./sequences/$i
    #mkdir domain_trees
    #mkdir mapping
    #for d in {1..3}; do
	#cp $cwd/domaintree_family_$d/domainTree_$i/domainTree/domainfamily_$d.pruned.tree ./domain_trees
	#cp $cwd/domaintree_family_$d/domainTree_$i/mapping/domainfamily_$d.pruned.leafmap ./mapping
    #done
    for s in {1..3}; do
	mkdir seq_$s
	cd seq_$s
	seed=$RANDOM
	echo "random number for family $i is $seed"
	python3 ~/researchMukul/pgtr/partial.py -m LG -z $seed -dl 100 -gl 150 -s 3 -ap $cwd/domaintree_family_$s/domainTree_$i/domainTree/ $cwd/genetrees/genetree_$i/prunedGeneTree/ $cwd/domaintree_family_$s/domainTree_$i/mapping/
	cd $cwd/sequences/$i/
    done
    #rm -r domain_trees
    #rm -r mapping
    cd $cwd
done
