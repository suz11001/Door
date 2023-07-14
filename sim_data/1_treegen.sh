#!/bin/bash
#SBATCH --job-name=treegen
#SBATCH --mail-user=sumaira.zaman@uconn.edu
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -o treegen_%j.out
#SBATCH -e treegen_%j.err
#SBATCH --partition=general
#SBATCH --qos=general

#generate species tree --> java -jar sagephy-1.0.0.jar HostTreeGen [options] [time interval] [birth rate] [death rate] [out prefix] 
mkdir -p species
for i in {1..100}; do
    seed=$RANDOM
    mkdir -p species/species_$i
    java -jar ./sagephy-1.0.0.jar HostTreeGen -s $seed -min 100 -max 100 -bi 1.0 5.0 2.0 species
    mv species.* species/species_$i
done

#generates 100 gene trees --> java -jar sagephy-1.0.0.jar GuestTreeGen [options] [species tree] [dup rate] [loss rate] [trans rate] [out prefix]
mkdir -p genetrees
for i in {1..100}; do
	seed=$RANDOM
	#echo gene tree: $i,$seed 
	java -jar ./sagephy-1.0.0.jar GuestTreeGen -s $seed -max 250 species/species_$i/species.pruned.tree 0.3 0.3 0 gene
	mkdir -p genetrees/genetree_$i
	mkdir -p genetrees/genetree_$i/prunedGeneTree
	mv gene.pruned.tree genetrees/genetree_$i/prunedGeneTree/
	mv gene.* genetrees/genetree_$i
done
#generate domain tree --> java -jar sagephy-1.0.0.jar DomainTreeGen [options] [species tree] [gene tree directory] [dup rate] [loss rate] [trans rate] [out prefix]
#generates domain trees with varying interspecies trasfer (-is) of 0.2, 0.5, and 0.8

for d in {1,2,3}; do
    mkdir -p domaintree_family_$d
    for i in {1..100}; do
	mkdir -p domaintree_family_$d/domainTree_$i
        mkdir -p domaintree_family_$d/domainTree_$i/domainTree
        mkdir -p domaintree_family_$d/domainTree_$i/mapping
	seed=$RANDOM
        java -jar ./sagephy-1.0.0.jar DomainTreeGen -s $seed -max 250 species/species_$i/species.pruned.tree genetrees/genetree_$i/prunedGeneTree/ 0.3 0.3 0 domainfamily_$d
        mv *.pruned.leafmap domaintree_family_$d/domainTree_$i/mapping
        mv *.pruned.tree domaintree_family_$d/domainTree_$i/domainTree
        mv domainfamily_$d.* domaintree_family_$d/domainTree_$i/
    done
done
