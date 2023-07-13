import itertools
import os
import ast
import linecache
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse
from Bio import AlignIO
import random
import subprocess
import re

def are_overlapping(r, s):
    #r and s are an ordered pair                                                                                                                                                                    
    print('checking overlap')
    return r[1] >= s[0] and s[1] >= r[0]


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
    prog='find_domainfamilies.py',
    usage='''python find_geneseqs.py --genes [path to pretransfer gene files] --out [name of output fasta file]''',
    description='''This program pulls out specific sequences from a fasta file, given the fasta file and a list of sequences saved in a text file''',
    epilog='''It requires numpy and biopython libraries''')

    parser.add_argument('--genes', type=str, help='path to pretransfer gene files', required=True)
    parser.add_argument('--out', type=str, help='name of output fasta file', required=False)

    args=parser.parse_args()
    genesDir=args.genes
    output=args.out

    cwd="/home/FCAM/szaman/researchMukul/pgtr/biological_data/"
    print('working on gene family ', genesDir)

    handle_genes=cwd+"/aligns-pep/"+genesDir+".pep.align"
    if os.path.isfile(handle_genes):
        pass
    else:
        handle_genes=cwd+"/aligns-pep-filtered/"+genesDir+".pep.align"

    id_dict_genes = SeqIO.to_dict(SeqIO.parse(handle_genes, "fasta"))

    check={}

    output_file="/home/FCAM/szaman/researchMukul/pgtr/biological_data/sagephy_struct_all/01_initial_domains/"+genesDir+"/"+genesDir+".out"
    patter="domain families"
    grep_out = subprocess.Popen('grep -w families %s ' % output_file, stdout=subprocess.PIPE, shell=True)
    output = grep_out.stdout.read()
    #print(output)
    #domain_families= output.strip().split(",")[1:]
    domain_families= ast.literal_eval(output.strip()[21:-1])
    print(domain_families)

    for gene_rec in id_dict_genes:
        for k,v in domain_families.items():
            for values in v:
                fb_gn=values[1]
                if fb_gn == gene_rec:
                    if fb_gn in check.keys():
                        check[fb_gn].append((values[0],k,int(values[3]),int(values[4])))
                    else:
                        check[fb_gn]=[(values[0],k,int(values[3]),int(values[4]))]
                else:
                    pass
    print(check)


    for gene in check.keys():
        vals=check.get(gene)
        combs=list(itertools.combinations(vals, 2))
        for c in combs:
            print(c)
            x1=c[0][2:]
            x2=c[1][2:]
            #are_overlapping(x1,x2)
            if (are_overlapping(x1,x2)):
                print('oh no domains are overlaping')
                print(gene,c)
