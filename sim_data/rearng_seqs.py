import os
import linecache
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse
from Bio import AlignIO
import random

parser = argparse.ArgumentParser(
     prog='unorder_seqs.py',
     usage='''python unorder_seqs.py --genes [path to pretransfer gene files] --out [name of output fasta file]''',
     description='''This program pulls out specific sequences from a fasta file, given the fasta file and a list of sequences saved in a text file''',
     epilog='''It requires numpy and biopython libraries''')

parser.add_argument('--genes', type=str, help='path to pretransfer gene files', required=True)
parser.add_argument('--out', type=str, help='name of output fasta file', required=False)

args=parser.parse_args()
genesDir=args.genes
output=args.out

gene_seqs={}
domain_seqs_across={}

probs={10: 0.45, 25: 0.27, 50: 0.24, 75: 0.15, 100: 0.34, 1000: 0.22}
#probs={10: 0.5, 25: 0.30, 50: 0.29, 75: 0.23, 100: 0.35, 1000: 0.26}
#probs={10: 0.5, 25: 0.70, 50: 0.71, 75: 0.77, 100: 0.65, 1000: 0.74}

def reorderDomains(domain_sel_gene,domain_order):
    prb_tandem = 0.63
    new_domain_order = []
    x = 0 if random.random() <= prb_tandem else 1
    #if x== 1 then split tandem and shuffle
    if x == 1: 
        if len(domain_sel_gene) == 3:
            #find the index of domain family w/ multiple domain - just finding the first family
            idx=domain_sel_gene.index([x for x in domain_sel_gene if len(x) > 1][0])
            if domain_order==(2,1,3):
                # if A has duplicated domains 
                if len(domain_sel_gene[0]) > 1: 
                    #pick B1
                    pick=domain_sel_gene[1][0]
                    #insert B1 into the second to last position in A
                    domain_sel_gene[0].insert(-1, pick)
                    #remove B1 from B
                    domain_sel_gene[1].remove(pick)
                elif len(domain_sel_gene[0]) == 1 and len(domain_sel_gene[1]) > 1:
                    pick=domain_sel_gene[1][0]
                    domain_sel_gene[0].insert(0, pick)
                    domain_sel_gene[1].remove(pick)
            elif domain_order==(1,3,2):
                # if B has duplicated domains
                if len(domain_sel_gene[1]) > 1: 
                    #pick C1                                                                                       
                    pick=domain_sel_gene[2][0]
                    #insert C1 into the second to last position in B                                              
                    domain_sel_gene[1].insert(-1, pick)
                    #remove C1 from C                                                                             
                    domain_sel_gene[2].remove(pick)
                elif len(domain_sel_gene[1]) == 1 and len(domain_sel_gene[2]) > 1:
                    pick=domain_sel_gene[2][0]
                    domain_sel_gene[1].insert(0, pick)
                    domain_sel_gene[2].remove(pick)
        elif len(domain_sel_gene) == 2:
            if len(domain_sel_gene[0]) > 1 and len(domain_sel_gene[1]) > 1:
                pick1=domain_sel_gene[0][-1]
                pick2=domain_sel_gene[1][0]
                domain_sel_gene[0][-1]=pick2
                domain_sel_gene[1][0]=pick1
            elif len(domain_sel_gene[0]) > 1 and len(domain_sel_gene[1]) == 1:
                domain_sel_gene[0].insert(-1,domain_sel_gene[1][0])
                domain_sel_gene[1].remove(domain_sel_gene[1][0])
            elif len(domain_sel_gene[0]) == 1 and len(domain_sel_gene[1]) > 1:
                domain_sel_gene[1].insert(-1,domain_sel_gene[0][0])
                domain_sel_gene[0].remove(domain_sel_gene[0][0])
            else:
                #they are all just singular domains
                domain_sel_gene[0][0], domain_sel_gene[1][0] = domain_sel_gene[1][0], domain_sel_gene[0][0]
    # if x==0, not splitting tandems and rearranging as usual
    else:
        if len(domain_sel_gene)== 3:
            if domain_order==(2,1,3):
                domain_sel_gene=[domain_sel_gene[1],domain_sel_gene[0],domain_sel_gene[2]]
            elif domain_order==(1,3,2):
                domain_sel_gene=[domain_sel_gene[0],domain_sel_gene[2],domain_sel_gene[1]]
        elif len(domain_sel_gene) == 2:
            domain_sel_gene=[domain_sel_gene[1],domain_sel_gene[0]]

    new_domain_order=domain_sel_gene

    return new_domain_order

def rearrange(gene_names,domain_seqs_across,domain_order):
    #genes that contain at least 2 domains    
    gene_choices=gene_names
    sel_gene=random.choice(gene_choices)
    domain_sel_gene=domain_seqs_across.get(sel_gene)
    new_domain_order = []
    if d==(2,1,3):
        #check if A or B have multiple domains if so go to some other fxn
        if len(domain_sel_gene[1]) > 1 or len(domain_sel_gene[0]) > 1:
            new_domain_order=reorderDomains(domain_sel_gene,d)
        elif len(domain_sel_gene[1]) == 1 or len(domain_sel_gene[0]) == 1:
            new_domain_order=[domain_sel_gene[1],domain_sel_gene[0],domain_sel_gene[2]]
        elif len(domain_sel_gene[1]) == 0 or len(domain_sel_gene[0]) == 0:
            # remove the empty element domain families
            domain_sel_gene=[x for x in domain_sel_gene if x != []]
            # if any element is duplicated go to tandem duplicate fxn in reorderDomains
            if len([x for x in domain_sel_gene if len(x) > 1])> 0:
                new_domain_order=reorderDomains(domain_sel_gene,domain_order)
            # check that post removal you have at least two domain to shuffle
            elif len(domain_sel_gene)==2:
                new_domain_order=[domain_sel_gene[1],domain_sel_gene[0]]

    elif d==(1,3,2):
        if len(domain_sel_gene[1]) > 1 or len(domain_sel_gene[2]) > 1:
            new_domain_order=reorderDomains(domain_sel_gene,d)
        elif len(domain_sel_gene[1]) == 1 and len(domain_sel_gene[2]) == 1:
            new_domain_order=[domain_sel_gene[0],domain_sel_gene[2],domain_sel_gene[1]]
            domain_seqs_across[sel_gene]=new_domain_order
        elif len(domain_sel_gene[1]) == 0 or len(domain_sel_gene[2]) == 0:
            domain_sel_gene=[x for x in domain_sel_gene if x != []]
            if len([x for x in domain_sel_gene if len(x) > 1]) > 0:
                new_domain_order=reorderDomains(domain_sel_gene,domain_order)
            elif len(domain_sel_gene)==2:
                new_domain_order=[domain_sel_gene[1],domain_sel_gene[0]]

    domain_seqs_across[sel_gene]=new_domain_order
    gene_choices.remove(sel_gene)
    return (domain_seqs_across,gene_choices)

if __name__ == "__main__":

    print('working on gene ', genesDir)

    for i in range(1,4):
        i=str(i)
        gene_f="./sequences/"+genesDir+"/seq_"+i+"/seqs/genes/pre_transfer/gene.pruned.tree"
        gene_dict = SeqIO.to_dict(SeqIO.parse(gene_f, "fasta"))
        domain_f="./sequences/"+genesDir+"/seq_"+i+"/seqs/domains/domainfamily_"+i+".pruned.tree"
        domain_dict = SeqIO.to_dict(SeqIO.parse(domain_f, "fasta"))
        print('the length of the domain fasta files (including domain duplications', len(domain_dict))
        mapping = "./domaintree_family_"+i+"/domainTree_"+genesDir+"/mapping/domainfamily_"+i+".pruned.leafmap"
        map_dict = {}

        import csv
        with open(mapping, 'rb') as csv_file:
            for row in csv.reader(csv_file, delimiter='\t'):
                map_dict[row[0]] = row[1]
                #print(map_dict)

        for record in gene_dict:
            rec=gene_dict[record].id
            domain_names=[]
            domain_seqs=""
            try:
                domain_names=[key for key,value in map_dict.items() if value==rec]
                #print(domain_names)
                for d in domain_names:
                    domain_seqs=domain_seqs+str(domain_dict[d].seq)
                    #print(domain_seqs)
            except:
                print('no domains associated with this gene')                                                  

            if rec not in gene_seqs.keys():
                gene_seqs[rec]=["","",""]
                gene_seqs[rec][int(i)-1]=str(gene_dict[record].seq)
                domain_seqs_across[rec]=[[],[],[]]
            else:
                gene_seqs[rec][int(i)-1]=str(gene_dict[record].seq)

            if domain_seqs!="":
                domain_seqs_across[rec][int(i)-1].append(domain_seqs)
                #domain_seqs_across[rec][int(i)-1]=domain_seqs


    prob_rearng=probs.get(min([ k for k in filter(lambda x : x > len(gene_seqs), probs) ]))

    print('probability of a sequence not being rearranged is :', prob_rearng)

    y = 0 if random.random() <= 0.5 else 1
    d=()
    if y==0:
        d=(2,1,3)
    else:
        d=(1,3,2)

    print('ordering of sequence will be: ', d)

    gene_keys=gene_seqs.keys()
    
    #print(domain_seqs_across)

    for g in gene_keys:
        x = 0 if random.random() <= prob_rearng else 1
        if x==1:
            print('rearranging')
            (domain_seqs_across,gene_keys)=rearrange(gene_keys,domain_seqs_across,d)
            print(domain_seqs_across)
        else:
            print('not rearranging')
        print(len(gene_keys))

unordered="./sequences/"+genesDir+"/rearnged_domains.fa"
with open (unordered, 'w') as out:
    for a in gene_seqs:
        seq=""
        for i in range(3):
            domain_seq=""
            try:
                domain_seq="".join(domain_seqs_across.get(a)[i])
            except:
                print('this gene does not have any domain families at i')
            seq=seq+gene_seqs.get(a)[i]+domain_seq
        out.write(">"+a+"\n")
        out.write(seq+"\n")

out.close()



        
