import os
import linecache
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse
from Bio import AlignIO
import re

parser = argparse.ArgumentParser(
     prog='04_cat_seqs.py',
     usage='''python 04_cat_seqs.py  --genes [path to pretransfer gene files] --out [name of output fasta file]''',
     description='''This program pulls out specific sequences from a fasta file, given the fasta file and a list of sequences saved in a text file''',
     epilog='''It requires numpy and biopython libraries''')

parser.add_argument('--genes', type=str, help='path to pretransfer gene files', required=True)


args=parser.parse_args()
genesDir=args.genes


gene_seqs={}
domain_seqs={}


domain_families=os.listdir('./domain_fam_alns')

gene_seq="./onlygenes.aln"
gene_dict = SeqIO.to_dict(SeqIO.parse(gene_seq, "fasta"))

for record in gene_dict:
    rec=gene_dict[record].id
    gene_seqs[rec]=str(gene_dict[record].seq)

for df in domain_families: 
    df_num=df.split("_")[2]
    print('working on domain family ', df_num) 


    # set all the files
    handle_domain="./domain_fam_alns/"+df
    ###

    domain_dict={}
    for a in AlignIO.parse(handle_domain, "fasta"):
        for record in a:
            domain_seq=str(record.seq)
            domain_dict[record.id]=domain_seq


    for record in gene_dict:
        rec=gene_dict[record].id
        domain_rec=""
        try:
            domain_rec=domain_dict.get(record)
        except:
            print('no domain record for this gene ', record)

        if rec not in domain_seqs.keys():
            domain_seqs[rec]=domain_rec
        else:
            domain_seqs[rec]=domain_seqs.get(rec)+domain_rec
            
print(domain_seqs)

order_out="./doorA_seqs.aln"

with open (order_out, 'w') as out:
    for a in gene_seqs:
        #print(a)
        #print(gene_seqs.get(a))
        seqs=gene_seqs.get(a)+domain_seqs.get(a)
        out.write(">"+a+"\n")
        out.write(seqs+"\n")
out.close()

