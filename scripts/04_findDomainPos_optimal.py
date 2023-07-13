from Bio.Align import AlignInfo
import os
import linecache
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse
from Bio import AlignIO
import random
import subprocess
import re
import glob
from collections import OrderedDict


def mapGeneDomain(mapping_file):
    mapping={}
    lines=[line.rstrip('\n') for line in open(mapping_file)]
    for line in lines:
        cols=line.split("\t")
        mapping[cols[0]]=cols[1]
    return(mapping)

def are_overlapping(r, s):
    #r and s are an ordered pair
    print('checking overlap')
    #print(r)
    #print(s)
    return r[1] >= s[0] and s[1] >= r[0]

def align(s1,s2,d):
    try:
        blastp_path="/isg/shared/apps/blast/ncbi-blast-2.7.1+/bin/blastp"
        subprocess.call(blastp_path + " -outfmt 6 -query " + s1 + " -db " + s2 + " -evalue 0.00001 -out tmp.out -num_alignments 1", shell=True,stdout=subprocess.PIPE)
        aln=[line.rstrip('\n') for line in open("tmp.out")]
        aln_list=str("\t".join(aln))
        aln_list_tab=aln_list.split("\t")
        os.remove("tmp.out")
        return(aln_list_tab)
    except subprocess.CalledProcessError as e:
        print(e.output)

def  makeblastdb(domain):
    subprocess.call("/isg/shared/apps/blast/ncbi-blast-2.7.1+/bin/makeblastdb -in " + domain + " -dbtype prot", shell=True,stdout=open(os.devnull, 'wb'))

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


    path = cwd+'sagephy_struct_all_keepOrder/01_initial_domains/'+genesDir+"/*.fa"
    files = glob.glob(path)

    domain_fam_sequences={}
    gene_seqs={}

    id_dict_genes = SeqIO.to_dict(SeqIO.parse(handle_genes, "fasta"))
    
    for gene_rec in id_dict_genes:
        gene_rep={}
        test_list=["dgri", "dmoj", "dvir", "dwil", "dpse","dper","dana","dyak","dere","dsim","dsec","dmel"]
        res = [ele for ele in test_list if(ele in gene_rec)]
        if len(res) == 0:
            continue #continue to the next gene_rec if it is not a fly id
        print('gene record is: ', gene_rec)
        for domain_fam in files: 
            print(domain_fam)
            domain_fam_num = os.path.basename(domain_fam).split(".")[0]
            id_dict_prots = SeqIO.to_dict(SeqIO.parse(domain_fam, "fasta"))
            map_file=cwd+'sagephy_struct_all_keepOrder/01_initial_domains/'+genesDir+"/"+domain_fam_num+"_mapping.txt"

            print(gene_rec)
            mapDG=mapGeneDomain(map_file)
            domains_to_align=[k for k, v in mapDG.items() if v == gene_rec]
            if domains_to_align == []:
                only_gene=str(id_dict_genes[gene_rec].seq).strip().replace("-","")
                gene_seqs[gene_rec]=only_gene

            print('domains to align in this domain family ', domains_to_align)
            gene_seq=str(id_dict_genes[gene_rec].seq).strip().replace("-","")

            #write out the gene sequence to a temp file                                    
            tempGeneFile=open("gene.fa","w")
            tempGeneFile.write(">%s\n%s\n" % (gene_rec,gene_seq))
            tempGeneFile.close()

            domain_positions={} ###key alignment start position and value is the domain family           
            start_stop={}
            for domain in domains_to_align:
                ###write  out the domain to a temp file
                tempDomFile=open("domain.fa","w")
                tempDomFile.write(">%s\n" % (id_dict_prots[domain].id))
                print("domain sequence:", str(id_dict_prots[domain].seq))
                tempDomFile.write("%s\n" % (str(id_dict_prots[domain].seq)))
                tempDomFile.close()
                ###create a blastdb of the domain sequence                                               
                makeblastdb("domain.fa")
                ###align sequence                                                                        
                alignment=align("gene.fa","domain.fa",id_dict_prots[domain].id)
                ###analyze alignment
                if len(alignment)==1:
                    print 'domain ', domain_record, 'could not align to ', gene_rec 
                #elif the alignment is not unique i.e. two alignments start at the same position = int(alignment[6])
                elif int(alignment[6]) in domain_positions.keys():
                    print('alignment is over lapping')
                    #if the new alignment is longer, ONLY then overwrite                                 
                    if int(alignment[7]) > start_stop.get(int(alignment[6])):
                        domain_positions[int(alignment[6])]=domain
                        start_stop[int(alignment[6])]=int(alignment[7])
                    else:
                        print("new alignment is shorter")
                #else alignment does not already exist, put it in the dictionary                           
                else:
                    domain_positions[int(alignment[6])]=domain
                    start_stop[int(alignment[6])]=int(alignment[7])

            print("domain position: ", domain_positions)

            #done iterating through all domains in the domain family that correspond to this gene_record
            if domain_fam_num not in domain_fam_sequences.keys():
                domain_fam_sequences[domain_fam_num]=[]
            else:
                pass
            ###sort the dictionary by ascending keys                                                               
            domain_positions_ordered = OrderedDict(sorted(domain_positions.items(), key=lambda x: x[0]))
            domain_names_ordered=sorted(domain_positions_ordered.items(), key=lambda x:x[0])
            print("domain positions ordered: ", domain_positions_ordered)
            print("domain names ordered: ", domain_names_ordered)
            print("start & stop:", start_stop)
            sorted_start_stops=OrderedDict(sorted(start_stop.items(), key=lambda x: x[0]))
            starts=list(sorted_start_stops.keys())
            ends=list(sorted_start_stops.values())
            
            ###construct the gene sequences again with the domains appended to the end                   
            untouched_gene=gene_seq
            if domain_names_ordered!=[]:
                subgene=""
                for domain in domain_names_ordered:
                    print('working on domain copy:', domain)
                    print('domain start position:', domain[0])
                    #print(type(domain_start.))                                                                     
                    domain_start=domain[0]
                    end=sorted_start_stops.get(domain_start)
                    start=domain_start-1
                    subgene=subgene+untouched_gene[start:end]
                    #appeard the domain name and sequence to the appropriate domain family                           
                domain_fam_sequences[domain_fam_num].append((gene_rec,subgene))

    
    for df in domain_fam_sequences.keys():
        with open(df+".fa", 'w') as out:
            sequences=domain_fam_sequences.get(df)
            for s in sequences:
                out.write(">"+s[0]+"\n"+s[1]+"\n")
        out.close()

    for df_num in domain_fam_sequences.keys():
        muscle_path="/isg/shared/apps/muscle/3.8.31/muscle"
        outf="./"+df_num+".temp.aln"
        in_f=df_num+".fa"
        subprocess.call(muscle_path + " -in " + in_f + " -out " + outf + " -maxiters 12", shell=True)
        aln_dict_dom = SeqIO.to_dict(SeqIO.parse(df_num+".temp.aln", "fasta"))
        aln_dict_gene = SeqIO.to_dict(SeqIO.parse("onlygenes.aln", "fasta"))
        aln_len=AlignIO.read(df_num+".temp.aln", "fasta").get_alignment_length()
        gene_dom_pos={}
        for rec in aln_dict_gene:
            if rec in aln_dict_dom.keys():
                gene_dom_pos[rec]=str(id_dict_genes[rec].seq)
            else:
                gene_dom_pos[rec]="-"*aln_len
    
        with open('domainSeqs_fam_'+df_num+'_'+genesDir+'.aln', 'w') as f1:
            for k in gene_dom_pos.keys():
                records_domain_copies="".join(gene_dom_pos.get(k)) 
                f1.write(">"+k+"\n")
                f1.write(records_domain_copies+"\n")
        f1.close()
