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

def matchIDs():
    global mapping
    mapping={}
    entries="/home/FCAM/szaman/researchMukul/pgtr/biological_data/EntriesWithIDs.txt"
    lines=[line.rstrip('\n') for line in open(entries)]
    for line in lines:
        cols=line.split("\t")
        mapping[cols[0]]=cols[1]

def are_overlapping(r, s):
    #r and s are an ordered pair
    print('checking overlap')
    print(r)
    print(s)
    return r[1] >= s[0] and s[1] >= r[0]

def align(s1,s2,d):
    try:
        print('trying to align')
        blastp_path="/isg/shared/apps/blast/ncbi-blast-2.7.1+/bin/blastp"
        subprocess.call(blastp_path + " -outfmt 6 -query " + s1 + " -db " + s2 + " -evalue 0.00001 -out tmp.out -num_alignments 1", shell=True,stdout=subprocess.PIPE)
        aln=[line.rstrip('\n') for line in open("tmp.out")]
        aln_list=str("\t".join(aln))
        aln_list_tab=aln_list.split("\t")
        #os.remove("tmp.out")
        return(aln_list_tab)
    except subprocess.CalledProcessError as e:
        print(e.output)

def  makeblastdb(domain):
    subprocess.call("/isg/shared/apps/blast/ncbi-blast-2.7.1+/bin/makeblastdb -in " + domain + " -dbtype prot", shell=True)

def prepAln(gene_name,id_dict_genes, id_dict_prots, domain ):
    #align domain sequence to gene to find actual position                                                                                                                
    gene_seq=str(id_dict_genes[gene_name].seq).replace("-","")
    tempGeneFile=open("gene.fa","w")
    tempGeneFile.write(">%s\n%s\n" % (gene_name,gene_seq))
    tempGeneFile.close()

    tempDomFile=open("domain.fa","w")
    tempDomFile.write(">%s\n%s\n" % (id_dict_prots[domain].id, str(id_dict_prots[domain].seq).replace("-","")))
    tempDomFile.close()
    
    ###create a blastdb of the domain sequence                                                                                                                            
    makeblastdb("domain.fa")
    ###align sequence                                                                                                                                                     
    alignment=align("gene.fa","domain.fa",id_dict_prots[domain].id)

    return(alignment)

if __name__ == "__main__":

    matchIDs()


    parser = argparse.ArgumentParser(
    prog='find_domainfamilies.py',
    usage='''python find_domainfamilies.py --genes [path to pretransfer gene files] --out [name of output fasta file]''',
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


    domain_fams={}

    id_dict_genes = SeqIO.to_dict(SeqIO.parse(handle_genes, "fasta"))

    for a in AlignIO.parse(handle_genes, "fasta"):
        for record in a:
            res=[record]
            if len(res) > 0:
                rec_id=record.id
                #gene_dict[record.id]=str(record.seq)
                fbgn_id=rec_id.split("_")[1]
                try:
                    #grep fbgn_id from EntriesWithIDs.txt
                    pfam_id=mapping[fbgn_id]
                    print(pfam_id)
                    #grep pfamID from domain-names-all.txt to find all associated domain -> B4KP41_DROMO-44-89, B4KP41_DROMO-1176-1268, B4KP41_DROMO-1146-1238, B4KP41_DROMO-1082-1202
                    grep_out = subprocess.Popen('grep %s ~/researchMukul/pgtr/biological_data/domain-names-all.txt' % pfam_id, stdout=subprocess.PIPE, shell=True)
                    output = grep_out.stdout.read()
                    domain_sequences= list(output.split("\n"))
                    print('domain seqs', domain_sequences)
                    #grep each domain  from domain_sequences 
                    for ds in domain_sequences[:-1]:
                        print(ds)
                        ds=str(ds)
                        length=None
                        start=None
                        stop=None

                        grep_out2=subprocess.Popen('grep %s ~/researchMukul/pgtr/biological_data/withLength/DomainFamilySeq*' % ds, stdout=subprocess.PIPE, shell=True)
                        output2 = grep_out2.stdout.read()
                        domain_family=output2.split("\n")
                        #print('domain family ', domain_family)

                        if len(domain_family) == 2:
                            result = re.search('DomainFamilySequ(.*):', domain_family[0])
                            df=str(result.group(1))
                            domain_fam="/home/FCAM/szaman/researchMukul/pgtr/biological_data/withLength/DomainFamilySequ"+df
                            id_dict_prots = SeqIO.to_dict(SeqIO.parse(domain_fam, "fasta"))
                            if df in domain_fams.keys():
                                #domain_fams[df].append(ds)
                                print(domain_fams)
                                gene_names=map(lambda x: x[1], domain_fams[df])
                                print(gene_names)

                                alignment=prepAln(rec_id,id_dict_genes, id_dict_prots, ds)
                                if len(alignment)==1:
                                    print 'domain ', ds , 'could not align to ', rec_id
                                else:
                                    start=int(alignment[6])
                                    stop=int(alignment[7])
                                    length=stop-start
                                    
                                    if rec_id in gene_names:
                                        print('code is working')
                                        idx=gene_names.index(rec_id)
                                        #check if they are overlapping
                                        if are_overlapping((domain_fams.get(df)[idx][3],domain_fams.get(df)[idx][4]),(start,stop)):
                                            print('they are overlapping')
                                            existing_len=domain_fams.get(df)[idx][2]
                                            if length > existing_len:
                                                domain_fams[df][idx]=(ds,rec_id,length,start,stop)
                                            else:
                                                pass
                                        else:
                                            domain_fams[df].append((ds,rec_id,length,start,stop))
                                    else:
                                        domain_fams[df].append((ds,rec_id,length,start,stop))
                            else:
                                alignment=prepAln(rec_id,id_dict_genes, id_dict_prots, ds)
                                if len(alignment)==1:
                                    print 'domain ', domain_record, 'could not align to ', gene_rec
                                else:
                                    start=int(alignment[6])
                                    stop=int(alignment[7])
                                    length=stop-start
                                    domain_fams[df]=[(ds,rec_id,length,start,stop)]

                        elif len(domain_family) > 2:
                            print('this domain sequence belongs to multiple domain families - i am confused ', ds)
                        else:
                            print('this domain sequence belongs to no  domain families - i am confused ', ds)
                except:
                    print('cannot find pfam id for ', rec_id)
                    #print('something going wrong')

    print('domain families ', domain_fams)

        
    for k in domain_fams.keys():
        print('domain family ', k)
        with open(str(k)+".fa", 'w') as out, open(str(k)+"_mapping.txt", 'w') as mapF:
            domain_fasta = cwd+"/withLength/DomainFamilySequ"+k
            id_dict_prots = SeqIO.to_dict(SeqIO.parse(domain_fasta, "fasta"))
            for v in domain_fams[k]:
                ###save the domain sequence to a list
                domain_sequence_aln=str(id_dict_prots[v[0]].seq).strip()
                sequence=domain_sequence_aln.replace('-','')
                out.write(">"+v[0]+"\n")
                out.write(sequence+"\n")
                mapF.write(v[0]+"\t"+v[1]+"\n")
        out.close()
        mapF.close()

        if len(domain_fams.get(k))==1:
            print('we have a singleton ', domain_fams.get(k))
        else:
            print('not a singleton')