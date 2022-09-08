'''
A pipeline to calculate the Ka, Ks and Ka/Ks
note: This program call on align2yn00.py. 
'''
import os
import sys
from Bio import SeqIO
sys.path.append('.\bin')


def search_orthologs (query_sequence_file, reference_file):
    #qfpath,qfile=os.path.split(query_sequence_file)
    #rfpath,rfile=os.path.split(reference_file)
    database = '.\\temp_data\\blast_result\\blastdb'
    blast_file_name=".\\temp_data\\blast_result\\blast.result"
    
    os.popen("makeblastdb -in %s -dbtype nucl -parse_seqids -out %s" % (reference_file, database))
    
    command = "blastn -num_alignments 1 -db %s -query %s -out %s -outfmt 6" % (database, query_sequence_file, blast_file_name)
    os.system(command)
    return blast_file_name


def blast_table2orthologs_sequence(blast_file_name,query_sequence_file,\
    reference_file):
    query_sequence_record=SeqIO.parse(open(query_sequence_file),'fasta')
    query_list={}
    for i in query_sequence_record:
        query_list[str(i.id).split("|")[0]]=str(i.seq)

    reference_sequence_record=SeqIO.parse(open(reference_file),'fasta')
    reference_list={}
    for j in reference_sequence_record:
        reference_list[str(j.id)]=str(j.seq)
    blast_result=open(blast_file_name).read().split("\n")
    paire_orthologs=[]
    gene_name=[]
    for k in blast_result:
        if k!="":
            gene_name.append(k.split("\t")[0:2])
            paire_orthologs.append(".\\temp_data\\pair_orthologs\\"+k.split("\t")[0].split("|")[0]+"_"+k.split("\t")[1]+".fasta")
            out_file=open(".\\temp_data\\pair_orthologs\\"+k.split("\t")[0].split("|")[0]+"_"+k.split("\t")[1]+".fasta", "w")
            out_file.write(">"+k.split("\t")[0]+"\n"+query_list[k.split("\t")[0].split("|")[0]]+"\n")
            out_file.write(">"+k.split("\t")[1]+"\n"+reference_list[k.split("\t")[1]]+"\n")
            out_file.close()
    return paire_orthologs, gene_name
def align_orthologs(orthologs_file):
    print ".\\bin\\MEGA\\megacc.exe -a .\\bin\\muscle_align_coding.mao -d "+\
    "\\temp_data\\pair_orthologs\\"+os.path.split(orthologs_file)[1]+" -f Fasta -o .\\temp_data\\align_result\\"+os.path.split(orthologs_file)[1]
    os.system(".\\bin\\MEGA\\megacc.exe -a .\\bin\\muscle_align_coding.mao -d "+\
    "\\temp_data\\pair_orthologs\\"+os.path.split(orthologs_file)[1]+" -f Fasta -o .\\temp_data\\align_result\\"+os.path.split(orthologs_file)[1])
    os.popen("C:\\Python27\\python.exe .\\bin\\align2yn00.py -i "+\
    ".\\temp_data\\align_result\\"+os.path.split(orthologs_file)[1]+" -o .\\temp_data\\align_result\\"+os.path.split(orthologs_file)[1]+".phylip")
if __name__ == "__main__":
    len_argv=len(sys.argv)
    for i in range (len_argv):
        if sys.argv[i] == "-q":
            query_sequence_file = sys.argv[i+1]
        if sys.argv[i] == "-r":
            reference_file = sys.argv[i+1]
    print "I blast search between query and reference species and aligned the sequence......"
    blast_file_name=search_orthologs(query_sequence_file, reference_file)
    print "II extracting sequence and transform format into phlip......"
    paire_orthologs, gene_name=blast_table2orthologs_sequence(blast_file_name,query_sequence_file, reference_file)
    print "calculate Ka and Ks"
    RESULT=open("result.summary","w")
    j=0
    for i in paire_orthologs:
        #print i+"  check check"
        Phylip_file=open(".\\temp_data\\KaKs_result\\"+os.path.split(i)[1]+".ctl","w")
        Phylip_file.write("      seqfile = "+".\\temp_data\\align_result\\"+os.path.split(i)[1]+".phylip"+"\n")
        Phylip_file.write("      outfile = "+".\\temp_data\\KaKs_result\\"+os.path.split(i)[1]+".result"+"\n")
        Phylip_file.write("      verbose = 0  "+"\n")
        Phylip_file.write("        icode = 0  "+"\n")
        Phylip_file.write("    weighting = 0  "+"\n")
        Phylip_file.write("   commonf3x4 = 0  "+"\n")
        Phylip_file.close()
        align_orthologs(".\\temp_data\\pair_orthologs\\"+os.path.split(i)[1])
        os.popen(".\\bin\\paml4.9c\\bin\\yn00.exe .\\temp_data\\KaKs_result\\"+os.path.split(i)[1]+".ctl")
        RESULT.write("    ".join(gene_name[j])+"\t")
        print "    ".join(gene_name[j])
        j=j+1
        dn=open("2YN.dN").read().split("\n")[2].split(" ")[-1]
        ds=open("2YN.dS").read().split("\n")[2].split(" ")[-1]
        dn_ds=float(dn)/float(ds)
        YN_t=open("2YN.t").read().split("\n")[2].split(" ")[-1]
        RESULT.write(dn.strip()+" "+ds.strip()+" "+str(dn_ds).strip()+" "+YN_t.strip()+"\n")
        os.remove("2YN.dN")
        os.remove("2YN.dS")
        os.remove("2YN.t")
        os.remove("rst1")
        os.remove("rub")
        os.remove("rst")
    RESULT.close()
    raw=open("result.summary").read().split("\n")
    final_result=open("evolutionary_fv","w")
    temp_id=[]
    for m in raw:
        if m!="":
            if m.split("    ")[0] not in temp_id:
                final_result.write(m+"\n")
        temp_id.append(m.split("    ")[0])
    os.remove("result.summary")
    final_result.close()       
    

    
    
    
    
    
    











 

 
