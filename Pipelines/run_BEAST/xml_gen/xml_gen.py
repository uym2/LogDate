#! /usr/bin/env python

from sys import argv
from sequence_lib import read_fasta
import re

template_file = argv[1]
sampling_time = argv[2] 
sequence_file = argv[3]
tree_file = argv[4]
log_file = argv[5]
output_file = argv[6]


def write_tree(tree_file,fout):
    with open(tree_file) as fin:
        tree = fin.readline();
        fout.write("\t\t"+tree)

def write_alignment(sequence_file,fout):
        taxon_names, seq_aln = read_fasta(sequence_file)
        for (taxon,seq) in zip(taxon_names,seq_aln):
            if taxon != "out":
                fout.write("\t\t<sequence>\n")
                fout.write("\t\t\t<taxon idref=\"" + taxon + "\"/>\n")
                fout.write("\t\t\t"+seq + "\n")
                fout.write("\t\t</sequence>\n")

def write_sampling_time(infile,fout,write_date=True):
    with open(infile,'r') as fin:
        fin.readline()
        for line in fin:
            taxon, time = line.strip().split()
            fout.write("\t\t<taxon id=\"" + taxon + "\">\n")
            if time and write_date:
                fout.write("\t\t\t<date value=\"" + time + "\" direction=\"forwards\" units=\"years\"/>\n")
            fout.write("\t\t</taxon>\n")

with open(template_file,'r') as fin:
    with open(output_file,'w') as fout:
        while 1:
            line = fin.readline()
            fout.write(line)
            if line.strip()[:8] == "<taxa id":
                taxa_start = line
                break
        
        has_date = True if (sampling_time[-1] == '1') else False        
        write_sampling_time(sampling_time[:-1],fout,write_date=has_date)
        
        while 1:
            line = fin.readline()
            if line.strip()[:6] == "</taxa":
                fout.write(line)
                break
        
        while 1:
            line = fin.readline()
            if not line.strip()[:2] == "<!":
                fout.write(line)
            if line.strip()[:10] == "<alignment":
                break
        
        write_alignment(sequence_file,fout)
        
        while 1:
            line = fin.readline()
            if line.strip()[:11] == "</alignment":
                fout.write(line)
                break
        
        while 1:
            line = fin.readline()
            if not line.strip()[:2] == "<!":
                fout.write(line)
            if line.strip()[:7] == "<newick":
                break
        write_tree(tree_file,fout)
        while 1:
            line = fin.readline()
            if line.strip()[:8] == "</newick":
                fout.write(line)
                break

        for line in fin:
            if line.strip()[:17] == "<log id=\"fileLog\"":
                fout.write(re.sub("fileName=.*txt","fileName=\""+log_file+".log.txt",line))
            elif line.strip()[:8] == "<logTree":
                fout.write(re.sub("fileName=.*txt","fileName=\""+log_file+".tree.txt",line))
            else:    
                fout.write(line)    
