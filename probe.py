import os
import pandas as pd
from Bio import Entrez

# 0 Prepare for further analysis
# make a dictionary for all sequences with accession number as key and systematic name as value
df=pd.read_csv('targets.csv')
sys_name = df.set_index(['Accession'])['Systematic_name'].to_dict()

# make a dictionary for full genome with accession number as key and sequence as value
genome = {}
acc = ''
full_genome_file = open("full_genome.txt", 'r')

for line in full_genome_file:
    if line.startswith('>'):   
        acc=line[1:].split()[0].split(":")[0]   
        genome[acc]=line
    else:
        genome[acc]+=line

full_genome_file.close() 

# 1 Download sequences and change names
# Download sequences
Entrez.email = "example@ed.ac.uk"
sequence_file = open("all_sequences.fasta", "w")
for acc_id in list(sys_name.keys()):
	if acc_id not in list(genome.keys()):
		print("Download sequence for "+acc_id)
		handle = Entrez.efetch(db="nucleotide", id=acc_id, rettype="fasta", retmode="text")
		record = handle.read()
		sequence_file.write(record)
	else: 
		record=genome[acc_id]
		sequence_file.write(record)

sequence_file.close()


print("The number of sequences is:")
os.system('grep -o ">" all_sequences.fasta | wc -l')


# Change the title of downloaded sequences into systematic name
all_sequence_sys_file = open("all_sequences_sys.fasta", "w")
all_sequences_file = open("all_sequences.fasta", 'r')

for line in all_sequences_file:
    if line.startswith('>'):   
        acc=line[1:].split()[0].split(":")[0]
        line='>'+sys_name[acc]+'\n' 
        all_sequence_sys_file.write(line)
    else:
        all_sequence_sys_file.write(line)

all_sequences_file.close()
all_sequence_sys_file.close()

# 2 Make 64-mer probes
all_sequence_sys_file = open("all_sequences_sys.fasta", "r")
raw_probes_file = open("raw_probes.fasta", 'w')

num=0 # Count sequences being processed
for line in all_sequence_sys_file:
    if line.startswith('>'):
    	name=line.rstrip()[1:]
    	counter=0 #Count sequence lines for each sequence
    	probe_counter=0 # Count probes for each sequence
    	num+=1
    	print("Currently proccessing sequence "+str(num)+" ("+name+")")
    elif len(line.rstrip()) >= 64: # Only keep the sequence line with more than 64 bp
    	counter+=1
    	if (counter & 1) != 0: # Only keep the lines in odd row
    		probe_counter+=1
    		raw_probes_file.write(">"+name+"_"+str(probe_counter)+"\n") # Print title for each probe
    		raw_probes_file.write(line[:64]+"\n") # Print 64-mer probe sequence

all_sequence_sys_file.close()
raw_probes_file.close()

# 3 Online blast against human sequences and remove the probes that have matches
# BLAST against human
bla_homo='blastn -db nt -query raw_probes.fasta -remote -entrez_query "Homo sapiens [Organism]" -outfmt 6 -out raw_probes_blastn_homo.txt'
os.system(bla_homo)

# Get non-redundant probe names that have matches to human sequences
df_blast=pd.read_csv('raw_probes_blastn_homo.txt',sep='\t', header=None)

# the number of probes that have matches with human sequences
num_probes=df_blast[0].nunique()
print(str(num_probes) +" probes have matches with human sequences")

# make a dictionary for raw probes with probe id as key and sequence as value
raw_probes_file = open("raw_probes.fasta", 'r')
probe_sequence = {}

for line in raw_probes_file:
    if line.startswith('>'):   
        name=line[1:].strip()
        probe_sequence[name]=''
    else:
        probe_sequence[name]+=line.strip()

raw_probes_file.close()

# Filter the probes and save the rest into new file
filtered_probes_file = open("filtered_probes.fasta", "w")

for key in probe_sequence.keys():
	if key not in df_blast[0].unique():
		filtered_probes_file.write(">"+key+"\n")
		filtered_probes_file.write(probe_sequence[key]+'\n')

filtered_probes_file.close()

print("The number of filtered probes is:")
os.system('grep -o ">" filtered_probes.fasta | wc -l')


# Blast again to check the result
os.system('blastn -db nt -query filtered_probes.fasta -remote -entrez_query "Homo sapiens [Organism]" -outfmt 6')


