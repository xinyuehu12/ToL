import os
import sys
import pandas as pd

pd.set_option('display.max_rows', 4)
ndata=sys.argv[1]

# 1 Filter duplicates of BLAST output
# Read the BLAST output with pandas
df = pd.read_csv('./'+ndata+'/filtered_probesv2.0_'+ndata+'.txt',sep='\t', names=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'])

# Check and filter the BLAST result (E value < 0.01) 
df_Evalue=df.loc[df['evalue'] > 0.01]

if df_Evalue.empty:
    print('The E values for all matches are less than 0.01')
    df1=df
else:
	print('The matches with > 0.01 E values will be removed:\n')
	print(df_Evalue)
	df1=df.loc[df['evalue'] < 0.01]


# Sort the BLAST result with bitscore in descending order
# Remove duplicated RNAseq sequences by bitscore
df2=df1.sort_values(by='bitscore' , ascending=False).drop_duplicates('sseqid').sort_values(by='sseqid' , ascending=True)
df2.to_csv('./'+ndata+'/filtered_probesv2.0_'+ndata+'_filtered.txt', sep='\t', index=False)

# 2 Retrieve the sequences of matches
# make a list of match names
matched_id=list(df2['sseqid'])

# make a dictionary for dataset sequences with sequence name as key and sequence as value
dataset = open(''+ndata+'.fasta', 'r')
dataset_sequence = {}

for line in dataset:
    if line.startswith('>'):   
        name=line[1:].strip().split()[0]
        dataset_sequence[name]=''
    else:
        dataset_sequence[name]+=line.strip()

dataset.close()

# Retrieve the sequences of matches with dictionary
match_sequences_file = open('./'+ndata+'/'+ndata+'_matches.fasta', 'w')

for name in matched_id:
	match_sequences_file.write('>'+name+'\n')
	match_sequences_file.write(dataset_sequence[name]+'\n')

match_sequences_file.close()
