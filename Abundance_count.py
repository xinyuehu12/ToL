import os
import sys
import pandas as pd

pd.set_option('display.max_rows', 4)
ndata=sys.argv[1]

# 1 ToL probes
# 1.1 Check and remove the sequences of matches that are highly similar to human sequences
# If bitscore > 160 (MSBB); 100 (Miami); 126 (Rockefeller); 150 (EBB), the match is human sequences
print('The sequences of ToL matches that are highly similar to human sequences will be checked and removed')
# Load the blast output with match sequences as query and human sequences as subject in pandas
df_blast=pd.read_csv('./'+ndata+'/'+ndata+'_matches_blast_homo.txt',sep='\t', names=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'])
# Save the blast output whose bitscore is > 160
df_bits=df_blast.loc[df_blast['bitscore'] > 160].sort_values(by='bitscore' , ascending=False)

# Load filtered blast result in pandas
df=pd.read_csv('./'+ndata+'/filtered_probesv2.0_'+ndata+'_filtered.txt',sep='\t')
# Get the columns of qseqid and sseqid and change column names into probes and matches
df_match=df.loc[:, ['qseqid', 'sseqid']].sort_values(by='qseqid' , ascending=True)
df_match.columns = ['probes','matches']

# Remove the duplicates by only keeping the hit with highest bitscore
df_num=df_blast.sort_values(by='bitscore' , ascending=False).drop_duplicates('qseqid')
# Get the number of hits in total
num_match=df_num['qseqid'].nunique()
print('The number of matches is '+str(num_match))

# Check and filter the hits
if df_bits.empty:
    print('The bit score for all matches are less than 160 (less or not similar to human sequences)')
    df_match=df_match
else:
	print('The matches with > 160 bit score (highly similar to human sequences) are:\n')
	print(df_bits)
	num_matches=df_bits['qseqid'].nunique()
	print(str(num_matches) +' unique matches are highly similar to human sequences')
	matches_homo=df_bits['qseqid'].unique()
	df_match=df_match[~df_match.matches.isin(matches_homo)]

# Save the filtered result in file
df_match.to_csv('./'+ndata+'/filtered_probesv2.0_matches.txt', sep='\t', index=False)

# 1.2 Abundance analysis
print('Abundance analysis for ToL probes will be performed')
# Count the number of matches for each probe and save in dictionary
probe_match_number=pd.value_counts(df_match['probes']).to_dict()

# Save the result in file
print('The number of matches will be saved.')
abundance_filtered_probes_dataset_file = open('./'+ndata+'/abundance_filtered_probesv2.0_'+ndata+'.txt', 'w')
filtered_probes_file = open('filtered_probesv2.0.fasta', 'r')
abundance_filtered_probes_dataset_file.write('probe\tmatch_count\n')

for line in filtered_probes_file:
	if line.startswith('>'):
		probe_name = line.rstrip()[1:]
		if probe_name in probe_match_number:
			abundance_filtered_probes_dataset_file.write(probe_name+'\t'+str(probe_match_number[probe_name])+'\n')
		else:
			abundance_filtered_probes_dataset_file.write(probe_name+'\t'+'0'+'\n')

abundance_filtered_probes_dataset_file.close()
filtered_probes_file.close()

# 2 Probes for housekeeping genes, mitochondria, contaminations, viruses and retroelements...

# 2.1 filter the duplicates of BLAST output
print('The duplicates of control probelist will be checked and removed')
# Read the BLAST output with pandas
df = pd.read_csv('./'+ndata+'/probes_control_'+ndata+'.txt',sep='\t', names=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'])

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
print('Remove duplicated RNAseq sequences by bitscore')
df_blast_control=df1.sort_values(by='bitscore' , ascending=False).drop_duplicates('sseqid').sort_values(by='sseqid' , ascending=True)
df_blast_control.to_csv('./'+ndata+'/probes_control_'+ndata+'_filtered.txt', sep='\t', index=False)

# Save the matches for each probe in file
df_match_control=df_blast_control[['qseqid','sseqid']]
df_match_control.columns= ['probes','matches']

df_match_control.to_csv('./'+ndata+'/probes_control_matches.txt', sep='\t', index=False)

# 2.2 Abundance analysis
print('Abundance analysis for control probes will be performed')
# Count the number of matches for each probe and save in dictionary
control_match_number=pd.value_counts(df_match_control['probes']).to_dict()

# Save the result in file
abundance_probe_control_dataset_file = open('./'+ndata+'/abundance_probe_control_'+ndata+'_filtered.txt', 'w')
probes_control_file = open('probes_control.txt', 'r')
abundance_probe_control_dataset_file.write('probe\tmatch_count\n')

for line in probes_control_file:
	if line.startswith('>'):
		probe_name = line.rstrip()[1:]
		if probe_name in control_match_number:
			abundance_probe_control_dataset_file.write(probe_name+'\t'+str(control_match_number[probe_name])+'\n')
		else:
			abundance_probe_control_dataset_file.write(probe_name+'\t'+'0'+'\n')

abundance_probe_control_dataset_file.close()
probes_control_file.close()
