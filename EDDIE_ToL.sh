# BLAST v2.11.0 always has memory error, so use v2.9.0 to make BLAST database
module load roslin/blast+/2.9.0	

#load conda into the node
module load anaconda

#unzip the file
gzip -d SD001-17-AMYG_R1_001.fastq.gz

#convert from fastq to fasta
sed -n '1~4s/^@/>/p;2~4p' SD001-17-AMYG_R1_001.fastq > SD001-17-AMYG_R1_001.fasta

#make the file into a blast-able database
makeblastdb -in SD001-17-AMYG_R1_001.fasta -dbtype nucl &

module load roslin/blast+/2.11.0

# make a directory to save the result
mkdir SD001-17-AMYG_R1_001

# BLAST against RNA-seq dataset with ToL and control probe list
# The ToL probes were made with probe.py
blastn -db SD001-17-AMYG_R1_001.fasta -query filtered_probesv2.0.fasta -outfmt 6 -num_threads 16 -max_target_seqs 60000 -out ./SD001-17-AMYG_R1_001/filtered_probesv2.0_SD001-17-AMYG_R1_001.txt &
blastn -db SD001-17-AMYG_R1_001.fasta -query probes_control.txt -outfmt 6 -num_threads 16 -max_target_seqs 2800000 -out ./SD001-17-AMYG_R1_001/probes_control_SD001-17-AMYG_R1_001.txt &

# Filter the ToL probes and produce a fasta file containing sequences of matches (reads)
python Abundance_ToL.py SD001-17-AMYG_R1_001

# BLAST against human sequences with sequences of matches
#blastn -db nt -query ./SD001-17-AMYG_R1_001/SD001-17-AMYG_R1_001_matches.fasta -remote -entrez_query "Homo sapiens [Organism]" -outfmt 6 -max_target_seqs 50 -out ./SD001-17-AMYG_R1_001/SD001-17-AMYG_R1_001_matches_blast_homo.txt
#the remote option is too slow, use the downloaded nt database instead
blastn -db nt -query ./SD001-17-AMYG_R1_001/SD001-17-AMYG_R1_001_matches.fasta -taxids 9606 -outfmt 6 -max_target_seqs 10 -out ./SD001-17-AMYG_R1_001/SD001-17-AMYG_R1_001_matches_blast_homo.txt

# For ToL probes, remove the matches which are highly homologous with human, then count the number of matches for each probe
# For control probes, remove the duplicates, then count the number of matches for each probe
python Abundance_count.py SD001-17-AMYG_R1_001
