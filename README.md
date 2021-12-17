# The ToL approach
The scripts were generated for electronic tree of life (eTOL) approach to identify the presence of microbes in human brains. 

### Make 64-mer probes
The probe.py was used to download the rRNA sequences of ToL organisms and make 64-mer probes. The raw 64-mer probe list was filtered to remove the probes that are homologous to human sequences by BLAST.

### Abundance analysis 
The EDDIE_ToL.sh uses the RNAseq dataset, SD001-17-AMYG_R1_001.fastq, as an example to show how the microbial abundance analysis was performed. The fastq format was converted to fasta format to make a blast-able database. Then, BLAST was performed against RNA-seq dataset with ToL and control probe list, respectively. The RNAseq sequence may have matches to different ToL probes, so the Abundance_ToL.py was used to remove the duplicates and one RNAseq sequence can only be allocated to one probe. The matched RNAseq sequences were retrieved for BLAST against human sequences. The Abundance_count.py was generated to remove the human homologous matches and count the number of matches for each ToL probe. For the control probes, the duplicates were removed and the number of matches for each probe was counted.
