# Setup script that creates all required directories and softlinks to files on the Leunissen server HPC

# Create directories and sub directories
mkdir samples
mkdir filtered_reads
mkdir mapped_reads
mkdir read_counts
mkdir read_couts/counts
mkdir read_counts/merged_counts
mkdir read_counts/quantified_replicates
mkdir differential_expression
mkdir pathway_clustering
mkdir genome
mkdir genome/hisat2index

# Create softlinks to the Arabidopsis genome
ln -s /prog/BIF30806/genomes/Arabidopsis_thaliana/genes.gtf genome
ln -s /prog/BIF30806/genomes/Arabidopsis_thaliana/genome.fa genome
ln -s /prog/BIF30806/genomes/Arabidopsis_thaliana/genome.fa.fai genome

# Create softlink to samples in the samples directory
ln -s /prog/BIF30806/project/data/transcriptome/*.gz samples

# Rename the samples to remove the date in front of the sample names
cd samples
for f in 20180605.*; do mv "$f" "$(echo "$f" | sed s/20180605.//)"; done
for f in 20180613.*; do mv "$f" "$(echo "$f" | sed s/20180613.//)"; done

# Exclude anomalous sample 21, fourth replicate
rm A-21R4_p2663_4256_189_R1.fastq.gz
# Exclude anomalous sample 38
rm A-38*
