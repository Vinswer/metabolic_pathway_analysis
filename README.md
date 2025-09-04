# metabolic_pathway_analysis
This repository was created in an effort to produce a pipeline that starts with raw RNA-seq reads and ends with an analysis of  differentially expressed genes in pathways. The pipeline was designed specifically for A. thaliana and its bacterial pathogens  but could, in theory, be applied to any RNA-seq analysis. This project is an advanced bioinformatics course project, completed by members of the NucleoVibes team (I was mainly responsible for the differential gene expression part). The data used in the project came from the article: A general non-self response as part of plant immunity, data source: https://www.nature.com/articles/s41477-021-00913-1#Sec27

The pipeline we have developed works in 7 steps:
 
1. Filter low quality/short reads using fastp
2. Map reads using hisat2
3. Compress mapped reads sam --> bam using samtools
4. Quantify the reads using stringtie
5. Produce individual read counts per replicate using a modified stringtie 
   prepDE.py3 python script
6. Concatenate the five replicates per sample to the five control group replicates
   using our own merge.py python script
7. Feed the concatenated read count files into our own run_edgeR.R Rscript to 
   calculate differential expression and pathway enrichment
