if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(dplyr)
library(ggplot2)
library(pheatmap)
library(DESeq2)

#input the data of counts, this csv file just for test
expression_data=read.table(
  "http://www.bioinformatics.nl/courses/RNAseq/ST_vs_HT.csv", row.names=1, header=TRUE, sep =",", stringsAsFactors=FALSE) 

#clean, choose 10 counts as a threshold
mx = apply( expression_data, 1, max )
expression_data = expression_data[ mx > 10, ] 

#create a DESeqDataSet
condition = factor(c("ST","ST","ST","HT","HT","HT"),c("ST","HT"))  #just for testing
col_data = data.frame(condition)
dds = DESeqDataSetFromMatrix(expression_data, col_data, ~condition)

#Normalization
dds = estimateSizeFactors(dds) 
rld = rlog(dds)

#Differential expression
dds = estimateDispersions(dds)
dds = nbinomWaldTest(dds)
res = results(dds)
res$padj = ifelse(is.na(res$padj), 1, res$padj)
resOrdered = res[order(res$padj), ]  #order padj
DEG = as.data.frame(resOrdered)

#filter
DEG_deseq2 = na.omit(DEG) #delete NA

#choose log2FC>1 and FDR(padj)<0.01
filtered_results = subset(DEG_deseq2, log2FoldChange > 1 & padj < 0.01) 

#sorted as p value
sorted_results = filtered_results[order(filtered_results$pvalue), ]
save(sorted_results, file = './DEG_deseq2.Rdata') #save results, position can be changed

#results file
write.table(sorted_results, col.names=NA, row.names=T, file ="expressions_DESeq2.tsv", sep ="\t")

#plot
load("./DEG_deseq2.Rdata")
logFC = 1 #same as paper
P.Value = 0.01  
k1 = (DEG_deseq2$pvalue < P.Value) & (DEG_deseq2$log2FoldChange < -logFC)
k2 = (DEG_deseq2$pvalue < P.Value) & (DEG_deseq2$log2FoldChange > logFC)
DEG_deseq2 = mutate(DEG_deseq2, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
p = ggplot(data = DEG_deseq2, 
            aes(x = log2FoldChange, 
                y = -log10(pvalue))) +
  geom_point(alpha = 0.4, size = 3.5, 
             aes(color = change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values = c("blue4", "grey", "red3"))+
  geom_vline(xintercept = c(-logFC, logFC), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(P.Value), lty = 4, col = "black", lwd = 0.8) +
  theme_bw()
p

ggsave(filename = "./volcano_plot_deseq2.pdf", plot = p, device = "pdf", width = 6, height = 5)
dev.off()

deg_opt = DEG_deseq2 %>% filter(DEG_deseq2$change != "stable")
exp_brca_heatmap = expression_data %>% filter(rownames(expression_data) %in% rownames(deg_opt))
rownames(col_data) = colnames(exp_brca_heatmap) 
p1 = pheatmap(exp_brca_heatmap, show_colnames = F, show_rownames = F,
               scale = "row",
               cluster_cols = F,
               col_data = col_data,
               breaks = seq(-3, 3, length.out = 100)) 
p1
ggsave(filename = "./heatmap_plot_deseq2.pdf", plot = p1, device = "pdf", width = 5, height = 6)
dev.off()



