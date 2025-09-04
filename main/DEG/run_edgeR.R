# load edgeR
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("edgeR")) BiocManager::install("edgeR")
library(edgeR)
library(ggplot2)
library(dplyr)
library(pheatmap)

# assign arguments
args = commandArgs(trailingOnly=T)
ctr_file = args[1]
trt_file = args[2]
output_file1 = args[3]
output_file2 = args[4]
output_file3 = args[5]

# ctr_file = 'input/Axenic_count.tsv'
# trt_file = 'input/Leaf048_count.tsv'
# output_file = 'edgeR_output.tsv'

# print args
print(ctr_file)
print(trt_file)
print(output_file1)
print(output_file2)
print(output_file3)

# import data from files
ctr_data = read.csv(ctr_file, check.names = FALSE, row.names=1, sep=',')
trt_data = read.csv(trt_file, check.names = FALSE, row.names=1, sep=',')

# rename columns according to condition
rename_columns = function(df, condition){
  for (i in 1:length(colnames(df))) {
    colnames(df)[i] = sprintf("%s_%s", condition, i)
  }
  return(df)
}

ctr_data = rename_columns(ctr_data, 'ctr')
trt_data = rename_columns(trt_data, 'trt')

# merge data
exp_data <- transform(merge(ctr_data, trt_data, by = 0, all = TRUE), 
                      row.names=Row.names, Row.names=NULL)

# condition data frame
col_data = matrix(ncol=length(colnames(exp_data)), nrow=1)
for (i in 1:length(colnames(exp_data))){
  col = colnames(exp_data)[i]
  col_data[i] = strsplit(col,split = '_')[[1]][1]
}
col_data = factor(col_data)


# create a DGEList
d = DGEList(counts=exp_data, samples=col_data)

# filter the gene use cpm value
keep = rowSums(cpm(d) > 0) >= 1
d = d[keep, , keep.lib.sizes = FALSE]
d$samples$lib.size <- colSums(d$counts)

# Normalization, using TMM
d = calcNormFactors(d)
dge = d

#Differential expression
design = model.matrix(~0 + col_data)
rownames(design) = colnames(dge)
colnames(design) = levels(col_data)

dge = estimateGLMCommonDisp(dge, design)
dge = estimateGLMTrendedDisp(dge, design)
dge = estimateGLMTagwiseDisp(dge, design)

# The generalized linear model (GLM) is fitted on the basis of the estimated model
fit = glmFit(dge, design)

# -1 for control sample(A-1)，1 for other samples(A2-40)
lrt = glmLRT(fit, contrast = c(-1, 1))

# get top DEGs from LRT results
nrDEG = topTags(lrt, n = nrow(dge))
DEG_edgeR = as.data.frame(nrDEG)

#choose log2FC>1 and FDR<0.01
# filtered_results = subset(DEG_edgeR, logFC > 1 & FDR < 0.01)
filtered_results = DEG_edgeR

#sorted as p value
sorted_results = filtered_results[order(filtered_results$PValue), ]

#results file
write.table(sorted_results, col.names=NA, row.names=T, file=output_file1, sep ="\t")

#save data for better use
save(sorted_results, file = './DEG_edgeR.Rdata')

#plot
load("./DEG_edgeR.Rdata")
logFC = 1
P.Value = 0.01
k1 = (DEG_edgeR$PValue < P.Value) & (DEG_edgeR$logFC < -logFC)
k2 = (DEG_edgeR$PValue < P.Value) & (DEG_edgeR$logFC > logFC)
DEG_edgeR = mutate(DEG_edgeR, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))

p = ggplot(data = DEG_edgeR, 
            aes(x = logFC, 
                y = -log10(PValue))) +
  geom_point(alpha = 0.4, size = 3.5, 
             aes(color = change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values = c("blue4", "grey", "red3"))+
  geom_vline(xintercept = c(-logFC, logFC), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(P.Value), lty = 4, col = "black", lwd = 0.8) +
  theme_bw()
p

ggsave(filename = output_file2, plot = p, device = "png", width = 6, height = 5)
dev.off()


deg_subset <- DEG_edgeR[, 1:5]
p1 = pheatmap(deg_subset, show_colnames = F, show_rownames = F,
              scale = "row",
              cluster_cols = F,
              col_data = col_data,
              breaks = seq(-3, 3, length.out = 100)) 
p1

ggsave(filename = output_file3, plot = p1, device = "png", width = 5, height = 6)
dev.off()
