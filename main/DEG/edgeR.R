#input merged counts
args <- commandArgs(trailingOnly = TRUE)
samples_str <- args[1]
output_file <- args[2]  
data <- read.csv(input_file, row.names = 1)


#import edgeR
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
library(edgeR)
library(ggplot2)
library(dplyr)
library(pheatmap)

# merge data from files
list_counts <- list()
for (sample in SAMPLES) {
  filename <- paste0(sample, "_merged.csv")
  counts <- read.csv(filename, row.names = 1)
  list_counts[[sample]] <- counts
}



#create DEGList
condition = factor(c("A-1", "A-2", "A-6", "A-9", "A-14", "A-15", "A-18", "A-24", "A-36", "A-39", "A-40")) 
expression_data <- do.call(cbind, list_counts)
col_data = data.frame(condition)
d = DGEList(counts = expression_data, col_data=col_data)

# filter the gene use cpm value
keep = rowSums(cpm(d) > 1) >= 2
d = d[keep, , keep.lib.sizes = FALSE]
d$samples$lib.size <- colSums(d$counts)

# Normalization, using TMM
d = calcNormFactors(d)
dge = d

#Differential expression
design = model.matrix(~0 + condition)
rownames(design) = colnames(dge)
colnames(design) = levels(condition)

dge = estimateGLMCommonDisp(dge, design)
dge = estimateGLMTrendedDisp(dge, design)
dge = estimateGLMTagwiseDisp(dge, design)

# The generalized linear model (GLM) is fitted on the basis of the estimated model
fit = glmFit(dge, design)

# -1 for ST，1 for HT
lrt = glmLRT(fit, contrast = c(-1, 1))

# get top DEGs from LRT results
nrDEG = topTags(lrt, n = nrow(dge))
DEG_edgeR = as.data.frame(nrDEG)

#choose log2FC>1 and FDR<0.01
filtered_results = subset(DEG_edgeR, logFC > 1 & FDR < 0.01)


#sorted as p value
sorted_results = filtered_results[order(filtered_results$PValue), ]
save(sorted_results, file = './DEG_edgeR.Rdata')

#results file
write.csv(sorted_results, file = output_file)

# plot
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

ggsave(filename = "./volcano_plot_edgeR.pdf", plot = p, device = "pdf", width = 6, height = 5)
dev.off()


deg_opt = DEG_edgeR %>% filter(DEG_edgeR$change != "stable")
exp_brca_heatmap = expression_data %>% filter(rownames(expression_data) %in% rownames(deg_opt))
rownames(col_data) = colnames(exp_brca_heatmap) 

p1 = pheatmap(exp_brca_heatmap, show_colnames = F, show_rownames = F,
               scale = "row",
               cluster_cols = F,
               col_data = col_data,
               breaks = seq(-3, 3, length.out = 100)) 

p1
ggsave(filename = "./heatmap_plot_edgeR.pdf", plot = p1, device = "pdf", width = 5, height = 6)
dev.off()
