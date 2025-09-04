if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("topGO")
BiocManager::install("org.At.tair.db")

library(topGO)

# get data form file
infile <- "araport11.txt"
tmp <- read.delim(infile)

# create topGOData object
geneList <- tmp$P.Value
names(geneList) <- tmp$Gene
GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              geneSelectionFun = function(x)x,
              annot = annFUN.org, mapping = "org.At.tair.db")

# link go terms to annotation
# sel.terms <- sample(usedGO(GOdata), 10)
sel.terms <- usedGO(GOdata)
ann.genes <- genesInTerm(GOdata, sel.terms) ## get the annotations

# turn list of genes into singular string
for (i in 1:length(ann.genes)) {
  ann.genes[[i]] <- paste(unlist(ann.genes[[i]]),collapse=",")
}
ann.genes <- t(t(ann.genes))  # double transpose to prevent write.csv error
write.csv(ann.genes, 'GOlist.csv')  # write to csv file
