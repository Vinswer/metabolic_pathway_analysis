if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("KEGGREST", quietly = TRUE))
  BiocManager::install("KEGGREST")

library(KEGGREST)

# Pull all pathways for AT  
pathways.list <- keggList("pathway", "ath")

# Pull all genes for each pathway
pathway.codes <- sub("path:", "", names(pathways.list)) 
genes.by.pathway <- sapply(pathway.codes,
                           function(pwid){
                             pw <- keggGet(pwid)
                             if (is.null(pw[[1]]$GENE)) return(NA)
                             pw2 <- pw[[1]]$GENE[c(TRUE,FALSE)] # may need to modify this to c(FALSE, TRUE) for other organisms
                             pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
                             return(pw2)
                           }
)

# import data
infile <- "araport11.txt"
tmp <- read.delim(infile)
geneList <- rep(0, length.out=length(tmp$Gene))
names(geneList) <- tmp$Gene


pathways <- sapply(names(genes.by.pathway), 
                   function(pathway) {
                     pathway.genes <- genes.by.pathway[[pathway]]
                     an.genes <- intersect(names(geneList), pathway.genes)
                     an.genes <- paste(an.genes,collapse=";")
                     return(c(description = pathways.list[pathway], genes = an.genes))
                     }
                   )

pathways <- t(pathways)  ## transpose so each pathway is a row
write.table(pathways, 'KEGGlist.tsv', sep="\t") ## write to tsv file

