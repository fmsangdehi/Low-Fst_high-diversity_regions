setwd("input_to_create_gene2GO_map")
library(topGO)
library(Rgraphviz)
library(tidyverse)

GO_herring <- read.table("mart_export_herring-GO_protein-coding.txt", header=T, sep="\t")
length(unique(GO_herring$Gene.stable.ID))
GO_NBIS <- read.table("gene2GO_NBIS_Teleostei.map", header=F, sep="\t")
colnames(GO_NBIS) <- colnames(GO_herring)
GO_ortho <- read.table("ortho_GO.txt", header=T, sep="\t")
GO_ortho <- GO_ortho[GO_ortho$Gene.stable.ID %in% GO_herring$Gene.stable.ID,]

GO_combined <- rbind(GO_herring, GO_NBIS, GO_ortho)
length(unique(GO_combined$Gene.stable.ID))

# First, combine all GO per gene
gene2GO <- data.frame(unique_genes = unique(GO_combined$Gene.stable.ID))
nrow(gene2GO)

for (i in 1:nrow(gene2GO)){
    gene2GO$GO[i] <- paste(GO_combined$GO.term.accession[GO_combined$Gene.stable.ID == gene2GO$unique_genes[i] & GO_combined$GO.term.accession != ""], collapse=", ")
}

# Then, remove duplicate GO terms per gene
gene2GO$GO_unique <- sapply(gene2GO$GO, function(x) paste(unique(unlist(str_split(x, ", "))), collapse = ", "))
gene2GO$GO <- NULL

write.table(gene2GO,"gene2GO_herring_NBIS_ortho.map",col.names=F,row.names=F,quote=F, sep="\t")
