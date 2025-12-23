low_Fst_genes <- read.delim("low-Fst_genes_ID_geneName.txt", header=T)
gene2GO_herring <- read.delim("gene2GO_herring.map", header=F)

library(tidyverse)
join_herring <- left_join(low_Fst_genes, gene2GO_herring, by = c("Gene.stable.ID" = "V1"))
nrow(join_herring[join_herring$V2 != "",])

gene2GO_herring_NBIS_ortho <- read.delim("gene2GO_herring_NBIS_ortho.map", header=F)

join_herring_NBIS_ortho <- left_join(low_Fst_genes, gene2GO_herring_NBIS_ortho, by = c("Gene.stable.ID" = "V1"))
nrow(join_herring_NBIS_ortho[join_herring_NBIS_ortho$V2 != "",])
