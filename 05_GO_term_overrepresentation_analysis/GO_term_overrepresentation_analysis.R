setwd("GO_term_overrepresentation_analysis")
library(topGO)
library(Rgraphviz)
library(tidyverse)

########## loading gene2GO
geneID2GO<-readMappings("gene2GO_herring_NBIS_ortho.map")
str(head(geneID2GO))

########## creating a list of interesting genes
genes<-read.delim("low-Fst_genes_ID_geneName.txt",header=T)
geneNames <- names(geneID2GO)
geneList <- factor(as.integer(geneNames %in% genes$Gene.stable.ID))
names(geneList) <- geneNames

########## creating topGOdata object and running the enrichment test
########## BP
GOdata_BP <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
GOdata_BP
head(termStat(GOdata_BP))

result_classic_fisher_BP<-runTest(GOdata_BP, algorithm = "classic", statistic = "fisher")
table_classic_fisher_BP<-GenTable(GOdata_BP, classicFisher = result_classic_fisher_BP,topNodes = 50, numChar = 1000)

result_weight01_fisher_BP<-runTest(GOdata_BP, algorithm = "weight01", statistic = "fisher")
table_weight01_fisher_BP<-GenTable(GOdata_BP, weight01Fisher = result_weight01_fisher_BP,topNodes = 50, numChar=1000)
table_weight01_fisher_BP_adjusted <- table_weight01_fisher_BP %>% 
    mutate(weight01Fisher = as.numeric(table_weight01_fisher_BP$weight01Fisher) * 6690)
sigNodes_weight01 <- nrow(table_weight01_fisher_BP_adjusted[as.numeric(table_weight01_fisher_BP_adjusted$weight01Fisher) < 0.05,])
printGraph(GOdata_BP, result_weight01_fisher_BP, firstSigNodes = sigNodes_weight01, fn.prefix = "BP", useInfo = "all", pdfSW = TRUE)

result_parentchild_fisher_BP<-runTest(GOdata_BP, algorithm = "parentchild", statistic = "fisher")
table_parentchild_fisher_BP<-GenTable(GOdata_BP, parentchildFisher = result_parentchild_fisher_BP,topNodes = 50)
sigNodes_parentchild <- length(table_parentchild_fisher_BP$parentchildFisher[as.numeric(table_parentchild_fisher_BP$parentchildFisher) < 0.05])

### comparison of the results of different methods for BP
table_BP<-GenTable(GOdata_BP, parentchildFisher = result_parentchild_fisher_BP,classicFisher = result_classic_fisher_BP,weight01Fisher = result_weight01_fisher_BP,orderBy = "weight01Fisher", ranksOf = "parentchildFisher", topNodes = 50, numChar=1000)
table_BP$Expected <- round(table_BP$Expected, digits=1)
table_BP$weight01Fisher <- formatC(as.numeric(table_BP$weight01Fisher), format = "e", digits = 1)
table_BP$parentchildFisher <- formatC(as.numeric(table_BP$parentchildFisher), format = "e", digits = 1)
table_BP$classicFisher <- formatC(as.numeric(table_BP$classicFisher), format = "e", digits = 1)

### retrieve significant genes for each GO ID
genes$Gene.name[genes$Gene.name == ""] <- genes$Gene.stable.ID[genes$Gene.name == ""]

for(i in 1:nrow(table_BP)){
AnnotatedGenes <- as.character(unlist(genesInTerm(object = GOdata_BP, whichGO = table_BP$GO.ID[i])))
SignificantGenes <- intersect(AnnotatedGenes, genes$Gene.stable.ID)
table_BP[i, "SignificantGenes"] <- paste(SignificantGenes, collapse=";")
SigGeneName <- inner_join(as.data.frame(SignificantGenes), genes, by=c("SignificantGenes" = "Gene.stable.ID"))
table_BP[i, "SigGeneName"] <- paste(SigGeneName$Gene.name, collapse=";")
}
write.table(table_BP, "table_BP.txt", quote=F, row.names=F, sep="\t")


########## MF
GOdata_MF <- new("topGOdata", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
GOdata_MF
head(termStat(GOdata_MF))

result_classic_fisher_MF<-runTest(GOdata_MF, algorithm = "classic", statistic = "fisher")
table_classic_fisher_MF<-GenTable(GOdata_MF, classicFisher = result_classic_fisher_MF,topNodes = 50)

result_weight01_fisher_MF<-runTest(GOdata_MF, algorithm = "weight01", statistic = "fisher")
table_weight01_fisher_MF<-GenTable(GOdata_MF, weight01Fisher = result_weight01_fisher_MF,topNodes = 50)
sigNodes_weight01 <- length(table_weight01_fisher_MF$weight01Fisher[as.numeric(table_weight01_fisher_MF$weight01Fisher) < 0.05])
printGraph(GOdata_MF, result_weight01_fisher_MF, firstSigNodes = sigNodes_weight01, fn.prefix = "MF", useInfo = "all", pdfSW = TRUE)

result_parentchild_fisher_MF<-runTest(GOdata_MF, algorithm = "parentchild", statistic = "fisher")
table_parentchild_fisher_MF<-GenTable(GOdata_MF, parentchildFisher = result_parentchild_fisher_MF,topNodes = 50)
sigNodes_parentchild <- length(table_parentchild_fisher_MF$parentchildFisher[as.numeric(table_parentchild_fisher_MF$parentchildFisher) < 0.05])

### comparison of the results of different methods for MF
table_MF<-GenTable(GOdata_MF, parentchildFisher = result_parentchild_fisher_MF,classicFisher = result_classic_fisher_MF,weight01Fisher = result_weight01_fisher_MF,orderBy = "weight01Fisher", ranksOf = "parentchildFisher", topNodes = 50, numChar=1000)
table_MF$weight01Fisher <- formatC(as.numeric(table_MF$weight01Fisher), format = "e", digits = 1)
table_MF$Expected <- round(table_MF$Expected, digits=1)

### retrieve significant genes for each GO ID (automated)
for(i in 1:nrow(table_MF)){
AnnotatedGenes <- as.character(unlist(genesInTerm(object = GOdata_MF, whichGO = table_MF$GO.ID[i])))
SignificantGenes <- intersect(AnnotatedGenes, genes$Gene.stable.ID)
SigGeneName <- inner_join(as.data.frame(SignificantGenes), genes, by=c("SignificantGenes" = "Gene.stable.ID"))
table_MF[i, "SigGeneName"] <- paste(SigGeneName$Gene.name, collapse=";")
}
write.table(table_MF, "table_MF.txt", quote=F, row.names=F, sep="\t")


########## CC
GOdata_CC <- new("topGOdata", ontology = "CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
GOdata_CC
head(termStat(GOdata_CC))

result_classic_fisher_CC<-runTest(GOdata_CC, algorithm = "classic", statistic = "fisher")
table_classic_fisher_CC<-GenTable(GOdata_CC, classicFisher = result_classic_fisher_CC,topNodes = 50)

result_weight01_fisher_CC<-runTest(GOdata_CC, algorithm = "weight01", statistic = "fisher")
table_weight01_fisher_CC<-GenTable(GOdata_CC, weight01Fisher = result_weight01_fisher_CC,topNodes = 50)
sigNodes_weight01 <- length(table_weight01_fisher_CC$weight01Fisher[as.numeric(table_weight01_fisher_CC$weight01Fisher) < 0.05])
printGraph(GOdata_CC, result_weight01_fisher_CC, firstSigNodes = sigNodes_weight01, fn.prefix = "CC", useInfo = "all", pdfSW = TRUE)

result_parentchild_fisher_CC<-runTest(GOdata_CC, algorithm = "parentchild", statistic = "fisher")
table_parentchild_fisher_CC<-GenTable(GOdata_CC, parentchildFisher = result_parentchild_fisher_CC,topNodes = 50)
sigNodes_parentchild <- length(table_parentchild_fisher_CC$parentchildFisher[as.numeric(table_parentchild_fisher_CC$parentchildFisher) < 0.05])

### comparison of the results of different methods for CC
table_CC<-GenTable(GOdata_CC, parentchildFisher = result_parentchild_fisher_CC,classicFisher = result_classic_fisher_CC,weight01Fisher = result_weight01_fisher_CC,orderBy = "weight01Fisher", ranksOf = "parentchildFisher", topNodes = 50, numChar=1000)
table_CC$weight01Fisher <- formatC(as.numeric(table_CC$weight01Fisher), format = "e", digits = 1)
table_CC$Expected <- round(table_CC$Expected, digits=1)

### retrieve significant genes for each GO ID (automated)
for(i in 1:nrow(table_CC)){
AnnotatedGenes <- as.character(unlist(genesInTerm(object = GOdata_CC, whichGO = table_CC$GO.ID[i])))
SignificantGenes <- intersect(AnnotatedGenes, genes$Gene.stable.ID)
SigGeneName <- inner_join(as.data.frame(SignificantGenes), genes, by=c("SignificantGenes" = "Gene.stable.ID"))
table_CC[i, "SigGeneName"] <- paste(SigGeneName$Gene.name, collapse=";")
}
write.table(table_CC, "table_CC.txt", quote=F, row.names=F, sep="\t")


#---------------------------------------------------------------------

########## dotplot
########## BP_weight01_fisher (cutoff P<0.01)

head(table_weight01_fisher_BP)
table_weight01_fisher_BP$weight01Fisher[table_weight01_fisher_BP$weight01Fisher == "< 1e-30"] <- "1e-30"
table_weight01_fisher_BP$weight01Fisher <- as.numeric(table_weight01_fisher_BP$weight01Fisher)
table_weight01_fisher_BP$weight01Fisher <- table_weight01_fisher_BP$weight01Fisher * 6690
table_weight01_fisher_BP<-table_weight01_fisher_BP[table_weight01_fisher_BP$weight01Fisher<0.05,]
table_weight01_fisher_BP$color <- round(1-(table_weight01_fisher_BP$Significant/table_weight01_fisher_BP$Annotated), digits=1)
require(ggplot2)

pdf("BP_weight01_dotplot.pdf", width = 8, height = 6)
ggplot(table_weight01_fisher_BP, aes(reorder(x = Term,-weight01Fisher), y = -log10(weight01Fisher), size = Significant)) +
    geom_point(shape = 21, colour = "black", fill = "#8BA0E4", alpha = .9) + # steelblue (alpha 0.7), chocolate1, orange1
    scale_size_area(max_size = 10, name="Gene count", breaks = seq(10, 70, 20), labels = seq(10, 70, 20)) +
        xlab("") +
    ylab(expression(bold(-log[10] (p.adjust)))) +
    ggtitle("") +
    scale_y_continuous(limits = c(0, 14), expand = c(0, 0)) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(color="grey90",size=0.3),
        panel.grid.major.x = element_line(color="grey90",size=0.3),
        axis.text.x=element_text(angle=0, size=10, face="bold", hjust=0.5, color = "black"),
        axis.text.y=element_text(angle=0, size=11, face="bold", vjust=0.5, color = "black"),
        axis.title=element_text(size=11, face="bold"),
        title=element_text(size=1)) +
    theme(legend.position = "right",
          legend.key.size = unit(.6, "lines"),
          legend.text = element_text(size = 8, face="bold"),
          legend.margin = margin(t = 0, unit='cm'),
          legend.title = element_text(size=9, face="bold")) +
    coord_flip()
dev.off()
