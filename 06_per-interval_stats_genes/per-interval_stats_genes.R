library(GenomicRanges)
library(rtracklayer)
library(tidyverse)

########## low-Fst intervals + stats
intervals <- read.csv("low-Fst_high-pi_intervals.csv", header=T)
intervals <- intervals %>%
  mutate(seqnames = str_replace(intervals$seqnames, pattern = "chr", replacement = ""))

stats <- read.csv("/path/to/diversity_AtlPac.csv", header=T)
stats <- stats %>%
  mutate(chromosome = str_replace(stats$chromosome, pattern = "chr", replacement = ""))

options(scipen = 999)

for (i in 1:nrow(intervals)){
  subset <- stats[stats$chromosome == intervals$seqnames[i] & stats$window_pos_1 >= intervals$start[i] & stats$window_pos_2 <= intervals$end[i],]
  intervals$pi_Atlantic[i] <- round(mean(subset$Atlantic, na.rm=T), digits=10)
  intervals$pi_Pacific[i] <- round(mean(subset$Pacific, na.rm=T), digits=10)
  intervals$Fst_AtlPac[i] <- round(mean(subset$fst_Atlantic.Pacific, na.rm=T), digits=10)
  intervals$Dxy_AtlPac[i] <- round(mean(subset$dxy_Atlantic.Pacific, na.rm=T), digits=10)
}


########## (low-Fst intervals + stats) + genes
intervals_GRanges <- GRanges(seqnames = intervals$seqnames, ranges = IRanges(start = intervals$start, end = intervals$end))
mcols(intervals_GRanges) <- intervals[, 5:8]

genes <- read.delim("merged_Ensembl_NBIS_geneName.txt", header=T)
genes_GRanges <- GRanges(seqnames = genes$seqnames, ranges = IRanges(start = genes$start, end = genes$end), strand = genes$strand)
mcols(genes_GRanges) <- genes[, c(1, 7:12)]

for (i in 1:nrow(intervals)){
  subset <- subsetByOverlaps(x=genes_GRanges, ranges=intervals_GRanges[i])
  intervals$coding_gene_count[i] <- length(subset)
  intervals$Gene.stable.ID[i] <- paste(data.frame(subset)$Gene.stable.ID, collapse="\n")
  intervals$Gene.name_NBIS.[i] <- paste(data.frame(subset)$Gene.name_NBIS., collapse="\n")
  intervals$gene.description[i] <- paste(data.frame(subset)$Description..product., collapse="\n")
  intervals$Uniprot.entry[i] <- paste(data.frame(subset)$Uniprot.entry, collapse="\n")
  intervals$Gene.start[i] <- paste(data.frame(subset)$start, collapse="\n")
  intervals$Gene.end[i] <- paste(data.frame(subset)$end, collapse="\n")
  intervals$Gene.width[i] <- paste(data.frame(subset)$width, collapse="\n")
  intervals$strand[i] <- paste(data.frame(subset)$strand, collapse="\n")
}

write.csv(intervals, "per-interval_stats_genes.csv", row.names=F)
