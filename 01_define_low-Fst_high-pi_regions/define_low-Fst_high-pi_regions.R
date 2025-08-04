library(tidyverse)
library(GenomicRanges)
library(rtracklayer)

diversity <- read.csv("diversity_AtlPac.csv", header=T)
diversity <- diversity[!is.na(diversity$fst_Atlantic.Pacific),]
diversity$fst_Atlantic.Pacific[diversity$fst_Atlantic.Pacific < 0] <- 0
diversity$ZFst <- scale(diversity$fst_Atlantic.Pacific)


quantile(diversity$fst_Atlantic.Pacific, probs=c(.01, .02, .03, .05), na.rm = T)
#quantile(diversity$dxy_Atlantic.Pacific, probs=c(.99, .98, .97, .95), na.rm = T)
quantile(diversity$Atlantic, probs=c(.99, .98, .97, .95), na.rm = T)
quantile(diversity$Pacific, probs=c(.99, .98, .97, .95), na.rm = T)


##### extract low-Fst_high-pi regions
low_Fst <- diversity %>%
           filter(fst_Atlantic.Pacific < 0.15990814 & Atlantic > 0.007347820 & Pacific > 0.006718794)


##### reduce low-Fst_high-pi regions
##### 1_without min.gapwidth & filtering
low_Fst_Granges <- GRanges(seqnames=low_Fst$chromosome, ranges=IRanges(start=low_Fst$window_pos_1, end=low_Fst$window_pos_2))
low_Fst_reduced<-reduce(low_Fst_Granges, with.revmap =T)
write.csv(low_Fst_reduced, "low-Fst_high-pi_intervals_reduced_nogapwidth_nofilter.csv", row.names=F)

##### 2_with min.gapwidth & filtering
low_Fst_Granges <- GRanges(seqnames=low_Fst$chromosome, ranges=IRanges(start=low_Fst$window_pos_1, end=low_Fst$window_pos_2))
low_Fst_reduced_2<-reduce(low_Fst_Granges, with.revmap =T, min.gapwidth = 12000)
low_Fst_reduced_filtered <- low_Fst_reduced_2[low_Fst_reduced_2@ranges@width > 12000]
write.csv(low_Fst_reduced_filtered, "low-Fst_high-pi_intervals_reduced_filtered.csv", row.names=F)

#----------------------------------------------------------------------

##### mark the low-Fst bins in the original dataframe
##### 1_without min.gapwidth & filtering
diversity.Fst.status <- diversity %>%
   mutate(Fst_status = ifelse(fst_Atlantic.Pacific < 0.15990814 & Atlantic > 0.009892208 & Pacific > 0.010216032, "LOW", "NORMAL"))
#diversity.Fst.status$Fst_status <- factor(diversity.Fst.status$Fst_status, levels = c("NORMAL", "LOW"))
write.csv(diversity.Fst.status,"diversity.Fst.status.csv",row.names=F)


##### 2_reduced and filtered low-Fst regions
all_bins<-diversity[,c(1,2,3)]
all_bins_Granges<-GRanges(seqnames=all_bins$chromosome, ranges=IRanges(start=all_bins$window_pos_1, end=all_bins$window_pos_2))
overlaps <- findOverlaps(all_bins_Granges, low_Fst_reduced_filtered)
low_Fst_bins <- all_bins[queryHits(overlaps), ]
low_Fst_bins$index <- paste(low_Fst_bins$chromosome, low_Fst_bins$window_pos_1, sep = "_")
diversity.Fst.status_filtered <- diversity %>%
                   mutate(index = paste(chromosome, window_pos_1, sep = "_")) %>% 
                   mutate(Fst_status = ifelse(index %in% low_Fst_bins$index, "LOW", "NORMAL"))
#diversity.Fst.status_filtered$Fst_status <- factor(diversity.Fst.status_filtered$Fst_status, levels = c("NORMAL", "LOW"))
diversity.Fst.status_filtered$index <- NULL ##or: use select() function from tidyverse package
write.csv(diversity.Fst.status_filtered,"diversity.Fst.status_filtered.csv",row.names=F)

#----------------------------------------------------------------------

##### manhattan plot
diversity.Fst.status <- diversity.Fst.status %>%
          separate(col=chromosome, into=c(NA, "chromosome"), sep=3) %>%
          mutate(chromosome= as.numeric(chromosome)) %>%
          arrange(chromosome, desc(Fst_status), window_pos_1) %>%
          mutate(color=ifelse(Fst_status == "LOW", "red2", ifelse(chromosome %in% seq(2, 26, 2), "gray30", "gray65")))

p <- ggplot(diversity.Fst.status, aes(window_pos_1, fst_Atlantic.Pacific)) + 
         facet_grid(. ~ chromosome, scales = "free_x", switch = "x", space = "free_x")  +
         geom_point(aes(colour = factor(color)), size=.6) +
         scale_colour_manual(values = c("#183059", "#276FBF","red2")) +
         labs(x="Chromosomes", y="Fst") +
         scale_y_continuous(expand = c(.01, .01)) +
         theme(legend.position="none",
              panel.border=element_blank(),
              panel.grid = element_blank(),
              axis.title.y = element_text(size=16),
              axis.title.x = element_text(size=14),
              axis.text = element_text(size=12),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank())

pdf("Manhattan_plot.pdf", width = 20, height = 5) 
print(p)
dev.off()
