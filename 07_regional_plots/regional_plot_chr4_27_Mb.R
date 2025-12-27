## Import NCBI gff file and extract the genes in the region of interest
library(tidyverse)
library(rtracklayer)

gff_NCBI <- import("GCF_900700415.2_Ch_v2.0.2_genomic.gff")
gff_NCBI_chr4 <- data.frame(gff_NCBI) %>%
  mutate(Parent = as.character(Parent)) %>%
  select(seqnames, start, end, width, strand, source, type, ID, Name, gene, gene_biotype, partial, Parent, product, transcript_id, standard_name, pseudo, description) %>%
  filter(seqnames == "NC_045155.1", end > 26660000, start < 27590000, type %in% c("gene", "mRNA", "pseudogene")) %>%
  mutate(seqnames = 4)
write.csv(gff_NCBI_chr4, "genes_chr4.csv", row.names = F)

#-------------------------------------------------------------------------------
## Generate plots

#install.packages("gggenes")
#install.packages("RColorBrewer")
library(tidyverse)
library(cowplot)
library(gggenes)
require(RColorBrewer)
head(example_genes)
head(example_subgenes)

genes_in <- read.csv("genes_chr4.csv", header=T)

genes <- genes_in %>% 
  filter(seqnames == 4, end > 26660000, start < 27590000, gene_biotype %in% c("protein_coding", "pseudogene"), type %in% c("gene", "pseudogene")) %>%
  mutate(start = ifelse(start < 26660000, 26660000, start)) %>%
  mutate(end = ifelse(end > 27590000, 27590000, end)) %>%
  mutate(strand = ifelse(strand=="+", 1, -1))


set.seed(123)
color_brewer_vector <- c(brewer.pal(n=11, name="Paired")[1:11], brewer.pal(n=8, name="Dark2")[c(8:6)])
color_brewer_sample <- sample(color_brewer_vector, 14,replace = F)
color_brewer <- c("#B2DF8A", "#FB9A99", "#1F78B4", "#A6761D", "#CAB2D6", "#FF7F00", "#666666", "#A6CEE3", "#33A02C","#FFFF99", "#6A3D9A", "#E31A1C", "#FDBF6F")

q <- ggplot(genes, aes(xmin = start/1000000, xmax = end/1000000, fill = gene_name_final, color = gene_name_final, y = y_axis_new, forward = strand)) +
  geom_gene_arrow() +
  scale_x_continuous(expand = c(0,0), limits = c(26.66, 27.59), breaks = seq(26.7, 27.5,.1)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 3)) +
  scale_color_manual(values = color_brewer, limits = c("CEP104", "LRRC47", "MAD2L2", "CCDC27", "Uncharacterized", "TNFRSF gene family", "Pseudogene", "L1CAM", "SMC1AL", "PFKFB3", "NDNL2", "ITIH6", "GDI1"), labels = c(expression(bolditalic("CEP104")), expression(bolditalic("LRRC47")), expression(bolditalic("MAD2L2")), expression(bolditalic("CCDC27")), expression(bold("Uncharacterized")), expression(paste(bolditalic("TNFRSF"), bold(" gene family"))), expression(bold("Pseudogene")), expression(bolditalic("L1CAM")), expression(bolditalic("SMC1AL")), expression(bolditalic("PFKFB3")), expression(bolditalic("NDNL2")), expression(bolditalic("ITIH6")), expression(bolditalic("GDI1")))) +
  scale_fill_manual(values = color_brewer, limits = c("CEP104", "LRRC47", "MAD2L2", "CCDC27", "Uncharacterized", "TNFRSF gene family", "Pseudogene", "L1CAM", "SMC1AL", "PFKFB3", "NDNL2", "ITIH6", "GDI1"), labels = c(expression(bolditalic("CEP104")), expression(bolditalic("LRRC47")), expression(bolditalic("MAD2L2")), expression(bolditalic("CCDC27")), expression(bold("Uncharacterized")), expression(paste(bolditalic("TNFRSF"), bold(" gene family"))), expression(bold("Pseudogene")), expression(bolditalic("L1CAM")), expression(bolditalic("SMC1AL")), expression(bolditalic("PFKFB3")), expression(bolditalic("NDNL2")), expression(bolditalic("ITIH6")), expression(bolditalic("GDI1")))) +
  labs(x = "Chromosome 4 (Mb)", y = "Genes") +
  geom_hline(yintercept=0, linewidth=0.5) +
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.key.size = unit(.5, "lines"),
        legend.text = element_text(size = 12),
        legend.margin = margin(0, 0, 0, 0),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.direction = "horizontal",
        panel.grid = element_blank(),
        axis.title.x = element_text(size=14, face = "bold"),
        axis.title.y = element_text(size=14, face = "bold"),
        axis.text.x = element_text(size=12),
        #axis.text = element_text(size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 2), color = guide_legend(nrow = 2))


#------------------------------------------

diversity_in <- read.csv("diversity_AtlPac.csv", header=T)
diversity_in <- diversity_in[!is.na(diversity_in$fst_Atlantic.Pacific),]
diversity_in$fst_Atlantic.Pacific[diversity_in$fst_Atlantic.Pacific < 0] <- 0

diversity <- diversity_in %>% 
  mutate(status = ifelse(fst_Atlantic.Pacific < 0.15990814 & Atlantic > 0.007347820 & Pacific > 0.006718794, "Selected", "Unselected"))

diversity <- diversity %>% filter(chromosome == "chr4" & window_pos_1 > 26660000 & window_pos_2 <= 27590003)

p1 <-  ggplot() + 
  geom_rect(aes(xmin=((diversity$window_pos_1[diversity$status == "Selected"] [1]) + 2500)/1000000, xmax=((diversity$window_pos_2[diversity$status == "Selected"] [length(diversity$window_pos_2[diversity$status == "Selected"])]) + 2500)/1000000), ymin=0, ymax=Inf, fill = "red", alpha = .15) +
  geom_line(data = diversity, aes(x = (window_pos_1 + 2500)/1000000, y = fst_Atlantic.Pacific)) +
  labs(y = expression(bolditalic(F[ST]))) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.25))+ 
  scale_x_continuous(expand = c(0,0), limits = c(26.66, 27.59), breaks = seq(26.7, 27.5,.1)) + 
  geom_hline(yintercept=mean(diversity_in$fst_Atlantic.Pacific), linewidth=0.4, linetype = "dashed", color ="black", alpha = 1) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.ticks.y = element_line(linewidth=0.2))


diversity_long <- diversity %>%
  pivot_longer(cols = c(dxy_Atlantic.Pacific, Pacific, Atlantic), names_to = "stat", values_to = "value") %>%
  mutate(stat = factor(stat, levels = c("dxy_Atlantic.Pacific", "Pacific", "Atlantic"), labels = c("dxy_Atlantic.Pacific", "Pacific", "Atlantic")))

p2 <-  ggplot() + 
  geom_rect(aes(xmin=((diversity$window_pos_1[diversity$status == "Selected"] [1]) + 2500)/1000000, xmax=((diversity$window_pos_2[diversity$status == "Selected"] [length(diversity$window_pos_2[diversity$status == "Selected"])]) + 2500)/1000000), ymin=0, ymax=Inf, fill = "red", alpha = .15) +
  geom_line(data = diversity_long, aes(x = (window_pos_1 + 2500)/1000000, y = value, color=stat)) +
  labs(y = expression(bolditalic(pi))) +
  scale_x_continuous(expand = c(0,0), limits = c(26.66, 27.59), breaks = seq(26.7, 27.5,.1)) + 
  scale_color_manual(values = c("dxy_Atlantic.Pacific" = "black", "Pacific" = "green3", "Atlantic" = "dodgerblue"), guide = guide_legend(reverse = TRUE), labels = c("dxy_Atlantic.Pacific"=expression(italic(d[xy])), "Pacific"=expression(italic(pi[Pacific])), "Atlantic"=expression(italic(pi[Atlantic])))) + 
  geom_hline(yintercept=mean(diversity_in$dxy_Atlantic.Pacific), linewidth=0.4, linetype = "dashed", color ="black", alpha = 1) +
  geom_hline(yintercept=mean(diversity_in$Pacific), linewidth=0.4, linetype = "dashed", color ="green3", alpha = 1) +
  geom_hline(yintercept=mean(diversity_in$Atlantic), linewidth=0.4, linetype = "dashed", color ="dodgerblue", alpha = 1) +
  theme_bw() +
  theme(legend.position = c(.99, .99),
        legend.justification = c("right", "top"),
        legend.margin = margin(0, 0, 0, 0),
        legend.title = element_blank(),
        legend.key.width = unit(.4, 'cm'),
        legend.text = element_text(size=13, face="bold"),
        legend.text.align = 0,
        legend.direction = "horizontal",
        legend.background=element_rect(fill = alpha("white", 0)),
        legend.key=element_rect(fill = alpha("white", 0)),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15, face="bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.ticks.y = element_line(linewidth=0.2))


#------------------------------------------

p <- plot_grid(p1, p2, q, ncol=1, rel_heights=c(1, 1, 1.6), align = "v") + theme(plot.margin = margin(2, 0, 0, 0))
png(filename="plot.png", width = 2000, height = 950, units = 'px',res=195) 
print(p)
dev.off()
