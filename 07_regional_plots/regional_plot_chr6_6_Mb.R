## Import NCBI gff file and extract the genes in the region of interest
library(tidyverse)
library(rtracklayer)

gff_NCBI <- import("GCF_900700415.2_Ch_v2.0.2_genomic.gff")
gff_NCBI_chr6 <- data.frame(gff_NCBI) %>%
  mutate(Parent = as.character(Parent)) %>%
  select(seqnames, start, end, width, strand, source, type, ID, Name, gene, gene_biotype, partial, Parent, product, transcript_id, standard_name, pseudo, description) %>%
  filter(seqnames == "NC_045157.1", end > 5800000, start < 7200000, type %in% c("gene", "mRNA")) %>%
  mutate(seqnames = 6)
write.csv(gff_NCBI_chr6, "genes_chr6.csv", row.names = F)

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

genes_in <- read.csv("genes_chr6.csv", header=T)

genes <- genes_in %>% 
  filter(seqnames == 6, end > 5930000, start < 7070000, gene_biotype == "protein_coding", type == "gene") %>%
  mutate(start = ifelse(start < 5930000, 5930000, start)) %>%
  mutate(end = ifelse(end > 7070000, 7070000, end)) %>%
  mutate(strand = ifelse(strand=="+", 1, -1))


q <- ggplot(genes, aes(xmin = start/1000000, xmax = end/1000000, fill = gene_name_final, color = gene_name_final, y = 1, forward = strand)) +
  geom_gene_arrow() +
  scale_x_continuous(expand = c(0,0), limits = c(5.93, 7.07), breaks = seq(6,7,.2)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 2)) +
  scale_color_manual(values = brewer.pal(n=8, name="Set1")[c(7,2,5,6,3,4,8,1)], limits = c("PIGG", "STAG3", "GPC2", "ITGA6", "NITR", "C1QL4", "ELA1", "RNASEL3"), labels = c(expression(bolditalic("PIGG")), expression(bolditalic("STAG3")), expression(bolditalic("GPC2")), expression(bolditalic("ITGA6")), expression(paste(bolditalic("NITR"), bold(" gene family"))), expression(bolditalic("C1QL4")), expression(bolditalic("ELA1")), expression(bolditalic("RNASEL3")))) +
  scale_fill_manual(values = brewer.pal(n=8, name="Set1")[c(7,2,5,6,3,4,8,1)], limits = c("PIGG", "STAG3", "GPC2", "ITGA6", "NITR", "C1QL4", "ELA1", "RNASEL3"), labels = c(expression(bolditalic("PIGG")), expression(bolditalic("STAG3")), expression(bolditalic("GPC2")), expression(bolditalic("ITGA6")), expression(paste(bolditalic("NITR"), bold(" gene family"))), expression(bolditalic("C1QL4")), expression(bolditalic("ELA1")), expression(bolditalic("RNASEL3")))) +
  labs(x = "Chromosome 6 (Mb)", y = "Genes") +
  geom_hline(yintercept=0, linewidth=0.5) +
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.key.size = unit(.5, "lines"),
        legend.text = element_text(size = 13),
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
  guides(fill = guide_legend(nrow = 1), color = guide_legend(nrow = 1))


#------------------------------------------

diversity_in <- read.csv("diversity_AtlPac.csv", header=T)
diversity_in <- diversity_in[!is.na(diversity_in$fst_Atlantic.Pacific),]
diversity_in$fst_Atlantic.Pacific[diversity_in$fst_Atlantic.Pacific < 0] <- 0

diversity <- diversity_in %>% 
  mutate(status = ifelse(fst_Atlantic.Pacific < 0.15990814 & Atlantic > 0.007347820 & Pacific > 0.006718794, "Selected", "Unselected"))

diversity <- diversity %>% filter(chromosome == "chr6" & window_pos_1 > 5930000 & window_pos_2 <= 7070000)

p1 <-  ggplot() + 
  geom_rect(aes(xmin=((diversity$window_pos_1[diversity$status == "Selected"] [1]) + 2500)/1000000, xmax=((diversity$window_pos_2[diversity$status == "Selected"] [length(diversity$window_pos_2[diversity$status == "Selected"])]) + 2500)/1000000), ymin=0, ymax=Inf, fill = "red", alpha = .15) +
  geom_line(data = diversity, aes(x = (window_pos_1 + 2500)/1000000, y = fst_Atlantic.Pacific)) +
  labs(y = expression(bolditalic(F[ST]))) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.25))+ 
  scale_x_continuous(expand = c(0,0), limits = c(5.93, 7.07), breaks = seq(6,7,.2)) + 
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
  scale_x_continuous(expand = c(0,0), limits = c(5.93, 7.07), breaks = seq(6,7,.2)) + 
  scale_color_manual(values = c("dxy_Atlantic.Pacific" = "black", "Pacific" = "green3", "Atlantic" = "dodgerblue"), guide = guide_legend(reverse = TRUE), labels = c("dxy_Atlantic.Pacific"=expression(italic(d[xy])), "Pacific"=expression(italic(pi[Pacific])), "Atlantic"=expression(italic(pi[Atlantic])))) + 
  geom_hline(yintercept=mean(diversity_in$dxy_Atlantic.Pacific), linewidth=0.5, linetype = "dashed", color ="black", alpha = 1) +
  geom_hline(yintercept=mean(diversity_in$Pacific), linewidth=0.3, linetype = "dashed", color ="green3", alpha = 1) +
  geom_hline(yintercept=mean(diversity_in$Atlantic), linewidth=0.3, linetype = "dashed", color ="dodgerblue", alpha = 1) +
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

p <- plot_grid(p1, p2, q, ncol=1, rel_heights=c(1, 1, 1.3), align = "v") + theme(plot.margin = margin(2, 0, 0, 0))
png(filename="plot.png", width = 2000, height = 1000, units = 'px',res=195) 
print(p)
dev.off()
