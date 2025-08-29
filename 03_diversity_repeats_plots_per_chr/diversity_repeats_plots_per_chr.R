library(tidyverse)
library(cowplot)

diversity_dir <- "/path/to/diversity_AtlPac.csv"
repeats_prop_dir <- "/path/to/Ch_v2.0.2.fasta.out.win_5000.superclasses.prop.csv"

diversity <- read.csv(diversity_dir, header=T)
diversity <- diversity[!is.na(diversity$fst_Atlantic.Pacific),]
diversity$fst_Atlantic.Pacific[diversity$fst_Atlantic.Pacific < 0] <- 0
diversity <- diversity %>%
     mutate(window_pos_1_bed = as.integer(plyr::round_any(window_pos_1, 5000, floor))) %>%
     mutate(index = paste(chromosome, window_pos_1_bed, sep = "_"))


df.super.p <- read.table(repeats_prop_dir, header = T)
df.super.p <- df.super.p %>%
     mutate(index = paste(CHROM, START, sep = "_"))

diversity_repeats <- merge(diversity, df.super.p, by="index",all.x=T)
diversity_repeats <- diversity_repeats %>%
     arrange(chromosome, window_pos_1)

#switch repeats data to long form
repeats <- diversity_repeats[, 10:30]
repeats_long <- repeats %>%
     pivot_longer(cols = -c("CHROM", "START", "GLOBAL_START","LENGTH"), values_to = "proportion", names_to = "class")

#find out what are the most common classes of repeats
repeats_long %>%
     group_by(class) %>% 
     summarise(mean.prop = mean(proportion, na.rm=T)) %>% 
     arrange(mean.prop)

#mostly it is Simple_repeat, DNA, LINE, LTR, and Low_complexity
repeats_long <- repeats_long %>% 
  mutate(simplified.class = ifelse(class != "Simple_repeat" & class!="DNA" & class!= "LINE" & class!="LTR" & class!="Low_complexity", "other", class))


repeats_long$simplified.class <- factor(repeats_long$simplified.class, levels=c("other", "Low_complexity", "LTR", "LINE", "DNA", "Simple_repeat"), labels=c("Other repeats", "Low complexity", "LTR", "LINE", "DNA transposon", "Simple repeat"))

#finding y-axis limits for proportion of repeats
repeats$sum_prop <- apply(repeats[,c(5:21)], 1, sum)
head(sort(repeats$sum_prop), 10)
tail(sort(repeats$sum_prop), 10)
tail(repeats[order(repeats$sum_prop),], 6)

#---------------------------------------------------------------------

#plots of diversity statistics and proportion of repeats per chromosome

diversity <- diversity %>% 
  mutate(status = ifelse(fst_Atlantic.Pacific < 0.15990814 & Atlantic > 0.007347820 & Pacific > 0.006718794, "Low Fst", "Background"))

plots <- list() 

pdf("plots_per_chr.pdf", width = 8.5, height = 11)
for (chr in paste("chr", 1:26, sep="")){
  
p1 <- diversity %>% filter(chromosome == chr) %>% 
    ggplot(aes(x = window_pos_1/1000000, y = fst_Atlantic.Pacific, color = status)) + 
    geom_point(size = .8,shape = 20,stroke = 0.2)+
    labs(y = expression(bolditalic(F[ST]))) +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25))+ 
    scale_x_continuous(expand = c(.02, .02)) + 
    scale_color_manual(values = c("Low Fst" = "red2", "Background" = "black")) + 
    theme_bw() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=16, face="bold"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.ticks.y = element_line(size=0.2))

p2 <- diversity %>% filter(chromosome == chr) %>% 
    ggplot(aes(x = window_pos_1/1000000, y = dxy_Atlantic.Pacific, color = status)) + 
    geom_point(size = .8,shape = 20,stroke = 0.2)+
    labs(y = expression(bolditalic(d[xy]))) +
    scale_y_continuous(limits = c(0,0.056))+ 
    scale_x_continuous(expand = c(.02, .02)) + 
    scale_color_manual(values = c("Low Fst" = "red2", "Background" = "black")) + 
    theme_bw() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=16, face="bold"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.ticks.y = element_line(size=0.2))

p3 <- diversity %>% filter(chromosome == chr) %>% 
    ggplot(aes(x = window_pos_1/1000000, y = Atlantic, color = status)) + 
    geom_point(size = .8,shape = 20,stroke = 0.2)+
    labs(y = expression(bolditalic(pi[Atlantic]))) +
    scale_y_continuous(limits = c(0,0.056))+ 
    scale_x_continuous(expand = c(.02, .02)) + 
    scale_color_manual(values = c("Low Fst" = "red2", "Background" = "black")) + 
    theme_bw() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=16, face="bold"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.ticks.y = element_line(size=0.2))

p4 <- diversity %>% filter(chromosome == chr) %>% 
    ggplot(aes(x = window_pos_1/1000000, y = Pacific, color = status)) + 
    geom_point(size = .8,shape = 20,stroke = 0.2)+
    labs(y = expression(bolditalic(pi[Pacific]))) +
    scale_y_continuous(limits = c(0,0.056))+ 
    scale_x_continuous(expand = c(.02, .02)) + 
    scale_color_manual(values = c("Low Fst" = "red2", "Background" = "black")) + 
    theme_bw() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=16, face="bold"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.ticks.y = element_line(size=0.2))

q1 <- repeats_long %>% filter(CHROM == chr) %>% 
    ggplot(aes(x = START/ 1000000, y = proportion, fill = simplified.class)) + 
    geom_area() + 
    labs(x = "Position (Mb)", y = "Repeats") +
    scale_y_continuous(limits = c(0, 2.11)) +
    scale_x_continuous(expand = c(.02, .02)) + 
    theme_bw() + 
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.title.x = element_text(size=16, face="bold"),
          axis.title.y = element_text(size=15, face="bold"),
          axis.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12),
          axis.ticks.y = element_line(size=0.2)) +
    theme(legend.position = "bottom",
          legend.key.size = unit(.6, "lines"),
          legend.text = element_text(size = 13, face="bold"),
          legend.margin = margin(t = 0, unit='cm')) +
    scale_fill_discrete(limits=c("Simple repeat", "DNA transposon", "LINE", "LTR", "Low complexity", "Other repeats")) +
    guides(fill = guide_legend(nrow = 1, title=NULL, override.aes = list(size = .5)))


chr_num <- str_sub(chr, start= 4)
title <- ggdraw() + draw_label(paste("Chromosome", chr_num, paste=" "), fontface = 'bold', size=18, vjust = .5, hjust = 0, x = 0)

plots[[chr]] <- plot_grid(title, p1, p3, p4, p2, q1, ncol=1, rel_heights=c(0.1, 1, 1, 1, 1, 1.3), align = "v") + theme(plot.margin = margin(7, 5, 7, 3))
print(plots[[chr]])
}
dev.off()
