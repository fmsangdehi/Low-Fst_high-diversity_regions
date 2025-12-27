library(tidyverse)

##### Prepare RR data frame
RR <- read.table("Baltic.LDhat.100K.txt", header = T, sep ="\t")
RR_5kb <- rep(RR$RR, each = 20)

list <- list()

for(i in 1:26){
  list[[i]] <- data.frame(chr = paste0("chr", i), bin_5kb = seq(1, nrow(RR[RR$CHR == i,])*100000 - 4999, by = 5000))
}

bin_5kb <- do.call(rbind, list)
RR_new <- cbind(bin_5kb, RR_5kb) %>%
  mutate(index = paste(chr, bin_5kb, sep = "_")) %>%
  select(index, RR_5kb)


##### Prepare stats data frame
stats <- read.csv("diversity_AtlPac.csv", header = T) %>%
  mutate(fst_Atlantic.Pacific = replace(fst_Atlantic.Pacific, fst_Atlantic.Pacific < 0, 0)) %>%
  subset(!is.na(fst_Atlantic.Pacific)) %>%
  mutate(status = ifelse(fst_Atlantic.Pacific < 0.15990814 & Atlantic > 0.007347820 & Pacific > 0.006718794, "Low_Fst", "Background")) %>%
  mutate(window_pos_1 = as.integer(plyr::round_any(window_pos_1, 5000, floor)) + 1) %>%
  mutate(index = paste(chromosome, window_pos_1, sep = "_"))

##### Merge RR data frame with stats data frame
stats_RR <- inner_join(stats, RR_new, by = join_by("index")) %>%
  select(-index)
any(is.na(stats_RR$RR_5kb))
any(is.na(stats_RR$fst_Atlantic.Pacific))
sum(as.integer(is.na(stats_RR$fst_Atlantic.Pacific)))
nrow(stats_RR[is.na(stats_RR$fst_Atlantic.Pacific) | is.na(stats_RR$dxy_Atlantic.Pacific) | is.na(stats_RR$Atlantic) | is.na(stats_RR$Pacific),])

##### calculate correlation
stats_RR_long <- stats_RR %>%
  pivot_longer(c(Atlantic, Pacific, fst_Atlantic.Pacific, dxy_Atlantic.Pacific), names_to = "stat", values_to = "stat_values")

cors <- stats_RR_long %>% group_by(stat) %>%
  summarise(cor = round(cor(RR_5kb, stat_values, use = "complete.obs", method = "pearson"), 2))

## Another way for doing this:
plyr::ddply(stats_RR_long, "stat", summarise, cor = round(cor(stat_values, RR_5kb, use="complete.obs"), 2))

##### Test the significance
cor.test(stats_RR$RR_5kb, stats_RR$Atlantic, method = "pearson")
cor.test(stats_RR$RR_5kb, stats_RR$Pacific, method = "pearson")
cor.test(stats_RR$RR_5kb, stats_RR$dxy_Atlantic.Pacific, method = "pearson")
cor.test(stats_RR$RR_5kb, stats_RR$fst_Atlantic.Pacific, method = "pearson")

##### scatter plots
stats_RR_long$stat <- factor(stats_RR_long$stat, levels=c("fst_Atlantic.Pacific", "dxy_Atlantic.Pacific", "Atlantic", "Pacific"), labels = c("Fst", "Dxy", "pi_Atlantic", "pi_Pacific"))
cors$stat <- factor(cors$stat, levels=c("fst_Atlantic.Pacific", "dxy_Atlantic.Pacific", "Atlantic", "Pacific"), labels = c("Fst", "Dxy", "pi_Atlantic", "pi_Pacific"))

list_2 <- list()

for(i in levels(stats_RR_long$stat)){
  list_2[[i]] <- ggplot(stats_RR_long[stats_RR_long$stat == i,], aes(RR_5kb, stat_values)) + 
    geom_point(aes(color=status), size=.3, alpha=.15) + 
    geom_smooth(method='lm', se = T) + 
    geom_smooth(method='loess', se = F, linetype = "dashed", colour = "orange1") + 
    scale_color_manual(values = c(Low_Fst = "red2", Background = "gray20")) + 
    annotate(geom="text", x = Inf, y = Inf, label=paste0("r=", cors[cors$stat == i, "cor"]), color="black", size=4, fontface = "bold", hjust = 1.1, vjust = 1.4) +
    labs(y = i) +
    theme_bw() + 
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12, color = "black", face = "bold")
    )
}

p <- cowplot::plot_grid(list_2[[1]], list_2[[2]], list_2[[3]], list_2[[4]],ncol=2, rel_heights=c(1, 1, 1, 1), align = "v") + theme(plot.margin = margin(7, 5, 7, 3))

x_axis_title <- cowplot::ggdraw() + cowplot::draw_label("Recombination rate", fontface = 'bold', size=15, vjust = .5, hjust = .5, x = .5)

p_final <- cowplot::plot_grid(p, x_axis_title, ncol=1, rel_heights=c(1, .05), align = "v")

png(filename="corrPlots_RR_vs_diversity.png", width = 2000, height = 2000, units = 'px', res=300) 
print(p_final)
dev.off()
