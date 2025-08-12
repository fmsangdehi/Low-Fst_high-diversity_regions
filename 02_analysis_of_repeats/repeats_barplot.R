library(tidyverse)
library(cowplot)

diversity_dir <- "/path/to/diversity_AtlPac.csv"
repeat_prop_dir <- "Ch_v2.0.2.fasta.out.win_5000.superclasses.prop.csv"

diversity <- read.csv(diversity_dir, header=T)
diversity <- diversity[!is.na(diversity$fst_Atlantic.Pacific),]
diversity$fst_Atlantic.Pacific[diversity$fst_Atlantic.Pacific < 0] <- 0
diversity <- diversity %>%
  mutate(window_pos_1_bed = as.integer(plyr::round_any(window_pos_1, 5000, floor))) %>%
  mutate(index = paste(chromosome, window_pos_1_bed, sep = "_"))


df.super.p <- read.table(repeat_prop_dir, header = T)
df.super.p <- df.super.p %>%
  mutate(index = paste(CHROM, START, sep = "_"))

diversity_repeats <- merge(diversity, df.super.p, by="index",all.x=T)
diversity_repeats <- diversity_repeats %>%
  arrange(chromosome, window_pos_1)

diversity_repeats <- diversity_repeats %>% 
  mutate(other_rep_prop = apply(diversity_repeats[,c(19:30)], 1, sum)) %>% 
  mutate(total_rep_prop = apply(diversity_repeats[,c(14:30)], 1, sum))


diversity_repeats <- diversity_repeats %>%
  mutate(status = ifelse(fst_Atlantic.Pacific < 0.15990814 & Atlantic > 0.007347820 & Pacific > 0.006718794, "High diversity", "Background"))

diversity_repeats_long <- diversity_repeats %>%
  select(c(chromosome, window_pos_1, window_pos_2, status, Atlantic, Pacific, fst_Atlantic.Pacific, dxy_Atlantic.Pacific, Simple_repeat, DNA, LINE, LTR, Low_complexity,other_rep_prop, total_rep_prop)) %>%
  pivot_longer(cols = c(Atlantic, Pacific, fst_Atlantic.Pacific, dxy_Atlantic.Pacific, Simple_repeat, DNA, LINE, LTR, Low_complexity,other_rep_prop, total_rep_prop), values_to = "values", names_to = "class")

diversity_repeats_long$class <- factor(diversity_repeats_long$class, levels=c("fst_Atlantic.Pacific", "Atlantic", "Pacific", "dxy_Atlantic.Pacific", "Simple_repeat", "DNA", "LINE", "LTR", "Low_complexity","other_rep_prop", "total_rep_prop"), labels=c("Fst", "pi-Atlantic", "pi-Pacific", "dxy", "Simple repeat", "DNA transposon", "LINE", "LTR", "Low complexity","Other repeats", "Total repeats"))
diversity_repeats_long$status <- factor(diversity_repeats_long$status, levels=c("High diversity", "Background"), labels=c("High diversity", "Background"))

summary <- diversity_repeats_long %>%
  group_by(class, status, .add = TRUE) %>% #or: group_by(groups = class + status) %>%
  summarise(mean = mean(values, na.rm = T), sd = sd(values, na.rm = T))


q2 <- summary %>% 
  filter(class %in% c("Simple repeat", "DNA transposon", "LINE", "LTR", "Low complexity","Other repeats")) %>%
  ggplot(aes(x = class, y = mean, fill = status, color = status)) +
  geom_bar(stat="identity", position=position_dodge(width = 0.8), width=.7, alpha = 0.3) +
#  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, fill = status, color = status), width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values=c("red2", "black")) +
  scale_color_manual(values=c("red2", "black")) +
  scale_y_continuous(limits = c(0, 0.091)) +
  labs(y = "Proportion", x = "Repetitive elements") +
  theme_bw() +
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(legend.title = element_blank(),
        legend.key.size = unit(.9, "lines"),
        legend.text = element_text(size = 13, face="bold"),
        legend.text.align = 0,
        legend.margin = margin(t = 0, unit='cm'),
        legend.direction = "vertical",
        legend.position = c(0.80, 0.90),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16, face="bold"),
        axis.text.x = element_text(angle = 55, size = 16, face = "bold", color = "black", hjust=1),
        axis.text.y = element_text(size=13),
  )

png(filename="repeats_barplot.png", width = 900, height = 900, units = 'px',res=170) 
print(q2)
dev.off()
