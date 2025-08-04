library(tidyverse)
library(cowplot)

diversity <- read.csv("diversity_AtlPac.csv", header=T)
diversity <- diversity[!is.na(diversity$fst_Atlantic.Pacific),]
diversity$fst_Atlantic.Pacific[diversity$fst_Atlantic.Pacific < 0] <- 0

diversity <- diversity %>%
  mutate(status = ifelse(fst_Atlantic.Pacific < 0.15990814 & Atlantic > 0.007347820 & Pacific > 0.006718794, "High diversity", "Background"))

diversity_long <- diversity %>%
  pivot_longer(cols = c("fst_Atlantic.Pacific", "Atlantic", "Pacific", "dxy_Atlantic.Pacific"), values_to = "values", names_to = "class")

diversity_long$class <- factor(diversity_long$class, levels=c("fst_Atlantic.Pacific", "Atlantic", "Pacific", "dxy_Atlantic.Pacific"), labels=c(expression(bolditalic(F[ST])), expression(bolditalic(pi[Atlantic])), expression(bolditalic(pi[Pacific])), expression(bolditalic(d[xy]))))
diversity_long$status <- factor(diversity_long$status, levels=c("High diversity", "Background"), labels=c("High diversity", "Background"))


p2 <- diversity_long %>% 
  filter(class == "bolditalic(F[ST])") %>% 
  ggplot(aes(x = class, y = values)) +
  geom_boxplot(aes(fill = status, color = status), position=position_dodge(.7), alpha = 0.3, outlier.size=1, outlier.shape=21, width = 0.5) +
  scale_fill_manual(values=c("red2", "black")) +
  scale_color_manual(values=c("red2", "black")) +
  scale_x_discrete(labels = ggplot2:::parse_safe) +
  labs(y = "Differentiation") +
#  facet_wrap(~class) +
  stat_summary(fun=mean, geom="point", aes(group=status), position=position_dodge(.7), color="#ECEC86", size=2) +
  theme_bw() +
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16, color = "black", face="bold"),
        axis.text.x = element_text(angle = 45, color = "black", size = 18, face = "bold", hjust = 1),
        axis.text.y = element_text(size=13),
        plot.margin = margin(t = 5, r=5, b=5, l=5, unit="pt")
  )

q2 <- diversity_long %>% 
  filter(class != "bolditalic(F[ST])") %>% 
  ggplot(aes(x = class, y = values)) +
  geom_boxplot(aes(fill = status, color = status), position=position_dodge(.7), alpha = 0.3, outlier.size=1, outlier.shape=21, width = 0.5) +
  scale_fill_manual(values=c("red2", "black")) +
  scale_color_manual(values=c("red2", "black")) +
  scale_x_discrete(labels = ggplot2:::parse_safe) +
  scale_y_continuous(breaks=c(0, 0.01, 0.02, 0.03, 0.04, 0.05)) +
  labs(y = "Nucleotide diversity") +
#  facet_wrap(~class) +
  stat_summary(fun=mean, geom="point", aes(group=status), position=position_dodge(.7), color="#ECEC86", size=2) +
  theme_bw() +
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(legend.title = element_blank(),
        legend.key.size = unit(0.9, "lines"),
        legend.text = element_text(size = 13, face="bold"),
        legend.text.align = 0,
        legend.margin = margin(t = 0, unit='cm'),
        legend.direction = "vertical", legend.position = c(0.75, 0.95),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16, color = "black", face="bold"),
        axis.text.x = element_text(angle = 45, color = "black", size = 18, face = "bold", hjust = 1),
        axis.text.y = element_text(size=13),
        plot.margin = margin(t = 5, r=5, b=5, l=10, unit="pt")
  )

plot2 <- plot_grid(p2, q2, ncol=2, rel_widths=c(1.5, 3), align = "h")
png(filename="diversity_boxplot.png", width = 900, height = 900, units = 'px',res=170) 
print(plot2)
dev.off()
