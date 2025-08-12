library(tidyverse)

df.super.p <- read.table("Ch_v2.0.2.fasta.out.win_5000.superclasses.prop.csv", header = T)
df.super.n <- read.table("Ch_v2.0.2.fasta.out.win_5000.superclasses.n.csv", header = T)

#switch data to long form
df.super.p.long <- df.super.p %>%
  pivot_longer(cols = -c("CHROM", "START", "GLOBAL_START","LENGTH"), values_to = "proportion", names_to = "class")
df.super.n.long <- df.super.n %>%
  pivot_longer(cols = -c("CHROM", "START", "GLOBAL_START","LENGTH"), values_to = "n", names_to = "class")

#add n to prop
df.super.p.long$n <- df.super.n.long$n

#find out what are the most common classes of repeats
df.super.p.long %>% group_by(class) %>% 
  summarise(mean.prop = mean(proportion),
            sum.n = sum(n)) %>% 
  arrange(desc(mean.prop))

#Based on proportion of repeats: Simple_repeat, DNA, LINE, LTR and Low_complexity
#Based on repeat counts: Simple_repeat, DNA, Low_complexity, LTR and LINE

#---------------------------------------------------------------------

##### Correlation between proportion and number of repeats

df.super.p <- read.table("Ch_v2.0.2.fasta.out.win_5000.superclasses.prop.csv", header = T)
df.super.n <- read.table("Ch_v2.0.2.fasta.out.win_5000.superclasses.n.csv", header = T)

df.super.p <- df.super.p %>% 
  mutate(total = apply(df.super.p[,5:21], 1, sum))
df.super.n <- df.super.n %>% 
  mutate(total = apply(df.super.n[,5:21], 1, sum))

#switch data to long form
df.super.p.long <- df.super.p %>%
  pivot_longer(cols = c(6,5,7,8,9,22), values_to = "proportion", names_to = "class")
df.super.n.long <- df.super.n %>%
  pivot_longer(cols = c(6,5,7,8,9,22), values_to = "n", names_to = "class")

#add n to prop
df.super.p.long$n <- df.super.n.long$n

df.super.p.long <- df.super.p.long %>%
  as.data.frame() %>%
  mutate(class = factor(class, levels=c("Simple_repeat", "DNA", "LINE","LTR","Low_complexity","total"), labels = c("Simple repeat", "DNA", "LINE","LTR","Low complexity","total")))

# Calculate correlation for each class
cors <- df.super.p.long %>% group_by(class) %>% 
    summarise(cor = round(cor(proportion, n, use="complete.obs"), 2))

cors_spearman <- df.super.p.long %>% group_by(class) %>% 
  summarise(cor = round(cor(proportion, n, use="complete.obs", method = "spearman"), 2))


p <- ggplot(df.super.p.long, aes(proportion, n)) + 
     geom_point(size=.5, alpha=.3, color="darkcyan") + 
     facet_wrap(~class, ncol = 3, scales = "free") + 
     geom_text(data=cors, aes(label=paste("r=", cor, sep=""), x = Inf, y = Inf), color="black", size=5, hjust = 1, vjust = 1) + 
     theme_bw() + 
     theme(legend.position = "none",
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           axis.title.x = element_blank(),
           axis.title.y = element_blank(),
)

png(filename="corrPlots_repeats_prop-vs-n.png", width = 1200, height = 800, units = 'px', res=170) 
print(p)
dev.off()
