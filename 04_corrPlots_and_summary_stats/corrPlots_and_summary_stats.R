library(plyr)
library(tidyverse)

diversity_dir <- "/path/to/diversity_AtlPac.csv"
repeats_prop_dir <- "/path/to/Ch_v2.0.2.fasta.out.win_5000.superclasses.prop.csv"

diversity <- read.csv(diversity_dir, header=T)
diversity <- diversity[!is.na(diversity$fst_Atlantic.Pacific),]
diversity$fst_Atlantic.Pacific[diversity$fst_Atlantic.Pacific < 0] <- 0
diversity <- diversity %>%
     mutate(window_pos_1_bed = as.integer(plyr::round_any(window_pos_1, 5000, floor))) %>%
     mutate(index = paste(chromosome, window_pos_1_bed, sep = "_"))

# there are some incompatibilities between plyr and dplyr. dplyr behaves differently if plyr is loaded. but there are some solutions:
# 1)[use detaching/reloading dplyr/plyr in each block as a workaround]: detach("package:plyr", unload=TRUE)
# 2)[if you load plyr first, and then dplyr, the masking will work the way you want, in that you'll still be able to use plyr's unique functions, but dplyr won't break.]
# If you need functions from both plyr and dplyr, load plyr first, then dplyr: library(plyr); library(dplyr)
# 3)[attach dplyr with library(dplyr) then use "::" to refer to functions from plyr]

df.super.p <- read.table(repeats_prop_dir, header = T)
df.super.p <- df.super.p %>%
     mutate(index = paste(CHROM, START, sep = "_"))

diversity_repeats <- merge(diversity, df.super.p, by="index",all.x=T)
diversity_repeats <- diversity_repeats %>%
     arrange(chromosome, window_pos_1)

#---------------------------------------------------------------------

##### mean values of diversity parameters and proportion of repeats in low-Fst regions and genome background
##### whole genome (in percent)
diversity_repeats <- diversity_repeats %>% 
     mutate(other_rep_prop = apply(diversity_repeats[,c(19:30)], 1, sum)) %>% 
     mutate(total_rep_prop = apply(diversity_repeats[,c(14:30)], 1, sum))

diversity_repeats_long <- diversity_repeats %>%
  pivot_longer(cols = c("fst_Atlantic.Pacific", "dxy_Atlantic.Pacific", "Atlantic", "Pacific", "Simple_repeat", "DNA", "LINE","LTR","Low_complexity","other_rep_prop","total_rep_prop"), values_to = "values", names_to = "class")

diversity_repeats_long %>% group_by(class) %>% 
  summarise(average = round(mean(values, na.rm=T)*100, 2))


##### low-Fst regions (in percent)
diversity_repeats_low_Fst <- diversity_repeats %>%
           filter(fst_Atlantic.Pacific < 0.15990814 & Atlantic > 0.007347820 & Pacific > 0.006718794)

diversity_repeats_low_Fst_long <- diversity_repeats_low_Fst %>%
  pivot_longer(cols = c("fst_Atlantic.Pacific", "dxy_Atlantic.Pacific", "Atlantic", "Pacific", "Simple_repeat", "DNA", "LINE","LTR","Low_complexity","other_rep_prop","total_rep_prop"), values_to = "values", names_to = "class")

diversity_repeats_low_Fst_long %>% group_by(class) %>% 
  summarise(average = round(mean(values, na.rm=T)*100, 2))


##### genome_background (in percent)
diversity_repeats_background <- diversity_repeats %>%
  filter(!(fst_Atlantic.Pacific < 0.15990814 & Atlantic > 0.007347820 & Pacific > 0.006718794))

diversity_repeats_background_long <- diversity_repeats_background %>%
  pivot_longer(cols = c("fst_Atlantic.Pacific", "dxy_Atlantic.Pacific", "Atlantic", "Pacific", "Simple_repeat", "DNA", "LINE","LTR","Low_complexity","other_rep_prop","total_rep_prop"), values_to = "values", names_to = "class")

diversity_repeats_background_long %>% group_by(class) %>% 
  summarise(average = round(mean(values, na.rm=T)*100, 2))


#---------------------------------------------------------------------
##### Mean comparison tests
diversity_repeats$status[diversity_repeats$fst_Atlantic.Pacific < 0.15990814 & diversity_repeats$Atlantic > 0.007347820 & diversity_repeats$Pacific > 0.006718794] <- "High diversity"
diversity_repeats$status[is.na(diversity_repeats$status)] <- "Background"

t.test(diversity_repeats[diversity_repeats$status == "High diversity", "fst_Atlantic.Pacific"], y = diversity_repeats[diversity_repeats$status == "Background", "fst_Atlantic.Pacific"], paired = FALSE, var.equal = FALSE)
t.test(diversity_repeats[diversity_repeats$status == "High diversity", "Atlantic"], y = diversity_repeats[diversity_repeats$status == "Background", "Atlantic"], paired = FALSE, var.equal = FALSE)
t.test(diversity_repeats[diversity_repeats$status == "High diversity", "Pacific"], y = diversity_repeats[diversity_repeats$status == "Background", "Pacific"], paired = FALSE, var.equal = FALSE)
t.test(diversity_repeats[diversity_repeats$status == "High diversity", "dxy_Atlantic.Pacific"], y = diversity_repeats[diversity_repeats$status == "Background", "dxy_Atlantic.Pacific"], paired = FALSE, var.equal = FALSE)

t.test(diversity_repeats[diversity_repeats$status == "High diversity", "Simple_repeat"], y = diversity_repeats[diversity_repeats$status == "Background", "Simple_repeat"], paired = FALSE, var.equal = FALSE)
t.test(diversity_repeats[diversity_repeats$status == "High diversity", "DNA"], y = diversity_repeats[diversity_repeats$status == "Background", "DNA"], paired = FALSE, var.equal = FALSE)
t.test(diversity_repeats[diversity_repeats$status == "High diversity", "LINE"], y = diversity_repeats[diversity_repeats$status == "Background", "LINE"], paired = FALSE, var.equal = FALSE)
t.test(diversity_repeats[diversity_repeats$status == "High diversity", "LTR"], y = diversity_repeats[diversity_repeats$status == "Background", "LTR"], paired = FALSE, var.equal = FALSE)
t.test(diversity_repeats[diversity_repeats$status == "High diversity", "Low_complexity"], y = diversity_repeats[diversity_repeats$status == "Background", "Low_complexity"], paired = FALSE, var.equal = FALSE)
t.test(diversity_repeats[diversity_repeats$status == "High diversity", "other_rep_prop"], y = diversity_repeats[diversity_repeats$status == "Background", "other_rep_prop"], paired = FALSE, var.equal = FALSE)
t.test(diversity_repeats[diversity_repeats$status == "High diversity", "total_rep_prop"], y = diversity_repeats[diversity_repeats$status == "Background", "total_rep_prop"], paired = FALSE, var.equal = FALSE)

#---------------------------------------------------------------------

##### correlation coefficients and correlation matrix for diversity statistics
##### in the following I use ggpairs() function from package GGally

library(GGally)

diversity_repeats$status[diversity_repeats$fst_Atlantic.Pacific < 0.15990814 & diversity_repeats$Atlantic > 0.007347820 & diversity_repeats$Pacific > 0.006718794] <- "High diversity"
diversity_repeats$status[is.na(diversity_repeats$status)] <- "Background"

diversity_repeats$status <- factor(diversity_repeats$status, levels=c("High diversity", "Background"), labels=c("High diversity", "Background"))

diversity_repeats <- diversity_repeats %>% 
    arrange(chromosome, window_pos_1)

colnames(diversity_repeats)[c(7, 5, 6, 8)] <- c(expression(bolditalic(F[ST])), expression(bolditalic(pi[Atlantic])), expression(bolditalic(pi[Pacific])), expression(bolditalic(d[xy])))

ggally_mysmooth <- function(data, mapping, ...){
  ggplot(data = data, mapping=mapping) +
   geom_density(mapping = aes(fill=status, colour = status), alpha = .3) #with this code, the height of the two density curves in each block is not proportional to the number of variants in each group
}

p <- ggpairs(diversity_repeats, c(7, 8, 5, 6),
     mapping = ggplot2::aes(color = status),
     labeller = "label_parsed",
     lower = list(continuous = wrap("points", alpha = 0.3, shape = 20, size = .7)),
     diag = list(continuous = ggally_mysmooth) #this is to solve the above problem
)
p <- p + theme_bw()
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               strip.text = element_text(margin = margin(t = 1.2, r = 1.2, b = 1.2, l = 1.2, unit = "pt"), size = 12, color = "black")
)

for(i in 1:p$nrow) {
  for(j in 1:p$ncol){
    p[i,j] <- p[i,j] + 
      scale_fill_manual(values=c("red2", "black")) +
      scale_color_manual(values=c("red2", "black"))
  }
}

png(filename="corrMatrix_diversity_method_1.png", width = 10000, height = 10000, units = 'px', res=1200) 
print(p)
dev.off()

#---------------------------------------------------------------------

##### However, if I want to add the peak values to the plots, the above method doesn't work. Therefore I used the following method.

library(GGally)

diversity_repeats$status[diversity_repeats$fst_Atlantic.Pacific < 0.15990814 & diversity_repeats$Atlantic > 0.007347820 & diversity_repeats$Pacific > 0.006718794] <- "High diversity"
diversity_repeats$status[is.na(diversity_repeats$status)] <- "Background"

diversity_repeats$status <- factor(diversity_repeats$status, levels=c("High diversity", "Background"), labels=c("High diversity", "Background"))

diversity_repeats <- diversity_repeats %>% 
     arrange(chromosome, window_pos_1)

options(scipen = 999)

peak_x_low <- list()
for (i in c(8, 5, 6)){
j <- match(i, c(8, 5, 6))
peak_y_index <- which.max(density(diversity_repeats[diversity_repeats$status == "High diversity", i], na.rm = T)$y)
peak_x_low[[j]] <- sprintf("%.4f", round(density(diversity_repeats[diversity_repeats$status == "High diversity", i], na.rm = T)$x[peak_y_index], 4))
}

peak_x_normal <- list()
for (i in c(8, 5, 6)){
j <- match(i, c(8, 5, 6))
peak_y_index <- which.max(density(diversity_repeats[diversity_repeats$status == "Background", i], na.rm = T)$y)
peak_x_normal[[j]] <- sprintf("%.4f", round(density(diversity_repeats[diversity_repeats$status == "Background", i], na.rm = T)$x[peak_y_index], 4))
}


plot1 <- ggplot(data = diversity_repeats, mapping = aes(x = fst_Atlantic.Pacific)) +
    geom_density(mapping = aes(fill=status, colour = status), alpha = .3) #not proportionate

plot2 <- ggplot(data = diversity_repeats, mapping = aes(x = dxy_Atlantic.Pacific)) +
    geom_density(mapping = aes(fill=status, colour = status), alpha = .3) + #not proportionate
    annotate(geom="text", x = Inf, y = Inf, label=paste("", paste("peak=", peak_x_normal[1], " ", sep=""), sep="\n"), color="black", size=4, hjust = 1.5, vjust = 1) +
    annotate(geom="text", x = Inf, y = Inf, label=paste("", "", paste("peak=", peak_x_low[1], " ", sep=""), sep="\n"), color="red2", size=4, hjust = 1.5, vjust = 1)

plot3 <- ggplot(data = diversity_repeats, mapping = aes(x = Atlantic)) +
    geom_density(mapping = aes(fill=status, colour = status), alpha = .3) + #not proportionate
    annotate(geom="text", x = Inf, y = Inf, label=paste("", paste("peak=", peak_x_normal[2], " ", sep=""), sep="\n"), color="black", size=4, hjust = 1.6, vjust = 1) +
    annotate(geom="text", x = Inf, y = Inf, label=paste("", "", paste("peak=", peak_x_low[2], " ", sep=""), sep="\n"), color="red2", size=4, hjust = 1.6, vjust = 1)

plot4 <- ggplot(data = diversity_repeats, mapping = aes(x = Pacific)) +
    geom_density(mapping = aes(fill=status, colour = status), alpha = .3) + #not proportionate
    annotate(geom="text", x = Inf, y = Inf, label=paste("", paste("peak=", peak_x_normal[3], " ", sep=""), sep="\n"), color="black", size=4, hjust = 1.6, vjust = 1) +
    annotate(geom="text", x = Inf, y = Inf, label=paste("", "", paste("peak=", peak_x_low[3], " ", sep=""), sep="\n"), color="red2", size=4, hjust = 1.6, vjust = 1)

colnames(diversity_repeats)[c(7, 5, 6, 8)] <- c(expression(bolditalic(F[ST])), expression(bolditalic(pi[Atlantic])), expression(bolditalic(pi[Pacific])), expression(bolditalic(d[xy])))

p <- ggpairs(diversity_repeats, c(7, 8, 5, 6),
     mapping = ggplot2::aes(color = status),
     labeller = "label_parsed",
     lower = list(continuous = wrap("points", alpha = 0.3, shape = 20, size = .7)),
     diag = "blank",
)

p[1, 1] <- plot1
p[2, 2] <- plot2
p[3, 3] <- plot3
p[4, 4] <- plot4

p <- p + theme_bw()
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               strip.text = element_text(margin = margin(t = 1.2, r = 1.2, b = 1.2, l = 1.2, unit = "pt"), size = 12, color = "black")
)

for(i in 1:p$nrow) {
  for(j in 1:p$ncol){
    p[i,j] <- p[i,j] + 
      scale_fill_manual(values=c("red2", "black")) +
      scale_color_manual(values=c("red2", "black"))
  }
}

png(filename="corrMatrix_diversity_method_2.png", width = 10000, height = 10000, units = 'px', res=1200) 
print(p)
dev.off()
