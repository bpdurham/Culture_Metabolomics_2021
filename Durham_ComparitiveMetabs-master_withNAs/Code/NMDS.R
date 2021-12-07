source('SourceCode/biostats.R')
library("vegan") 
library("cluster")
library("pvclust")
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(here)
library(RColorBrewer)

#Name inputs-----
#Read in PA dat, mudge to get into a matrix, toss any compounds that weren't seen ever, make NAs 0s ----
PA.dat <- cbind(rownames(dat6), data.frame(dat6))
colnames(PA.dat)[1] = "Organism"
meta.dat <- PA.dat %>% dplyr::select(Organism, Domain, Group)

matrix_noNA = matrix
matrix_noNA[is.na(matrix_noNA)] <- 0

#Run NMDS, extract point location
nmds.raw<-metaMDS(matrix_noNA, distance='jaccard', k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)
monte.pvalue.raw <-nmds.monte(t(matrix_noNA), distance='jaccard', k=2, autotransform=FALSE, trymax=20)
monte.pvalue.result.raw <- monte.pvalue.raw[[2]]
print(paste(monte.pvalue.result.raw, "= pvalue of nmds"))
pointlocation.nmds.raw <- nmds.raw[['points']] %>% as.data.frame() %>%
  mutate(Organism = rownames(nmds.raw[['points']])) %>%
  left_join(meta.dat, by = "Organism") 


#Plot out the point location for the raw NMDS----
pal = c("gold", "purple", "red", "darkorange", "deepskyblue3", "navy", "green2")
d.raw <- ggplot(data = pointlocation.nmds.raw, aes(x =MDS1, y =  MDS2, group = Group,
                                              colour = Group,
                                              fill = Group))+
  geom_point(position = position_jitter(width = 0.02, height = 0.02), size = 2) +
  geom_text(stat='identity', position='identity', label=pointlocation.nmds.raw$Organism, size=3) +
  scale_fill_manual(values = pal)+
  scale_color_manual(values = pal)+
  annotate("text", x = min(pointlocation.nmds.raw$MDS1)*.9, y = max(pointlocation.nmds.raw$MDS2)*.9, 
           label = paste0("Stress = ", 
                          round(nmds.raw[['stress']], digits = 2), 
                          "\n p = ", 
                          round(monte.pvalue.result.raw, digits = 4)), size = 3)+
    theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

d.raw  
print(d.raw)
save_plot("NMDS_PA.pdf", d.raw, base_height = 6, base_width = 7)
