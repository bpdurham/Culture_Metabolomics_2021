#modified from original code from Katherine Heal kheal.github.io
#This analysis attempts to use the variables to be able to predict an assignment to a group.  Right now its just in and out of group (set it as test.group below).  There's a few things you'd want to report along the way for good measure, but the analysis works really well.  To see which compounds are best at discriminating an organism within (or not within), look at the coeff.vars.2 data frame paired with the final plot.  The compounds with the largest loadings (either pos or neg, check plot to see which) are the best at saying whether an organism is in a group

source('SourceCode/biostats.R')
library("vegan") 
library("cluster")
library("pvclust")
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(here)
library(RColorBrewer)
library(MASS)
library(gplots)
library(ggpubr)


#Name inputs-----
#Read in PA dat, mudge to get into a matrix, toss any compounds that weren't seen ever, make NAs 0s ----
PA.dat <- cbind(rownames(dat6), data.frame(dat6))
colnames(PA.dat)[1] = "Organism"
meta.dat <- PA.dat %>% dplyr::select(Organism, Domain, Group) 
matrix_trim_noNA = matrix_trim
matrix_trim_noNA[is.na(matrix_trim_noNA)] <- 0
wide.matrix = matrix_trim_noNA

#remove universal metabolites
wide.matrix = matrix_trim_noNA[, -which(colnames(matrix_trim_noNA) == "cAMP-HILIC" | 
                                        colnames(matrix_trim_noNA) == "NADH-HILIC" | 
                                        colnames(matrix_trim_noNA) == "Thymidine-HILIC" | 
                                        colnames(matrix_trim_noNA) == "Uracil-HILIC" | 
                                        colnames(matrix_trim_noNA) == "cGMP-HILIC" | 
                                        colnames(matrix_trim_noNA) == "GMP-HILIC" | 
                                        colnames(matrix_trim_noNA) == "ADP-HILIC" | 
                                        colnames(matrix_trim_noNA) == "NADP-HILIC" | 
                                        colnames(matrix_trim_noNA) == "ATP-HILIC" | 
                                        colnames(matrix_trim_noNA) == "GTP-HILIC" | 
                                        colnames(matrix_trim_noNA) == "Uridine-HILIC" | 
                                        colnames(matrix_trim_noNA) == "NAD-HILIC" | 
                                        colnames(matrix_trim_noNA) == "AMP-HILIC" | 
                                        colnames(matrix_trim_noNA) == "FAD-HILIC" | 
                                        colnames(matrix_trim_noNA) == "Guanine-HILIC" | 
                                        colnames(matrix_trim_noNA) == "Adenine-HILIC" | 
                                        colnames(matrix_trim_noNA) == "Thymine-HILIC" | 
                                        colnames(matrix_trim_noNA) == "Adenosine-HILIC" |
                                        colnames(matrix_trim_noNA) == "Cytidine-HILIC" |
                                        colnames(matrix_trim_noNA) == "Guanosine-HILIC" |
                                        colnames(matrix_trim_noNA) == "Glycerol-3-Phosphate-HILIC" |
                                        colnames(matrix_trim_noNA) == "Ribose-5-Phosphate-HILIC" |
                                        colnames(matrix_trim_noNA) == "Acetyl CoA-HILIC" |
                                        colnames(matrix_trim_noNA) == "Glucose-6-Phosphate-HILIC" |
                                        colnames(matrix_trim_noNA) == "Fructose-6-Phosphate-HILIC" |
                                        colnames(matrix_trim_noNA) == "UDP-Glucosamine-HILIC" |
                                        colnames(matrix_trim_noNA) == "UDP-Glucose-HILIC" |
                                        colnames(matrix_trim_noNA) == "Citric Acid-HILIC" |
                                        colnames(matrix_trim_noNA) == "Ketoglutaric Acid-HILIC" |
                                        colnames(matrix_trim_noNA) == "Isocitric Acid-HILIC" |
                                        colnames(matrix_trim_noNA) == "PEP-HILIC" |
                                        colnames(matrix_trim_noNA) == "Malic Acid-HILIC" |
                                        colnames(matrix_trim_noNA) == "Aconitic Acid-HILIC" |
                                        colnames(matrix_trim_noNA) == "DHAP-HILIC" |
                                        colnames(matrix_trim_noNA) == "Phosphoglyceric Acid-HILIC" |
                                        colnames(matrix_trim_noNA) == "Fumaric Acid-HILIC"
                                        ) ]

#Perform the LDA here
group.lda<-lda(wide.matrix, meta.dat$Group)
#group.lda
#plot(group.lda)
group.lda$svd^2/sum(group.lda$svd^2)

#This gives us the predicted group membership and the probability of membership in each group 
#according to our LDA, which we can use to evaluate the efficacy of the discrimination
group.lda.pred<-predict(group.lda)
scores<-group.lda.pred$x

#To test if the LDA can predict the groups
result<-as.data.frame(cbind(meta.dat$Group,scores))
env.table<-table(meta.dat$Group,group.lda.pred$class)
overall.ccr <- sum(diag(env.table))/sum(env.table)
print(paste0(overall.ccr*100, " = overall classification rate (%)"))

#Now to look at the coefficients,
#We can define a canonical function on the basis of the structure coefficients 
#by noting the variables that have the largest loadings. 
coeff.vars <- lda.structure(scores, wide.matrix) 
coeff.vars.2 <- coeff.vars %>% as.data.frame() %>%
  mutate(compound = row.names(coeff.vars))

#To tell if its positive or negative loading that's associated with in or out group assignment, check here
plot(scores[,c("LD1", "LD2")],type='p', pch=16, col=colors)
text(scores[,c("LD1", "LD2")],labels=meta.dat$Organism, cex=.6, col=colors)


#ordination of LD1 and LD2
main_lda = 
  ggplot(as.data.frame(scores[,c("LD1", "LD2")])) + 
  aes(LD1, LD2, fill=as.factor(meta.dat$Group), color=as.factor(meta.dat$Group)) +
  scale_x_continuous(name="LD1", limits=c(-16, 16), breaks=seq(-16, 16, 4)) +
  scale_y_continuous(name="LD2", limits=c(-16, 16), breaks=seq(-16, 16, 4)) +
  geom_point() +
  scale_color_manual(values = pal) +
  geom_text(aes(label=meta.dat$Organism), cex=3) +
  theme(legend.title = element_blank())

#boxplots of lda groups

lda1_box = ggplot(aes(x = reorder(as.factor(V1), as.numeric(LD1), median), y = as.numeric(LD1), 
           fill=as.factor(V1), color=as.factor(V1)), data=result) +
  geom_boxplot() +
  scale_y_continuous(limits=c(-16, 16), breaks=seq(-16, 16, 4)) +
  scale_x_discrete(breaks = NULL) +
  scale_fill_manual(values = alpha(pal, 0.5))+
  scale_color_manual(values = pal)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank()) +
  rotate()


lda2_box = ggplot(aes(x = reorder(as.factor(V1), as.numeric(LD2), median), y = as.numeric(LD2), 
           fill=as.factor(V1), color=as.factor(V1)), data=result) +
  geom_boxplot() +
  scale_y_continuous(limits=c(-16, 16), breaks=seq(-16, 16, 4)) +
  scale_x_discrete(breaks = NULL) +
  scale_fill_manual(values = alpha(pal, 0.5)) +
  scale_color_manual(values = pal) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank())

#plot all together
pdf(file="LDA.pdf")
ggarrange(lda2_box, main_lda, NULL, lda1_box,
          ncol = 2, nrow = 2,  align = "hv", 
          widths = c(1, 2), heights = c(2, 1), legend="none")
dev.off()



#table of LD coeff values
LDs = coeff.vars.2
write.table(LDs, file="structure_correlations_1.csv", sep=",")







#remove metabolites that aren't important in LDs
coeff.vars.2_cutoff = coeff.vars
coeff.vars.2_cutoff = as.data.frame(sapply(coeff.vars.2_cutoff, as.numeric), row.names=rownames(coeff.vars.2_cutoff))
coeff.vars.2_cutoff$abs_max = apply(abs(coeff.vars.2_cutoff), 1, max)
coeff.vars.2_cutoff$compound = rownames(coeff.vars.2_cutoff)
observed_2 <- coeff.vars.2_cutoff[which(coeff.vars.2_cutoff$abs_max >0.5), ]
good.compounds_2 <- observed_2$compound
wide.matrix_2 <- wide.matrix[, good.compounds_2]


#Perform the LDA here
group.lda2<-lda(wide.matrix_2, meta.dat$Group)
#group.lda
#plot(group.lda)
group.lda2$svd^2/sum(group.lda2$svd^2)

#This gives us the predicted group membership and the probability of membership in each group according to our LDA, which we can use to evaluate the efficacy of the discrimination
group.lda2.pred<-predict(group.lda2)
scores2<-group.lda2.pred$x

#To test if the LDA can predict the groups
result2<-as.data.frame(cbind(meta.dat$Group,scores2))
env.table2<-table(meta.dat$Group,group.lda2.pred$class)
overall.ccr2 <- sum(diag(env.table2))/sum(env.table2)
print(paste0(overall.ccr2*100, " = overall classification rate (%)"))

#Now to look at the coefficients,
#We can define a canonical function on the basis of the structure coefficients by noting the variables that have the largest loadings. 
coeff.vars_red <- lda.structure(scores2, wide.matrix_2) 
coeff.vars.2_red <- coeff.vars_red %>% as.data.frame() %>%
  mutate(compound = row.names(coeff.vars_red))


#To tell if its positive or negative loading that's associated with in or out group assignment, check here
plot(scores2[,c("LD1", "LD2")],type='p', pch=16, col=colors)
text(scores2[,c("LD1", "LD2")],labels=meta.dat$Organism, cex=.6, col=colors)

#ordination of LD1 and LD2
main_lda_red = 
  ggplot(as.data.frame(scores2[,c("LD1", "LD2")])) + 
  aes(LD1, LD2, fill=as.factor(meta.dat$Group), color=as.factor(meta.dat$Group)) +
  scale_x_continuous(name="LD1", limits=c(-12,12), breaks=seq(-12, 12, 4)) +
  scale_y_continuous(name="LD2", limits=c(-12,12), breaks=seq(-12, 12, 4)) +
  geom_point() +
  scale_color_manual(values = pal) +
  geom_text(aes(label=meta.dat$Organism), cex=3) +
  theme(legend.title = element_blank())



#boxplots of lda groups
lda1_box_red = ggplot(aes(x = reorder(as.factor(V1), as.numeric(LD1), median), y = as.numeric(LD1), 
           fill=as.factor(V1), color=as.factor(V1)), data=result2) +
  geom_boxplot() +
  scale_y_continuous(limits=c(-12,12), breaks=seq(-12, 12, 4)) +
  scale_x_discrete(breaks = NULL) +
  scale_fill_manual(values = alpha(pal, 0.5))+
  scale_color_manual(values = pal)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank()) +
  rotate()  


lda2_box_red = ggplot(aes(x = reorder(as.factor(V1), as.numeric(LD2), median), y = as.numeric(LD2), 
           fill=as.factor(V1), color=as.factor(V1)), data=result2) +
  geom_boxplot() +
  scale_y_continuous(limits=c(-12,12), breaks=seq(-12, 12, 4)) +
  scale_x_discrete(breaks = NULL) +
  scale_fill_manual(values = alpha(pal, 0.5))+
  scale_color_manual(values = pal)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank())


#plot all together
pdf(file="LDA_red.pdf")
ggarrange(lda2_box_red, main_lda_red, NULL, lda1_box_red, 
          ncol = 2, nrow = 2,  align = "hv", 
          widths = c(1, 2), heights = c(2, 1), legend="none")
dev.off()


#table of LD coeff values
LDs_red = coeff.vars.2_red
write.table(LDs_red, file="structure_correlations_red.csv", sep=",")



LDs_red_plot = data.frame(sapply(LDs_red[,1:2], as.numeric), row.names = LDs_red[,7])

color.palette  <- colorRampPalette(c("red", "white", "blue"))

heatmap.2(as.matrix(LDs_red_plot), trace = "none", dendrogram = "row", key = T,
          cellnote=as.matrix(LDs_red_plot), notecol="black", col=color.palette,
          Colv=FALSE, margins = c(3,8), cexRow=0.6, cexCol = 1)

pdf(file='coefficientsofLDs_red_heatmap.pdf', width=6, height=6)
heatmap.2(as.matrix(LDs_red_plot), trace = "none", dendrogram = "row", key = T,
          cellnote=as.matrix(LDs_red_plot), notecol="black", col=color.palette,
          Colv=FALSE, margins = c(3,8), cexRow=0.6, cexCol = 1)
dev.off()




LDs.dataframe= data.frame(sapply(LDs[,1:2], as.numeric), row.names = LDs[,7])
set.compounds <- LDs_red$compound
set.compounds.LDs <- LDs.dataframe[set.compounds,]

color.palette  <- colorRampPalette(c("red", "white", "blue"))

pdf(file='coefficientsofLDs_heatmap.pdf', width=6, height=6)
heatmap.2(as.matrix(set.compounds.LDs), trace = "none", dendrogram = "row", key = T,
          cellnote=as.matrix(set.compounds.LDs), notecol="black", col=color.palette,
          Colv=FALSE, margins = c(3,8), cexRow=0.6, cexCol = 1)
dev.off()

