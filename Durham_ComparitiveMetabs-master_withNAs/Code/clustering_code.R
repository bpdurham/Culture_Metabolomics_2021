
#get the packages uploaded
install.packages("recluster")
install.packages("phytools")
install.packages("vegan")

library("vegan") 
library("cluster")
library("pvclust")
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(here)
library(RColorBrewer)
library(gplots)
library(factoextra)


#upload and transform data into presence/absence matrix
#Blk=present at background (media blank) level, 0=not present, 1=present in all replicates, fractions=fraction of replicates that contain metabolite
dat = read.table(file='~/Desktop/paper_files/All_PA_features-2021-12-3.csv', header = T, row.names = NULL, sep=',', stringsAsFactors = FALSE)
unique(dat$notes)
#remove rows with metabolites we don't want to use (e.g., provided in media or were not monitored in all samples) and remove unncessary extra columns
dat0 = dat[which(dat$notes == 'QE-TQS mix' | dat$notes == 'QE' | dat$notes == 'TQS'),]
unique(dat0$notes)
dat1 = dat0[ , -which(names(dat0) %in% c("notes"))]

#manually remove compounds that had high background (those with more than 20 Blk values)
dat2 = dat1[-which(dat1$Compound.Name == "Sulfoacetic Acid-HILIC" | dat1$Compound.Name == "Succinic Acid-HILIC" | dat1$Compound.Name == "4-Hydroxybenzaldehyde-RP"), ]

#Change Blk to NA since we could not detect the metabolite above background levels
dat3 = data.frame(mutate_if(dat2, 
                is.character, 
                str_replace_all, pattern = "Blk", replacement = "NA"), row.names=1)
dat4 = data.frame(sapply(dat3, as.numeric), row.names = rownames(dat3))
#Need to see the metabolite in more than half of the replicates to count as present
dat5 = data.frame(sapply(dat4, function(x) 
  ifelse(x <= 0.5, 0, ifelse(x > 0.5, 1, x))), row.names = rownames(dat4))

#add group information
Domain = c(
  "Eukaryote",
  "Eukaryote",
  "Archaea",
  "Bacteria",
  "Bacteria",
  "Bacteria",
  "Bacteria",
  "Bacteria",
  "Bacteria",
  "Bacteria",
  "Bacteria",
  "Bacteria",
  "Bacteria",
  "Bacteria",
  "Bacteria",
  "Bacteria",
  "Bacteria",
  "Bacteria",
  "Bacteria",
  "Bacteria",
  "Bacteria",
  "Bacteria",
  "Bacteria",
  "Eukaryote",
  "Eukaryote",
  "Eukaryote",
  "Eukaryote",
  "Eukaryote",
  "Eukaryote",
  "Eukaryote",
  "Eukaryote",
  "Eukaryote",
  "Eukaryote",
  "Eukaryote",
  "Eukaryote",
  "Eukaryote",
  "Eukaryote",
  "Eukaryote",
  "Eukaryote",
  "Eukaryote",
  "Bacteria",
  "Bacteria",
  "Eukaryote",
  "Eukaryote",
  "Eukaryote")

Group = c(
  "Dinoflagellate",
  "Prasinophyte",
  "Archaea",
  "Heterotroph",
  "Heterotroph",
  "Heterotroph",
  "Heterotroph",
  "Heterotroph",
  "Heterotroph",
  "Heterotroph",
  "Heterotroph",
  "Heterotroph",
  "Heterotroph",
  "Heterotroph",
  "Heterotroph",
  "Heterotroph",
  "Cyanobacteria",
  "Cyanobacteria",
  "Cyanobacteria",
  "Cyanobacteria",
  "Cyanobacteria",
  "Cyanobacteria",
  "Cyanobacteria",
  "Diatom",
  "Haptophyte",
  "Diatom",
  "Diatom",
  "Haptophyte",
  "Diatom",
  "Diatom",
  "Diatom",
  "Dinoflagellate",
  "Dinoflagellate",
  "Dinoflagellate",
  "Haptophyte",
  "Haptophyte",
  "Haptophyte",
  "Haptophyte",
  "Haptophyte",
  "Prasinophyte",
  "Heterotroph",
  "Heterotroph",
  "Diatom",
  "Diatom",
  "Diatom")

bind = rbind(Group, dat5)
rownames(bind)[1] = "Group"
bind2 = rbind(Domain, bind)
rownames(bind2)[1] = "Domain"
dat6 = t(bind2)
dat6 = as.data.frame(dat6)
meta.matrix <- dat6 %>% dplyr::select(Domain, Group)

#create matrix for analyses
dat7 = dat6[ , -which(names(dat6) %in% c("Domain", "Group"))]
dat8 = data.frame(sapply(dat7, as.numeric), row.names = rownames(dat7))
colnames(dat8) = colnames(dat7)
matrix = as.matrix(dat8)

#remove uniques where a compound was only seen in one or two organisms and remove compounds present in all
observed <- t(matrix) %>% as.data.frame() %>%
  mutate(Compound = colnames(matrix)) %>%
  dplyr::select(Compound)%>%
  mutate(total = rowSums(t(matrix), na.rm = TRUE)) %>%
  filter(total < length(row.names(matrix))) %>%
  filter(total > 2)

good.compounds <- observed$Compound
matrix_trim <- matrix[, good.compounds]



#metabolites in phytoplankton only
subset = matrix[which(rownames(matrix) == 'Gamma.SA7' | rownames(matrix) == 'Gamma.SA55' | 
                         rownames(matrix) == 'Alpha.SA30' | rownames(matrix) == 'Alpha.SA42' | 
                         rownames(matrix) == 'Alpha.SA33' | rownames(matrix) == 'Beta.SA59' | 
                         rownames(matrix) == 'CFB.SA60' | rownames(matrix) == 'Alpha.SA16' |
                         rownames(matrix) == 'Alpha.SA44' | rownames(matrix) == 'Alpha.SA36' |
                         rownames(matrix) == 'Alpha.SA48' | rownames(matrix) == 'Alpha.SA53' |
                         rownames(matrix) == 'Alpha.SA11' | rownames(matrix) == 'Alpha.DSS.3' |
                         rownames(matrix) == 'Alpha.Och114' | rownames(matrix) == 'Archaea.SCM1'),]

subset_trim <- subset[,which(colSums(subset, na.rm=TRUE) == 0)]
colnames(subset_trim)

#metabolites in eukaryotic phytoplankton only
subset1 = matrix[which(rownames(matrix) == 'Gamma.SA7' | rownames(matrix) == 'Gamma.SA55' | 
                         rownames(matrix) == 'Alpha.SA30' | rownames(matrix) == 'Alpha.SA42' | 
                         rownames(matrix) == 'Alpha.SA33' | rownames(matrix) == 'Beta.SA59' | 
                         rownames(matrix) == 'CFB.SA60' | rownames(matrix) == 'Alpha.SA16' |
                         rownames(matrix) == 'Alpha.SA44' | rownames(matrix) == 'Alpha.SA36' |
                         rownames(matrix) == 'Alpha.SA48' | rownames(matrix) == 'Alpha.SA53' |
                         rownames(matrix) == 'Alpha.SA11' | rownames(matrix) == 'Alpha.DSS.3' |
                         rownames(matrix) == 'Alpha.Och114' | rownames(matrix) == 'Archaea.SCM1' |
                         rownames(matrix) == 'Syn.7803' | rownames(matrix) == 'Croco.8501' | 
                         rownames(matrix) == 'Syn.8102' | rownames(matrix) == 'Pro.MED4' | 
                         rownames(matrix) == 'Pro.1314' | rownames(matrix) == 'Pro.AS9601' |
                         rownames(matrix) == 'Pro.Natl2A'),]
                  

subset1_trim <- subset1[,which(colSums(subset1, na.rm=TRUE) == 0)]
colnames(subset1_trim)

#metabolites in heterotrophic bacteria only

het_sum = colSums(matrix[which(rownames(matrix) == 'Gamma.SA7' | rownames(matrix) == 'Gamma.SA55' | 
                                  rownames(matrix) == 'Alpha.SA30' | rownames(matrix) == 'Alpha.SA42' | 
                                  rownames(matrix) == 'Alpha.SA33' | rownames(matrix) == 'Beta.SA59' | 
                                  rownames(matrix) == 'CFB.SA60' | rownames(matrix) == 'Alpha.SA16' |
                                  rownames(matrix) == 'Alpha.SA44' | rownames(matrix) == 'Alpha.SA36' |
                                  rownames(matrix) == 'Alpha.SA48' | rownames(matrix) == 'Alpha.SA53' |
                                  rownames(matrix) == 'Alpha.SA11' | rownames(matrix) == 'Alpha.DSS.3' |
                                  rownames(matrix) == 'Alpha.Och114'),], na.rm=TRUE)

phyto_archaea_sum = colSums(matrix[which(rownames(matrix) == 'Dino.1314' | rownames(matrix) == 'Green.1545' | 
                                   rownames(matrix) == 'Diatom.To' |
                                   rownames(matrix) == 'Hapto.2090' | rownames(matrix) == 'Diatom.Cm' |
                                   rownames(matrix) == 'Diatom.Tp' | rownames(matrix) == 'Hapto.371' |
                                   rownames(matrix) == 'Diatom.Np' | rownames(matrix) == 'Diatom.Pt' |
                                   rownames(matrix) == 'Diatom.Pc55x' | rownames(matrix) == 'Dino.449' |
                                   rownames(matrix) == 'Dino.1771' | rownames(matrix) == 'Dino.2021' |
                                   rownames(matrix) == 'Green.3430' | rownames(matrix) == 'Syn.7803' | 
                                   rownames(matrix) == 'Croco.8501' | rownames(matrix) == 'Archaea.SCM1' |
                                   rownames(matrix) == 'Syn.8102' | rownames(matrix) == 'Pro.MED4' | 
                                   rownames(matrix) == 'Pro.1314' | rownames(matrix) == 'Pro.AS9601' |
                                   rownames(matrix) == 'Pro.Natl2A'),], na.rm=TRUE)

sums1 = rbind(matrix, het_sum)
sums4 = rbind(sums1, phyto_archaea_sum)
colnames(sums4[,which(sums4[46,] >= 1 & sums4[47,] == 0)])

#metabolites in archaea and bacteria only

prok_archaea_sum = colSums(matrix[which(rownames(matrix) == 'Gamma.SA7' | rownames(matrix) == 'Gamma.SA55' | 
                                  rownames(matrix) == 'Alpha.SA30' | rownames(matrix) == 'Alpha.SA42' | 
                                  rownames(matrix) == 'Alpha.SA33' | rownames(matrix) == 'Beta.SA59' | 
                                  rownames(matrix) == 'CFB.SA60' | rownames(matrix) == 'Alpha.SA16' |
                                  rownames(matrix) == 'Alpha.SA44' | rownames(matrix) == 'Alpha.SA36' |
                                  rownames(matrix) == 'Alpha.SA48' | rownames(matrix) == 'Alpha.SA53' |
                                  rownames(matrix) == 'Alpha.SA11' | rownames(matrix) == 'Alpha.DSS.3' |
                                  rownames(matrix) == 'Alpha.Och114' | rownames(matrix) == 'Archaea.SCM1' |
                                  rownames(matrix) == 'Syn.7803' | rownames(matrix) == 'Croco.8501' | 
                                  rownames(matrix) == 'Syn.8102' | rownames(matrix) == 'Pro.MED4' | 
                                  rownames(matrix) == 'Pro.1314' | rownames(matrix) == 'Pro.AS9601' |
                                  rownames(matrix) == 'Pro.Natl2A'),], na.rm=TRUE)

euk_sum = colSums(matrix[which(rownames(matrix) == 'Dino.1314' | rownames(matrix) == 'Green.1545' | 
                                 rownames(matrix) == 'Diatom.To' |
                                 rownames(matrix) == 'Hapto.2090' | rownames(matrix) == 'Diatom.Cm' |
                                 rownames(matrix) == 'Diatom.Tp' | rownames(matrix) == 'Hapto.371' |
                                 rownames(matrix) == 'Diatom.Np' | rownames(matrix) == 'Diatom.Pt' |
                                 rownames(matrix) == 'Diatom.Pc55x' | rownames(matrix) == 'Dino.449' |
                                 rownames(matrix) == 'Dino.1771' | rownames(matrix) == 'Dino.2021' |
                                 rownames(matrix) == 'Green.3430'),], na.rm=TRUE)

sums = rbind(matrix, prok_archaea_sum)
sums2 = rbind(sums, euk_sum)
colnames(sums2[,which(sums2[46,] >= 1 & sums2[47,] == 0)])

#make summary table of these summaries
sums2 = rbind(sums, euk_sum)
sums3 = rbind(sums2, phyto_archaea_sum)
sums4 = rbind(sums3, het_sum)
write.table(sums4, file='summary_groups.csv', sep=',')

#metabolites different in chryso strains
subset_chry = matrix[which(rownames(matrix) == 'Chryso.Ca' | rownames(matrix) == 'Chryso.Cr' |
                             rownames(matrix) == 'Chryso.Cs' | rownames(matrix) == 'Chryso.116' | 
                             rownames(matrix) == 'Chryso.P3' | rownames(matrix) == 'Chryso.P5'),]
marine_sum = colSums(subset_chry[which(rownames(subset_chry) == 'Chryso.Ca' | rownames(subset_chry) == 'Chryso.Cr' | rownames(subset_chry) == 'Chryso.Cs'),], na.rm=TRUE)
fresh_sum = colSums(subset_chry[which(rownames(subset_chry) == 'Chryso.116' | rownames(subset_chry) == 'Chryso.P3' | rownames(subset_chry) == 'Chryso.P5'),], na.rm=TRUE)
comb_chry1 = rbind(subset_chry, marine_sum)
comb_chry2 = rbind(comb_chry1, fresh_sum)
colnames(comb_chry2[,which(comb_chry2[6,] == 0 & comb_chry2[7,] >= 1)])
colnames(comb_chry2[,which(comb_chry2[6,] >= 1 & comb_chry2[7,] == 0)])

write.table(comb_chry2, file='chryso.txt', sep="\t")



#convert NAs to 0 for comparison
matrix_noNA = matrix
matrix_noNA[is.na(matrix_noNA)] <- 0


new_matrix <- matrix[ order(row.names(matrix)), ]
new_matrix_noNA = new_matrix
new_matrix_noNA[is.na(new_matrix_noNA)] <- 0
jacd_ordered <- vegdist(new_matrix_noNA, method="jaccard", na.rm=TRUE)
fviz_dist(jacd_ordered, order=FALSE, lab_size=6)

pdf(file='dist_matrix_noNA.pdf', height=8, width=8)
fviz_dist(jacd_ordered, order=FALSE, lab_size=6)
dev.off()



#test out some different distance methods, jaccard seems most appropriate
#simpson_dist <- recluster.dist(matrix_trim, dist="simpson")
#jaccard_dist <- recluster.dist(matrix_noNA, dist="jaccard")
#jaccard_dist2 <- recluster.dist(t(matrix_noNA), dist="jaccard")
#sorensen_dist <- recluster.dist(matrix_trim, dist="sorensen")
#euclidian_dist <- vegdist(matrix_trim, method="euclidean")
#brayd <- vegdist(matrix_trim, method = "bray")
jacd <- vegdist(matrix, method="jaccard", na.rm=TRUE)
jacd2 <- vegdist(t(matrix), method="jaccard", na.rm=TRUE)

#test out some different cluster methods, most generally get the same clusters
#clust = hclust(euclidian_dist, method='ward.D')
clust = hclust(jacd, method='ward.D2')
plot(clust)

clust2 = hclust(jacd2, method='ward.D2')
plot(clust2)

#look at p-values, bootstraps for your clusters
clus.stab <- pvclust(t(matrix),
                     method.hclust = "ward.D2",
                     method.dist = "binary",
                     use.cor="na.or.complete",
                     nboot = 100)
plot(clus.stab,
     hang=-1)
pvrect(clus.stab,alpha=0.95)


pdf(file='dendogram.pdf', height=6, width=8)
plot(clus.stab,
     hang=-1)
pvrect(clus.stab,alpha=0.95)
dev.off()

#make a heatmap
colors = c('darkorange', 'green2', 'gold', 'navy', 'navy', 'navy', 'navy', 'navy', 'navy', 'navy', 
           'navy', 'navy', 'navy', 'navy', 'navy', 'navy', 'purple', 'purple', 'purple', 'purple', 'purple', 
           'purple', 'purple', 'red', 'deepskyblue3', 'red', 'red', 'deepskyblue3', 'red', 'red', 'red', 'darkorange', 'darkorange', 'darkorange', 
           'deepskyblue3', 'deepskyblue3', 'deepskyblue3', 'deepskyblue3', 'deepskyblue3', 'green2', 'navy', 'navy', 'red', 'red', 'red')

heatmap.2(matrix, trace = "none", dendrogram = "both", scale='none', Rowv=as.dendrogram(clust), Colv=as.dendrogram(clust2), key = F,
          colRow = colors, col=c('white', 'cadetblue4'), na.color='gray',
          margins = c(12,4), cexRow=0.4, cexCol = 0.3)

pdf(file='heatmap_PA_withNA.pdf', height=6, width=7)
heatmap.2(matrix, trace = "none", dendrogram = "both", Rowv=as.dendrogram(clust), Colv=as.dendrogram(clust2), key = F,
          colRow = colors, col=c('white','cadetblue4'), na.color='gray',
          margins = c(12,4), cexRow=0.4, cexCol = 0.25)
dev.off()

pdf(file='heatmap_PA_noNA.pdf', height=6, width=7)
heatmap.2(matrix_noNA, trace = "none", dendrogram = "both", Rowv=as.dendrogram(clust), Colv=as.dendrogram(clust2), key = F,
          colRow = colors, col=c('white','cadetblue4'), na.color='gray',
          margins = c(12,4), cexRow=0.4, cexCol = 0.25)
dev.off()

