
library(matrixStats)
library(FactoMineR)
library(factoextra)
library(ggpubr)

#Read in PA dat, mudge to get into a matrix, toss any compounds that weren't seen ever, make NAs 0s
PA.dat <- cbind(rownames(dat6), dat6)
colnames(PA.dat)[1] = "Organism"
meta.dat <- PA.dat %>% dplyr::select(Organism, Domain, Group) 
matrix_trim_noNA = matrix_trim
matrix_trim_noNA[is.na(matrix_trim_noNA)] <- 0

colnames(matrix_trim_noNA) = gsub("-HILIC", "", colnames(matrix_trim_noNA))
colnames(matrix_trim_noNA) = gsub("-RP", "", colnames(matrix_trim_noNA))

#remove universal metabolites
wide.matrix = matrix_trim_noNA[, -which(colnames(matrix_trim_noNA) == "cAMP" | 
                                          colnames(matrix_trim_noNA) == "NADH" | 
                                          colnames(matrix_trim_noNA) == "Thymidine" | 
                                          colnames(matrix_trim_noNA) == "Uracil" | 
                                          colnames(matrix_trim_noNA) == "cGMP" | 
                                          colnames(matrix_trim_noNA) == "GMP" | 
                                          colnames(matrix_trim_noNA) == "ADP" | 
                                          colnames(matrix_trim_noNA) == "NADP" | 
                                          colnames(matrix_trim_noNA) == "ATP" | 
                                          colnames(matrix_trim_noNA) == "GTP" | 
                                          colnames(matrix_trim_noNA) == "Uridine" | 
                                          colnames(matrix_trim_noNA) == "NAD" | 
                                          colnames(matrix_trim_noNA) == "AMP" | 
                                          colnames(matrix_trim_noNA) == "FAD" | 
                                          colnames(matrix_trim_noNA) == "Guanine" | 
                                          colnames(matrix_trim_noNA) == "Adenine" | 
                                          colnames(matrix_trim_noNA) == "Thymine" | 
                                          colnames(matrix_trim_noNA) == "Adenosine" |
                                          colnames(matrix_trim_noNA) == "Cytidine" |
                                          colnames(matrix_trim_noNA) == "Guanosine" |
                                          colnames(matrix_trim_noNA) == "Glycerol-3-Phosphate" |
                                          colnames(matrix_trim_noNA) == "Ribose-5-Phosphate" |
                                          colnames(matrix_trim_noNA) == "Acetyl CoA" |
                                          colnames(matrix_trim_noNA) == "Glucose-6-Phosphate" |
                                          colnames(matrix_trim_noNA) == "Fructose-6-Phosphate" |
                                          colnames(matrix_trim_noNA) == "UDP-Glucosamine" |
                                          colnames(matrix_trim_noNA) == "UDP-Glucose" |
                                          colnames(matrix_trim_noNA) == "Citric Acid" |
                                          colnames(matrix_trim_noNA) == "Ketoglutaric Acid" |
                                          colnames(matrix_trim_noNA) == "Isocitric Acid" |
                                          colnames(matrix_trim_noNA) == "PEP" |
                                          colnames(matrix_trim_noNA) == "Malic Acid" |
                                          colnames(matrix_trim_noNA) == "Aconitic Acid" |
                                          colnames(matrix_trim_noNA) == "DHAP" |
                                          colnames(matrix_trim_noNA) == "Phosphoglyceric Acid" |
                                          colnames(matrix_trim_noNA) == "Fumaric Acid"
                                          ) ]



test.dat <- merge(meta.dat, wide.matrix, by='row.names')
colnames(test.dat)[1] = "Organism"
rownames(test.dat) = test.dat[,1]
test.dat = test.dat[,-c(1:3)]
test.dat[] <- lapply( test.dat, factor)


res.mca <- MCA(test.dat, quali.sup=1, graph = FALSE)

res.mca$quali.sup

eig.val <- get_eigenvalue(res.mca)
fviz_screeplot(res.mca, addlabels = TRUE, ylim = c(0, 45))

fviz_mca_var(res.mca, choice = "mca.cor",
             repel = TRUE)

fviz_mca_var(res.mca, repel = TRUE,
             ggtheme= theme_minimal())


all_names <- list(name = c(
  "Acetyl-L-Carnitine_1", 
  "DHPS_1", 
  "Propionyl-L-carnitine_1", 
  "Carnitine_1", 
  "DMS-Acetate_1", 
  "trans-Hydroxylproline_1", 
  "Indole-3-Acetic Acid_1", 
  "Cystine_1", 
  "DMSP_1", 
  "Choline_1", 
  "Homarine_1", 
  "Coenzyme Q10_1", 
  "Glycerophosphocholine_1", 
  "Glucosylglycerol_1",
  "Proline Betaine_1",
  "2iP_1", 
  "Aminobutyric Acid_1", 
  "Chitobiose_1", 
  "Dimethyl Glycine_1", 
  "Methionine Sulfoxide_1", 
  "Gonyol_1", 
  "Taurine_1", 
  "Sarcosine_1", 
  "Carotene_1", 
  "Picolinic Acid_1", 
  "Sulfolactic Acid_1", 
  "Desthiobiotin_1", 
  "Oxalic Acid_1", 
  "Cystathionine_1", 
  "Sucrose_1", 
  "Trigonelline_1", 
  "Xanthine_1", 
  "Isethionic Acid_1", 
  "Glutathione Disulfide_1", 
  "Trehalose_1", 
  "Gluconic Acid_1", 
  "Betaine_1", 
  "Kynurenine_1", 
  "Tryptophol_1", 
  "Argininosuccinic Acid_1", 
  "Orotic Acid_1", 
  "Indole-3-Carboxylic Acid_1",
  "Glutathione_1", 
  "Cyclic diGMP_1", 
  "Ectoine_1", 
  "N-Acetyl-Lysine_1", 
  "Cysteic Acid_1", 
  "Homoserine_1", 
  "Ornithine_1", 
  "Vitamin C_1", 
  "Indole-3-Acetamide_1", 
  "Sulfurol_1", 
  "Citrulline_1"))




# Top active variables with the highest cos2

fviz_mca_biplot(res.mca, axes=1:2, col.var = "contrib",
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), label='var',
                ggtheme = theme_minimal())

fviz_mca_biplot(res.mca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             select.var= list(contrib = 54),
             select.ind = list(contrib = 45),
             color.ind=colors,
             ggtheme = theme_minimal())

fviz_mca_ind(res.mca,
             labels=FALSE,
             habillage = "Group", # color by groups 
             palette = pal,
             ggtheme = theme_minimal()) 





cos2_round1 = as.data.frame(res.mca[["var"]][["cos2"]])
cos2_round1$abs_max = apply(abs(cos2_round1), 1, max)
cos2_round1 = cos2_round1[-grep('_0', row.names(cos2_round1)),, drop = FALSE]
rownames(cos2_round1) = gsub("_1", "", rownames(cos2_round1))
cos2_round1 = as.data.frame(cos2_round1)
cos2_round1_keepers <- cos2_round1[which(cos2_round1$abs_max>.2), ]

test.dat2 <- test.dat[, rownames(cos2_round1_keepers)]
test.dat2 <- merge(meta.dat, test.dat2, by='row.names')
colnames(test.dat2)[1] = "Organism"
rownames(test.dat2) = test.dat2[,1]
test.dat2 = test.dat2[,-c(1:3)]
test.dat2[] <- lapply( test.dat2, factor)

res.mca2 <- MCA(test.dat2, quali.sup=1, graph = FALSE)

eig.val2 <- get_eigenvalue(res.mca2)
fviz_screeplot(res.mca2, addlabels = TRUE, ylim = c(0, 45))

fviz_mca_ind(res.mca2,
             labels=FALSE,
             habillage = "Group", # color by groups 
             palette = pal,
             ggtheme = theme_minimal()) 



cos2_round1_names = as.data.frame(res.mca[["var"]][["cos2"]])
cos2_round1_names$abs_max = apply(abs(cos2_round1_names[,1:2]), 1, max)
cos2_round1_names = cos2_round1_names[-grep('_0', row.names(cos2_round1_names)),, drop = FALSE]
cos2_round1_names = as.data.frame(cos2_round1_names)
cos2_round1_names <- cos2_round1_names[which(cos2_round1_names$abs_max>.2), ]

name <- list(name = rownames(cos2_round1_names))



cos2_round2 = as.data.frame(res.mca2[["var"]][["cos2"]])
cos2_round2$abs_max = apply(abs(cos2_round2), 1, max)
cos2_round2 = cos2_round2[-grep('_0', row.names(cos2_round2)),, drop = FALSE]
rownames(cos2_round2) = gsub("_1", "", rownames(cos2_round2))
cos2_round2 = as.data.frame(cos2_round2)

cos2_round2_names = as.data.frame(res.mca2[["var"]][["cos2"]])
cos2_round2_names$abs_max = apply(abs(cos2_round2_names[,1:2]), 1, max)
cos2_round2_names = cos2_round2_names[-grep('_0', row.names(cos2_round2_names)),, drop = FALSE]
cos2_round2_names = as.data.frame(cos2_round2_names)
cos2_round2_names <- cos2_round2_names[which(cos2_round2_names$abs_max>.2), ]

name2 <- list(name = rownames(cos2_round2_names))


colors1 = c('navy', 'navy', 'navy', 'navy', 'navy', 'navy', 'navy', 'navy', 'navy', 'navy', 'navy', 
           'gold', 'navy', 'navy',
           'deepskyblue3', 'deepskyblue3', 'deepskyblue3', 'deepskyblue3', 'deepskyblue3', 'purple', 
           'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 
           'darkorange', 'darkorange', 'darkorange', 'darkorange', 'navy', 'navy',
           'green2',  'green2', 'deepskyblue3',  'deepskyblue3',
           'purple', 'purple', 'purple', 'purple', 'purple', 'purple')

fviz_mca_ind(res.mca2,
             labels=FALSE,
             habillage = "Group", # color by groups 
             palette = pal,
             ggtheme = theme_minimal()) 

mca1_ind_df = data.frame(res.mca$ind$coord)
mca1_ind_df <- merge(meta.dat, mca1_ind_df, by='row.names')
rownames(mca1_ind_df) = mca1_ind_df[,1]

mca2_ind_df = data.frame(res.mca2$ind$coord)
mca2_ind_df <- merge(meta.dat, mca2_ind_df, by='row.names')
rownames(mca2_ind_df) = mca2_ind_df[,1]


mca1_box1 = ggplot(aes(x = reorder(as.factor(Group), as.numeric(Dim.1), median), y = as.numeric(Dim.1), 
                      fill=as.factor(Group), color=as.factor(Group)), data=mca1_ind_df[,-c(1:3)]) +
  geom_boxplot() +
  scale_y_continuous(limits=c(-1.1, 1.4), breaks=seq(-1.5, 1.5, 0.5)) +
  scale_x_discrete(breaks = NULL) +
  scale_fill_manual(values = alpha(pal, 0.5))+
  scale_color_manual(values = pal)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank()) + rotate()

mca1_box2 = ggplot(aes(x = reorder(as.factor(Group), as.numeric(Dim.2), median), y = as.numeric(Dim.2), 
                       fill=as.factor(Group), color=as.factor(Group)), data=mca1_ind_df[,-c(1:3)]) +
  geom_boxplot() +
  scale_y_continuous(limits=c(-0.9, 1.4), breaks=seq(-1, 1.5, 0.5)) +
  scale_x_discrete(breaks = NULL) +
  scale_fill_manual(values = alpha(pal, 0.5))+
  scale_color_manual(values = pal)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank())

mca2_box1 = ggplot(aes(x = reorder(as.factor(Group), as.numeric(Dim.1), median), y = as.numeric(Dim.1), 
                      fill=as.factor(Group), color=as.factor(Group)), data=mca2_ind_df[,-c(1:3)]) +
  geom_boxplot() +
  scale_y_continuous(limits=c(-1.1, 1.4), breaks=seq(-1.5, 1.5, 0.5)) +
  scale_x_discrete(breaks = NULL) +
  scale_fill_manual(values = alpha(pal, 0.5))+
  scale_color_manual(values = pal)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank()) + rotate()

mca2_box2 = ggplot(aes(x = reorder(as.factor(Group), as.numeric(Dim.2), median), y = as.numeric(Dim.2), 
                       fill=as.factor(Group), color=as.factor(Group)), data=mca2_ind_df[,-c(1:3)]) +
  geom_boxplot() +
  scale_y_continuous(limits=c(-1.4, 0.9), breaks=seq(-1.5, 1, 0.5)) +
  scale_x_discrete(breaks = NULL) +
  scale_fill_manual(values = alpha(pal, 0.5))+
  scale_color_manual(values = pal)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank())


plot1 = fviz_mca_biplot(res.mca, axes=1:2, select.var= name2, col.var = "cos2",
                        gradient.cols = c("gray90", "black"), label='var', labelsize=3,
                        ggtheme = theme_classic()) + 
                        scale_x_continuous(limits=c(-1.1, 1.4), breaks=seq(-1.5, 1.5, 0.5)) +
                        scale_y_continuous(limits=c(-0.9, 1.4), breaks=seq(-1, 1.5, 0.5))


plot2 = fviz_mca_biplot(res.mca2, axes=1:2, select.var= name2, col.var = "cos2",
                        gradient.cols = c("gray90", "black"), label='var', labelsize=3,
                        ggtheme = theme_classic()) + 
                        scale_x_continuous(limits=c(-1.1, 1.4), breaks=seq(-1.5, 1.5, 0.5)) +
                        scale_y_continuous(limits=c(-1.4, 0.9), breaks=seq(-1.5, 1, 0.5))


mca1_combo = fviz_add(plot1, mca1_ind_df[,5:9], addlabel=FALSE, col=colors1)
mca2_combo = fviz_add(plot2, mca2_ind_df[,5:9], addlabel=FALSE, col=colors1)

pdf(file="MCA1.pdf")
ggarrange(mca1_box2, mca1_combo, NULL, mca1_box1,
          ncol = 2, nrow = 2,  align = "hv", 
          widths = c(1, 2), heights = c(2, 1), legend="none")
dev.off()


pdf(file="MCA2.pdf")
ggarrange(mca2_box2, mca2_combo, NULL, mca2_box1,
          ncol = 2, nrow = 2,  align = "hv", 
          widths = c(1, 2), heights = c(2, 1), legend="none")
dev.off()



pdf(file="mca1-scale.pdf")
mca1_combo
dev.off()

pdf(file="mca2-scale.pdf")
mca2_combo
dev.off()
