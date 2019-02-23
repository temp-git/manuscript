library(edgeR)
library(limma)
library(DESeq2)
library(sva)
library(dplyr)
library(data.table)
library(ggplot2)
library(gridExtra)
library(ggsci) 
library(ggpubr)

###################################################################
# All differentiation stages & multiple experiments (Fig 2A)
###################################################################

#count matrix of all genes
##Perez-Alcantara (2018)
marta <- read.delim("~/OneDrive/oxford/summer_internship/counts/marta.gene.counts.tsv", header = TRUE, row.names = 1, check.names = FALSE)
marta <- marta[,-1] # take out gene name
marta.donors = c("Ad2.1", "Ad3.1", "Neo1.1")  # data comes from three donors
marta.stage = c("iPSC", "DE", "PGT", "PFG", "PE", "EP", "ENstage6", "ENstage7") # 8 stages
order_by_stages = function(counts, stage. = marta.stage) {
  nc = counts[, grepl("iPSC" , names(counts))]  # takes columns whose name matches x
  
  # do the same for the other stages
  
  for (s in stage.[-1])  {
    i = counts[, grepl(s , names(counts))]
    nc = cbind(nc, i)
  }
  
  return(nc)
}
marta = order_by_stages(marta) # order
marta.stage2 = c("iPSC", "DE", "PGT", "PFG", "PE", "EP", "EN", "BLC") # rename
colnames(marta) = paste(rep(marta.stage2, each = 3), rep(marta.donors, 8), sep = "_")  # rename columns

##Xie (2013)
xiecounts <- read.delim("~/OneDrive/oxford/summer_internship/counts/stembancc.xie.gene.counts.tsv", header = TRUE, check.names = FALSE, row.names = 1)
xiecounts <- subset(xiecounts, xiecounts$GeneType == "protein_coding") # select only protein coding genes
xiecounts = xiecounts[, c(1,2, 26:73)] # get rid of unwanted samples
xie <- c(
  "DP1_ES_1",
  "DP1_DE_1",
  "DP1_GT_1",
  "DP1_PF_1",
  "DP1_PE_1",
  "DP1_ES_2",
  "DP1_DE_2",
  "DP1_GT_2",
  "DP1_PF_2",
  "DP1_PE_2",
  "DP3_ES",
  "DP3_DE",
  "DP3_GT",
  "DP3_PF",
  "DP3_PE",
  "DP3-4-5_late_PE_1",
  "DP3-4-5_late_PE_2",
  "DP3-4-5_late_PE_3",
  "DP3-4-5_poly_1",
  "DP3-4-5_poly_2",
  "DP3-4-5_poly_3",
  "E2147_in-vivo_matured_1",
  "E2147_in-vivo_matured_2",
  "E2182_in-vivo_matured_1",
  "E2182_in-vivo_matured_2"
) 
xiecounts <- xiecounts[, xie]  #sort (25 samples)
Xie_stages = c("ES",
               "DE",
               "PGT",
               "PFG",
               "PE",
               "late_PE",
               "polyhormonal",
               "matured_in_vivo")
Xie_donors = c(rep("Xie1", 10), rep("Xie2", 5), rep("Xie3", 6), rep("Xie4", 4))

##HNF4A (Ours)
hnf4a <- read.delim("~/OneDrive/oxford/summer_internship/counts/hnf4a_tidyCounts.txt", header = TRUE, row.names = 1, check.names = FALSE)
rownames(hnf4a) = rownames(marta)
hnf4a.donors = c("SB", "EX1", "03A", "04A")  # data comes from four 'donors' (2 indiduals, 4 cell clones)
hnf4a.stage = c("DE", "PE", "BLC") # 3 stages
order_by_stages.hnf4a = function(counts, stage. = hnf4a.stage) {
  nc = counts[, grepl("DE" , names(counts))] 
  
  for (s in stage.[-1])  {
    i = counts[, grepl(s , names(counts))]
    nc = cbind(nc, i)
  }
  
  return(nc)
}
hnf4a = order_by_stages.hnf4a(hnf4a) 
order_by_donors.hnf4a = function(counts, donors. = hnf4a.donors) {
  nc = counts[, grepl("SB" , names(counts))]  
  
  for (s in donors.[-1])  {
    i = counts[, grepl(s , names(counts))]
    nc = cbind(nc, i)
  }
  
  return(nc)
}
hnf4a = order_by_donors.hnf4a(hnf4a) 


#combine data together
c1 <- merge(hnf4a, marta, by = 0, all = TRUE)
rownames(c1) <- c1$Row.names
c1 <- c1[,-1]

combined_commongenes <- merge(c1, xiecounts, by = 0, all = TRUE) 
rownames(combined_commongenes) <- combined_commongenes$Row.names
combined_commongenes <- combined_commongenes[,-1]
combined_commongenes = na.omit(combined_commongenes)   #remove rows that contain NA values, 19316 genes

#design matrix 
samples <- c(rep("SB", 9), rep("EX1", 9), rep("03A", 8), rep("04A", 9),
             rep(marta.donors, 8),
             rep(Xie_donors, 1))
samples = as.factor(samples) 

stages <- c(rep("DE", 3), rep("PE", 3), rep("BLC",3),
            rep("DE", 3), rep("PE", 3), rep("BLC",3),
            rep("DE", 2), rep("PE", 3), rep("BLC", 3),
            rep("DE", 3), rep("PE", 3), rep("BLC", 3),
            rep(marta.stage2, each = 3),
            rep(Xie_stages[1:5], 3),
            rep(Xie_stages[6], 3),
            rep(Xie_stages[7], 3),
            rep(Xie_stages[8], 4))
stages = as.factor(stages)

study <- c(rep("HNF4A", 35), rep("Perez-Alcantara", 24), rep("Xie", 25))
study <- as.factor(study)

design <- model.matrix( ~ stages + study) 

#filter counts
##get protein coding genes only
xie.old <- read.delim("~/OneDrive/oxford/summer_internship/counts/stembancc.xie.gene.counts.tsv", header = TRUE)
xie.old <- subset(xie.old, xie.old$GeneType == "protein_coding")
xie.old.names <- as.character(xie.old$GeneID)
combined_commongenes = combined_commongenes[xie.old.names,]
counts = as.matrix(combined_commongenes)

##remove unexpressed/low expressed genes 
###load the data into edgeR
dge <- DGEList(counts = counts)
isexpr <- rowSums(cpm(counts) > 1) >= 2
dge <- dge[isexpr, , keep.lib.sizes = TRUE] 
filtered_combined_commongenes <- calcNormFactors(dge) # calculate normalisation factors, TMM used by default (Robinson, 2010)

##normalisation 
v <- voom(filtered_combined_commongenes, design, plot =F) # voom normalize the read counts (normalisation factors stored in DGE object)

batch_corrected = removeBatchEffect(v$E, study) # remove batch effects 

#plot PCA
pca <- prcomp(t(batch_corrected)) # calculate

# for labeling x and y axis i.e. PCA1 and PCA2
pca_out <- as.data.frame(pca$x)
percentage <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 2)
percentage <- paste(colnames(pca_out), "(", paste( as.character(percentage), "%", ")", sep=""))

# bind with design matrix for labeling plots 
pca_plot <- as.data.frame(cbind(pca_out, stages, samples, study))

# remove unwanted samples after calculating pcs (for FHS write up)
some.stages <- filter(pca_plot, stages %in% c("iPSC", "DE", "PGT", "PFG", "PE", "EP", "EN", "BLC", "matured_in_vivo"))
some.stages <- filter(some.stages, !(study == "Xie" & stages %in% c("DE", "EN", "EP", "iPSC", "PE", "PFG", "BLC", "PGT")))

# order stages for plot 
some.stages$stages <- factor(some.stages$stages, 
                             levels = c("iPSC", "DE", "PGT", "PFG", "PE", "EP", "EN", "BLC", "matured_in_vivo"))

# simpson's palette but change colours that look too similar
simpsons = c("#FED439FF", "#709AE1FF", "#8A9197FF",
             "#D2AF81FF", "#FD7446FF", "#D5E4A2FF", "#197EC0FF",
             "#46732EFF", "#71D0F5FF")

# plot 
p <- ggplot(some.stages,
            aes(x = PC1, y = PC2, color = stages, shape = study)) +
  geom_point(size = 3, aes(fill = stages), stroke = 1) +
  scale_shape_manual(values = c(16, 17, 18)) +
  scale_color_manual(values = simpsons) +
  xlab(percentage[1]) +
  ylab(percentage[2]) +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.text.x=element_text(colour="black", size=10),
        axis.text.y=element_text(colour="black", size=10),
        axis.title.x=element_text(colour="black", size=15),
        axis.title.y=element_text(colour="black", size=15),
        axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=15))

pdf("pca.rna.seq.allDev.pdf", width = 8, height = 5)
plot(p)
dev.off()
###################################################################

###################################################################
# For each differentiation stage (Figure 2B)
###################################################################
hnf4a = order_by_stages.hnf4a(hnf4a) # order by stages again
de = hnf4a[,1:11]
pe = hnf4a[,12:23]
pe = pe[,-9] # remove outlier
blc = hnf4a[,24:35]
counts_list = list(de = de,
                   pe = pe, 
                   blc = blc)

#design
differentiation = substr(names(hnf4a), 1, 5)
de_diff = differentiation[1:11]
pe_diff = differentiation[12:22]
blc_diff = differentiation[23:34]

de_cell.line = c(rep("SB", 3), rep("EX1", 3), rep("03A", 2), rep("04A", 3))
pe_cell.line = c(rep("SB", 3), rep("EX1", 3), rep("03A", 2), rep("04A", 3))
blc_cell.line = c(rep("SB", 3), rep("EX1", 3), rep("03A", 3), rep("04A", 3))

de_disease = c(rep("0", 3), rep("0", 3), rep("1", 2), rep("1", 3))
pe_disease = c(rep("0", 3), rep("0", 3), rep("1", 2), rep("1", 3))
blc_disease = c(rep("0", 3), rep("0", 3), rep("1", 3), rep("1", 3))

de_des = data.frame(differentiation = de_diff,
                    cell.line = de_cell.line,
                    disease = de_disease)
pe_des = data.frame(differentiation = pe_diff,
                    cell.line = pe_cell.line,
                    disease = pe_disease)
blc_des = data.frame(differentiation = blc_diff,
                     cell.line = blc_cell.line,
                     disease = blc_disease)

design_list = list(de = de_des,
                   pe = pe_des,
                   blc = blc_des)

corr_vst = list()
pca = list()
dge = list()
v = list()
for(i in 1:3) {
  
  # filter counts 
  isexpr = rowSums(cpm(counts_list[[i]]) > 1) >= 2 # expressed at least in one sample group e.g. EX1
  counts_list[[i]] = counts_list[[i]][isexpr, ]
  
  # normalisation
  dge[[i]] = DGEList(counts = counts_list[[i]], samples = design_list[[i]])
  dge[[i]] = calcNormFactors(dge[[i]]) # TMM normalisation 
  design = model.matrix(~ disease + differentiation, design_list[[i]])
  v[[i]] = voom(dge[[i]], design, plot = F) # normalise reads taking acct of TMM factors
  
  # batch correction 
  batch = design_list[[i]]$differentiation
  batch_corrected = removeBatchEffect(v[[i]]$E, batch) # remove batch effects 
  
  # store 
  corr_vst[[i]] = batch_corrected
  
  # pca 
  pca[[i]] = prcomp(t(corr_vst[[i]]))
  
}


# get plot data
pca_out = list()
percentage = list()
all = list()
for(i in 1:3) {
  pca_out[[i]] = as.data.frame(pca[[i]]$x)
  percentage[[i]] = round(pca[[i]]$sdev^2 / sum(pca[[i]]$sdev^2) * 100, 2)
  percentage[[i]] = paste( colnames(pca_out[[i]]), "(", paste( as.character(percentage[[i]]), "%", ")", sep="") )
  all[[i]] = cbind(design_list[[i]], pca_out[[i]])
}

theme = theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              strip.background=element_blank(),
              axis.text.x=element_text(colour="black", size=12),
              axis.text.y=element_text(colour="black", size=12),
              axis.title.x=element_text(colour="black", size=15),
              axis.title.y=element_text(colour="black", size=15),
              axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"),
              legend.title=element_text(size=15), 
              legend.text=element_text(size=12))

# plot 
pdf("pca.rna.seq.eachDev.rough.pdf", width = 5, height = 5)
for(i in 1:3) {
  print(ggplot(all[[i]], aes(PC1, PC2, colour = differentiation, shape = cell.line)) +
          geom_point(size=5) +
          scale_shape_manual(values=c(17, 18, 0, 1)) +
          scale_color_simpsons() +
          theme +
          xlab(percentage[[i]][1]) +
          ylab(percentage[[i]][2]) + 
          coord_fixed())
}
dev.off()

# plot of equal size
g = list()
for(i in 1:3) {
  g[[i]] <- ggplot(all[[i]], aes(PC1, PC2, colour = differentiation, shape = cell.line)) +
    geom_point(size=5) +
    scale_shape_manual(values=c(17, 18, 0, 1)) +
    scale_color_simpsons() +
    theme +
    xlab(percentage[[i]][1]) +
    ylab(percentage[[i]][2]) + 
    coord_fixed() + 
    theme(legend.position="none")
}
pdf("pca.rna.seq.eachDev.same.size.pdf", width = 10)
ggarrange(g[[1]], g[[2]], g[[3]], heights = c(2, 2, 2), widths = c(2, 2, 2),
          ncol = 3, nrow = 1, align = "h")
dev.off()

###################################################################
