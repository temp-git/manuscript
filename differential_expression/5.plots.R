library(VennDiagram)
library(ggplot2)
library(ggsci)
library(pheatmap)
library(RColorBrewer)
##########################################################
# Venn Diagram (Figure 3A-C)
##########################################################

# simple venn diagram to see distribution 
# a fast (but bad plot) tool is venny: http://bioinfogp.cnb.csic.es/tools/venny/ 
# plot using VennDiagram package is better looking 

venn.pal = pal_simpsons("springfield", alpha = 0.5)(3) # palette

# ALL (not just HNF4A targets) deg across development stages 
pdf("all.deg.pdf")
all.deg <- draw.triple.venn(area1 = 33 + 12 + 33 + 11,
                            area2 = 931 + 12 + 33 + 375,
                            area3 = 1163 + 11 + 33 + 375,
                            n12 = 12 + 33,
                            n13 = 11 + 33,
                            n23 = 375 + 33,
                            n123 = 33,
                            category = c("DE", "PE", "BLC"),
                            fill = venn.pal,
                            alpha = 1,
                            fontface = "bold",
                            fontfamily = "sans",
                            cex = 3.5,
                            cat.fontface = "bold",
                            cat.fontfamily = "sans",
                            cat.cex = 4,
                            cat.default.pos = "outer",
                            lty = "blank",
                            overrideTriple = 1) 
dev.off()
grid.newpage()

# HNF4A DEG only
pdf("all.hnf4a.t.deg.pdf")
all.deg <- draw.triple.venn(area1 = 15 + 4 + 13 + 6,
                            area2 = 495 + 4 + 13 + 188,
                            area3 = 587 + 6 + 13 + 188,
                            n12 = 4 + 13,
                            n13 = 6 + 13,
                            n23 = 188 + 13,
                            n123 = 13,
                            category = c("DE", "PE", "BLC"),
                            fill = venn.pal,
                            alpha = 1,
                            fontface = "bold",
                            fontfamily = "sans",
                            cex = 5,
                            cat.fontface = "bold",
                            cat.fontfamily = "sans",
                            cat.cex = 6,
                            cat.default.pos = "outer",
                            lty = "blank",
                            overrideTriple = 1) #all deg across development stages (cut-off LFC = 0)
dev.off()

# HNF4A DEG downregulated
pdf("all.hnf4a.t.deg.down.pdf")
all.deg <- draw.triple.venn(area1 = 7 + 2 + 9 + 4,
                            area2 = 287 + 2 + 9 + 57,
                            area3 = 322 + 4 + 9 + 57,
                            n12 = 2 + 9,
                            n13 = 4 + 9,
                            n23 = 57 + 9,
                            n123 = 9,
                            category = c("DE", "PE", "BLC"),
                            fill = venn.pal,
                            alpha = 1,
                            fontface = "bold",
                            fontfamily = "sans",
                            cex = 5,
                            cat.fontface = "bold",
                            cat.fontfamily = "sans",
                            cat.cex = 6,
                            cat.default.pos = "outer",
                            lty = "blank",
                            overrideTriple = 1) #all deg across development stages (cut-off LFC = 0)
dev.off()

# HNF4A DEG upregulated 
pdf("all.hnf4a.t.deg.up.pdf")
all.deg <- draw.triple.venn(area1 = 8 + 2 + 4 + 2,
                            area2 = 238 + 2 + 4 + 101,
                            area3 = 295 + 2 + 4 + 101,
                            n12 = 2 + 4,
                            n13 = 2 + 4,
                            n23 = 101 + 4,
                            n123 = 4,
                            category = c("DE", "PE", "BLC"),
                            fill = venn.pal,
                            alpha = 1,
                            fontface = "bold",
                            fontfamily = "sans",
                            cex = 5,
                            cat.fontface = "bold",
                            cat.fontfamily = "sans",
                            cat.cex = 6,
                            cat.default.pos = "outer",
                            lty = "blank",
                            overrideTriple = 1) #all deg across development stages (cut-off LFC = 0)
dev.off()
##########################################################

##########################################################
# Heatmap (Figure 3D)
##########################################################
gene_ontology = read.delim("quickgo_pancreatic_dev.tsv") # from quickgo GO:0003323
go_panc_dev = unique(gene_ontology$SYMBOL)

pe_deg = read.csv("~/OneDrive/oxford/summer_internship/differential_gene_expression/dream.deg.pe.07122018.csv")
pe_heat = data.frame(log2FC = pe_deg[pe_deg$GeneSymbol %in% go_panc_dev, 4],
                     row.names = pe_deg[pe_deg$GeneSymbol %in% go_panc_dev, 3])
pe_heat = as.matrix(pe_heat)

blc_deg = read.csv("~/OneDrive/oxford/summer_internship/differential_gene_expression/dream.deg.blc.07122018.csv")
blc_heat = data.frame(log2FC = blc_deg[blc_deg$GeneSymbol %in% go_panc_dev, 4],
                      row.names = blc_deg[blc_deg$GeneSymbol %in% go_panc_dev, 3])
blc_heat = as.matrix(blc_heat)

# combine
combined_heat = merge(pe_heat, blc_heat, by = 0, all = TRUE)
rownames(combined_heat) = combined_heat[,1]
colnames(combined_heat) = c("delete", "PE", "BLC")
combined_heat = combined_heat[,-1]
order = c("HNF4A", "NEUROG3", "NKX2-2", "NKX6-1", "NKX6-2", "NEUROD1", "PAX6", "INSM1", "GIPR", "WNT5A")
combined_heat = combined_heat[match(order, rownames(combined_heat)),]
combined_mat = as.matrix(combined_heat)

# plot matrix 
breaksList = seq(-3, 3, by = 1)
pdf("l2fc.heatmap.pancreatic development.pdf", width = 3.5, height = 5)
pheatmap(combined_mat,
         cluster_rows = FALSE, # retain order of genes as given
         cluster_cols = FALSE, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList,
         display_numbers = TRUE,
         fontsize = 20)
dev.off()

channels = c("SLC2A2", "ABCC8", "ABCC9", "KCNQ1", "KCNH2", "KCNJ8")

pe_heat2 = data.frame(log2FC = pe_deg[pe_deg$GeneSymbol %in% channels, 4],
                      row.names = pe_deg[pe_deg$GeneSymbol %in% channels, 3])
pe_heat2 = as.matrix(pe_heat2)

blc_heat2 = data.frame(log2FC = blc_deg[blc_deg$GeneSymbol %in% channels, 4],
                       row.names = blc_deg[blc_deg$GeneSymbol %in% channels, 3])
blc_heat2 = as.matrix(blc_heat2)

combined_heat2 = merge(pe_heat2, blc_heat2, by = 0, all = TRUE)
rownames(combined_heat2) = combined_heat2[,1]
colnames(combined_heat2) = c("delete", "PE", "BLC")
combined_heat2 = combined_heat2[,-1]
combined_heat2 = combined_heat2[match(channels, rownames(combined_heat2)),]
combined_mat2 = as.matrix(combined_heat2)

breaksList = seq(-3, 3, by = 1)
pdf("l2fc.heatmap.channels.pdf", width = 3.2, height = 3)
pheatmap(combined_mat2,
         cluster_rows = FALSE, # retain order of genes as given
         cluster_cols = FALSE, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList,
         display_numbers = TRUE,
         fontsize = 20)
dev.off()

##########################################################