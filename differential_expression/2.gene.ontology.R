library(clusterProfiler)
library(org.Hs.eg.db)

# create gene universe of all tested genes with hnf4a motif
universe = as.character(all.fimo.unique$seq.name) #get entrezid of all genes with hnf4a motif
de.universe = subset(de.all, de.all$entrez.de %in% universe) #get all TESTED genes with hnf4a motif along with diff exp results 
pe.universe = subset(pe.all, pe.all$entrez.pe %in% universe)
blc.universe = subset(blc.all, blc.all$entrez.blc %in% universe)

universe.list <- list(de = universe,
                      pe = universe,
                      blc = universe) #all genes with hnf4a motif
# upregulated genes 
upreg.list <- list(de = de.hnf4a.t.up,
                   pe = pe.hnf4a.t.up,
                   blc = blc.hnf4a.t.up) #first rename to make colnames the same 
go.up = list()
for(i in 1:3) {
  resul <- enrichGO(gene = upreg.list[[i]]$entrez,
                    universe = universe.list[[i]]$seq.name,
                    OrgDb = org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "BP", 
                    minGSSize = 1, 
                    pAdjustMethod = "BH",
                    pvalueCutoff = 1,
                    qvalueCutoff = 1)
  go.up[[length(go.up)+1]] = resul
} 

# downregulated genes 
downreg.list <- list(de = de.hnf4a.t.down,
                     pe = pe.hnf4a.t.down,
                     blc = blc.hnf4a.t.down) #first rename to make colnames the same 
go.down = list()
for(i in 1:3) {
  resul <- enrichGO(gene = downreg.list[[i]]$entrez,
                    universe = universe.list[[i]]$seq.name,
                    OrgDb = org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "BP", 
                    minGSSize = 1, 
                    pAdjustMethod = "BH",
                    pvalueCutoff = 1,
                    qvalueCutoff = 1)
  go.down[[length(go.down)+1]] = resul
}

# quick visualisation 
emapplot(go.down[[3]])
cnetplot(go.down[[2]])
barplot(go.up[[1]], order = TRUE)

# export
go.up.de = as.data.frame(go.up[[1]])
go.up.pe = as.data.frame(go.up[[2]])
go.up.blc = as.data.frame(go.up[[3]])
write.csv(go.up.de, "go.up.de.csv")
write.csv(go.up.pe, "go.up.pe.csv")
write.csv(go.up.blc, "go.up.blc.csv")

go.down.de = as.data.frame(go.down[[1]])
go.down.pe = as.data.frame(go.down[[2]])
go.down.blc = as.data.frame(go.down[[3]])
write.csv(go.down.de, "go.down.de.csv")
write.csv(go.down.pe, "go.down.pe.csv")
write.csv(go.down.blc, "go.down.blc.csv")

