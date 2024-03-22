
# description -------------------------------------------------------------

# Final network

# WGCNA to identify the gene expression profiling associated with AIP, DRAP1, POLR2I, or PQBP1
# Functional enrichment and STRING PPI analyses
# Correlations of the Module Hub Gene Signature Scores and Macrophage Plaque Gene Signature Score with the Final Target Gene(s)
# PPI validation for the final targets



# settings ----------------------------------------------------------------

rm(list=ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(viridis)
  library(RColorBrewer)
  library(xlsx)
  library(stringr)
  library(GGally)
  library(GEOquery)
  library(WGCNA)
  library(tidyr)
  library(hgu133plus2.db)
  library(STRINGdb)
  library(ggpubr)
  library(gridExtra)
  library(pathfindR)
  library(singscore)
  library(NMF)
  library(eulerr)
  library(OmnipathR)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(DOSE)
  library("pathview")
  library(ReactomePA)
  
})

setwd()
list.files()

set.seed(444)


# Step IV: Functional enrichment ----

# AIP was found to be positively associated with the clinical phenotye
# in turn, selected modules (containing the genes of interest) were positively (turquoise) or negatively (brown) correlated with phenotype

modules <- read.csv("data/gene_modules_identifiers.csv", row.names = 1)
table(modules$colors)
head(modules)

module.list <- unstack(modules)
lapply(module.list, length)
names(module.list)

load("data/eset_GSE40231.RData")
exprs <- exprs(eset)
head(exprs)
universe <- unique(mapIds(hgu133plus2.db, keys= rownames(exprs), column= "ENTREZID", keytype='PROBEID', multiVals = "first"))
length(universe) # 21368

# ORA

for(i in 1:length(names(module.list))){
  
  genelist <- unique(mapIds(org.Hs.eg.db, keys= unique(module.list[[i]]), 
                            column='ENTREZID', keytype='SYMBOL', multiVals = "first")) # keep only the first gene symbol for pathway analyses
  
  # Gene Ontology over-representation (Biological Process)
  
  ora.res <- enrichGO(gene = genelist,  # Entrez Gene IDs
                      OrgDb = "org.Hs.eg.db",
                      ont = "BP", # using Biological Process terms
                      pAdjustMethod = "BH",
                      pvalueCutoff = 1,
                      qvalueCutoff = 1,
                      universe = universe, # using all genes in the dataset as universe
                      minGSSize = 10,
                      maxGSSize = 500)
  dotplot(ora.res, showCategory = 20, 
          title = paste0("Gene Ontology: Biological Process\nModule ", names(module.list)[i])
          )
  ggsave(filename = paste0("results/functional_analyses/ORA/dotplot_GO_", names(module.list)[i], ".pdf"))
  
  edox <- setReadable(ora.res, 'org.Hs.eg.db', 'ENTREZID')
  cnetplot(edox, categorySize="pvalue", foldChange=NULL, layout = "kk") +
    ggtitle(paste0("Gene-Concept Network Module ", names(module.list)[i], "\nGene Ontology : Biological Process"))
  ggsave(filename = paste0("results/functional_analyses/ORA/cnetplot_GO_", names(module.list)[i], ".pdf"))
  
  edox2 <- pairwise_termsim(edox)
  treeplot(edox2) + ggtitle(paste0("Gene Ontology: Biological Process\n", names(module.list)[i]))
  ggsave(filename = paste0("results/functional_analyses/ORA/treeplot_GO_", names(module.list)[i], ".pdf"))
  
  results.table <- as.data.frame(edox)
  write.csv(results.table, file = paste0("results/functional_analyses/ORA/Enrichment_GO.BP_", names(module.list)[i], ".csv"))
  
  rm(ora.res, edox, edox2, results.table)
  
  # Pathway over-representation (KEGG)
  
  ora.res <- enrichKEGG(gene         = genelist,
                        organism     = 'hsa',
                        pvalueCutoff = 1,
                        qvalueCutoff  = 1)
  dotplot(ora.res, showCategory = 20, 
          title = paste0("KEGG Pathways\nModule ", names(module.list)[i])
  )
  ggsave(filename = paste0("results/functional_analyses/ORA/dotplot_KEGG_", names(module.list)[i], ".pdf"))

  edox <- setReadable(ora.res, 'org.Hs.eg.db', 'ENTREZID')
  cnetplot(edox, categorySize="pvalue", foldChange=NULL, layout = "kk") +
    ggtitle(paste0("Gene-Concept Network Module ", names(module.list)[i], "\nKEGG"))
  ggsave(filename = paste0("results/functional_analyses/ORA/cnetplot_KEGG_", names(module.list)[i], ".pdf"))
  
  edox2 <- pairwise_termsim(edox)
  treeplot(edox2) + ggtitle(paste0("KEGG\n", names(module.list)[i]))
  ggsave(filename = paste0("results/functional_analyses/ORA/treeplot_KEGG_", names(module.list)[i], ".pdf"))
  
  results.table <- as.data.frame(edox)
  write.csv(results.table, file = paste0("results/functional_analyses/ORA/Enrichment_KEGG_", names(module.list)[i], ".csv"))
  
  rm(ora.res, edox, edox2, results.table)

  # Pathway over-representation (Reactome)
  
  ora.res <- enrichPathway(genelist, pvalueCutoff = 1, qvalueCutoff  = 1, readable=TRUE)
  dotplot(ora.res, showCategory = 20, 
          title = paste0("Reactome Pathways\nModule ", names(module.list)[i])
  )
  ggsave(filename = paste0("results/functional_analyses/ORA/dotplot_Reactome_", names(module.list)[i], ".pdf"))
  
  edox <- setReadable(ora.res, 'org.Hs.eg.db', 'ENTREZID')
  cnetplot(edox, categorySize="pvalue", foldChange=NULL, layout = "kk") +
    ggtitle(paste0("Gene-Concept Network Module ", names(module.list)[i], "\nReactome"))
  ggsave(filename = paste0("results/functional_analyses/ORA/cnetplot_Reactome_", names(module.list)[i], ".pdf"))
  
  edox2 <- pairwise_termsim(edox)
  treeplot(edox2) + ggtitle(paste0("Reactome\n", names(module.list)[i]))
  ggsave(filename = paste0("results/functional_analyses/ORA/treeplot_Reactome_", names(module.list)[i], ".pdf"))
  
  results.table <- as.data.frame(edox)
  write.csv(results.table, file = paste0("results/functional_analyses/ORA/Enrichment_Reactome_", names(module.list)[i], ".csv"))

  rm(ora.res, edox, edox2, results.table)
  
}



# Step V: STRING PPI ----

rm(list=ls())

####

load("data/networkConstruction.RData")

GOenr = GOenrichmentAnalysis(moduleColors, moduleLabels, organism = "human", nBestP = 10, ontologies = "BP")
GOenr$bestPTerms$BP$enrichment
tab = GOenr$bestPTerms[[2]]$enrichment
names(tab)
head(tab)
write.table(tab, file = "results/functional_analysis/GO_BP_EnrichmentTable.csv", sep = ",", quote = TRUE, row.names = FALSE)

####


STRINGdb$help("new")
STRINGdb$help("map")
STRINGdb$help("plot_network")

string_db <- STRINGdb$new( version="11", species=9606,
                           score_threshold=400, input_directory="")

geneModuleMembership <- read.csv("data/geneModuleMembership.csv", row.names = 1)
head(geneModuleMembership)

# duplicate symbols have been replaced in geneModuleMembership
# use the ones from "modules" and save for downstream analyses

modules <- read.csv("data/gene_modules_identifiers.csv", row.names = 1)
table(modules$colors)
head(modules)
head(modules$gene_id, 100)
head(rownames(geneModuleMembership), 100)
tail(modules$gene_id, 100)
tail(rownames(geneModuleMembership), 100)

geneModuleMembership$Symbol <- sapply(strsplit(modules$gene_id, "\\|"), `[`, 1)
#write.csv(geneModuleMembership, file = "data/geneModuleMembership.csv")

turquoise <- geneModuleMembership[order(geneModuleMembership$MMturquoise, decreasing = T), ]
#turquoise <- geneModuleMembership[order(geneModuleMembership$MMturquoise, decreasing = F), ] # negative memberships in turquoise are captured by the blue module

# consider unsigned network
#turquoise$MMturquoise <- abs(turquoise$MMturquoise)
#turquoise <- turquoise[order(turquoise$MMturquoise, decreasing = T), ]

turquoise <- turquoise[1:60, ]
head(turquoise)
tail(turquoise)

turquoise <- string_db$map(turquoise, "Symbol", removeUnmappedRows = TRUE )
head(turquoise)
hits <- turquoise$STRING_id[]

par(mfrow=c(1,1))
string_db$plot_network(hits)

#enrichment <- string_db$get_enrichment(hits)
#colnames(enrichment)
#head(enrichment$description, n=10)
#enrichment[1:3,c(1:5,10)]

brown <- geneModuleMembership[order(geneModuleMembership$MMbrown, decreasing = T), ]
brown <- brown[1:60, ]
brown <- string_db$map(brown, "Symbol", removeUnmappedRows = TRUE )
hits <- brown$STRING_id[]
string_db$plot_network(hits)

blue <- geneModuleMembership[order(geneModuleMembership$MMblue, decreasing = T), ]
blue <- blue[1:60, ]
blue <- string_db$map(blue, "Symbol", removeUnmappedRows = TRUE )
hits <- blue$STRING_id[]
string_db$plot_network(hits)



# AIP-interacting network in STRING PPI:
AIP_network <- c("PRAMEF17","PRAMEF10","STAT4","AIP","NUB1","AIPL1","HSP90AA1","HSP90AB1","PTGES3","ARNT")

test <- read.csv("data/gene_modules_identifiers.csv", row.names = 1)
head(test)
test2 <- test[test$gene_id%in%AIP_network, ]
head(test2)

test[grep("PRAMEF", test$gene_id), ]
test[grep("STAT", test$gene_id), ]
test[grep("PTGES", test$gene_id), ]



# Step VI: correlations ----

rm(list=ls())

modules <- read.csv("data/gene_modules_identifiers.csv", row.names = 1)
head(modules)
turquoise <- modules[modules$colors == "turquoise", ]
brown <- modules[modules$colors == "brown", ]
blue <- modules[modules$colors == "blue", ]

load("data/input_mat.RData")
input_mat = as.data.frame(input_mat)
input_mat <- input_mat[, !duplicated(colnames(input_mat))]
dim(input_mat)
input_mat[1:5,1:10]  

AIP_network <- c("PRAMEF17","PRAMEF10","STAT4","NUB1","AIPL1","HSP90AA1","HSP90AB1","PTGES3","ARNT")
AIP_network <- AIP_network[AIP_network%in%colnames(input_mat)]

p <- list()
for(i in 1:length(AIP_network)){
  p[[i]] <- ggscatter(input_mat, x = "AIP", y = AIP_network[i], 
                add = "reg.line", conf.int = TRUE, size = 1, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "AIP", ylab = paste0(AIP_network[i]))
  #jpeg(filename = paste0("results/functional_analyses/AIP_correlation_known.interactor_", AIP_network[i]), width = 480, height = 480, quality = 100)
  #plot(p1)
  #dev.off()
}

do.call("grid.arrange", c(p, ncol=2))

geneModuleMembership <- read.csv("data/geneModuleMembership.csv", row.names = 1)
head(geneModuleMembership)

test <- cor(input_mat)[,"AIP"]
test["AIP"]
test["NUB1"]
test["HSP90AB1"]

head(test)
tail(test)
test <- sort(test, decreasing = T)

int <- intersect(turquoise$gene_id, colnames(input_mat))
p <- list()
for(i in 1:length(int)){
  p[[i]] <- ggscatter(input_mat, x = "AIP", y = int[i], 
                      add = "reg.line", conf.int = TRUE, size = 1, 
                      cor.coef = TRUE, cor.method = "pearson",
                      xlab = "AIP", ylab = paste0(int[i]))
  #jpeg(filename = paste0("AIP_correlation_known.interactor_", AIP_network[i]), width = 480, height = 480, quality = 100)
  #plot(p1)
  #dev.off()
}
do.call("grid.arrange", c(p, ncol=6))

int <- intersect(brown$gene_id, colnames(input_mat))
p <- list()
for(i in 1:length(int)){
  p[[i]] <- ggscatter(input_mat, x = "AIP", y = int[i], 
                      add = "reg.line", conf.int = TRUE, size = 1, 
                      cor.coef = TRUE, cor.method = "pearson",
                      xlab = "AIP", ylab = paste0(int[i]))
  #jpeg(filename = paste0("AIP_correlation_known.interactor_", AIP_network[i]), width = 480, height = 480, quality = 100)
  #plot(p1)
  #dev.off()
}
do.call("grid.arrange", c(p, ncol=6))

int <- intersect(blue$gene_id, colnames(input_mat))
p <- list()
for(i in 1:length(int)){
  p[[i]] <- ggscatter(input_mat, x = "AIP", y = int[i], 
                      add = "reg.line", conf.int = TRUE, size = 1, 
                      cor.coef = TRUE, cor.method = "pearson",
                      xlab = "AIP", ylab = paste0(int[i]))
  #jpeg(filename = paste0("AIP_correlation_known.interactor_", AIP_network[i]), width = 480, height = 480, quality = 100)
  #plot(p1)
  #dev.off()
}
do.call("grid.arrange", c(p, ncol=6))



p1 <- ggscatter(input_mat, x = "AIP", y = "CDC14A", 
          add = "reg.line", conf.int = TRUE, size = 1, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "AIP", ylab = "CDC14A")

p2 <- ggscatter(input_mat, x = "AIP", y = "CACNA2D3", 
          add = "reg.line", conf.int = TRUE, size = 1, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "AIP", ylab = "CACNA2D3")

grid.arrange(p1,p2)



# Step VII: Association with Macrophage Plaque Formation ----

rm(list=ls())

gse <- getGEO('GSE11138',GSEMatrix=TRUE)
show(gse)
# 21888 features, 66 samples
show(pData(phenoData(gse[[1]]))[1:5,c(1,6,8)])

pdata <- pData(phenoData(gse[[1]]))
pdata <- pdata[, c(1,2,8,11:14,89:93)]
head(pdata)

table(pdata$`tissue:ch1`)
table(pdata$`patient number:ch1`)
table(pdata$`plasma lipid/stenosis:ch1`)
table(pdata$`characteristics_ch1.4`)

eset <- gse$GSE11138_series_matrix.txt.gz
eset <- eset[, rownames(pdata)]

# check gene annotations
test <- featureData(eset)
str(test)
featureNames(test)
varLabels(test)
test$Gene_Symbol
test$UniGene_ID

# check expression distribution

boxplot(log2(exprs(eset)[1:1000, ]), las = 2)
# count data was RMA-normalized

save(eset, file="data/eset_GSE11138.RData")


# get gene lists

rm(list = ls())

load(file="data/eset_GSE11138.RData")

geneModuleMembership <- read.csv("data/geneModuleMembership.csv", row.names = 1)
head(geneModuleMembership)

# plaque signature (doi: 10.1161/CIRCGENETICS.111.960773)

plaque.sig <- c("ZNF710", "DENND1A", "WDFY4", "FAM78A", "SYK", "IL4I1", "ZMYND15", "FAM20A",
               "LY86", "CCDC88B", "DISC1", "LILRB3", "NOD2", "RUNX1", "NLRC4", "KIF21B", "ADRBK2", "CPVL", "CORO7")

table(plaque.sig%in%featureData(eset)$Gene_Symbol)
# FALSE  TRUE 
#  13     6

exprs <- exprs(eset)
rownames(exprs) <- featureData(eset)$Gene_Symbol
#rownames(exprs) <- featureData(eset)$UniGene_ID
head(exprs)

intersect(plaque.sig, rownames(exprs))
setdiff(plaque.sig, rownames(exprs))
# "ZNF710"  "DENND1A" "WDFY4"   "FAM78A"  "IL4I1"   "ZMYND15" "FAM20A"  "LY86"    "CCDC88B" "NOD2"    "NLRC4"   "KIF21B"  "CORO7" 
# replacements:
# ZNF10...
# DENND1A by KIAA1608
# WDFY4 by KIAA1607
# FAM78A by FLJ00024
# "IL4I1" ...
# "ZMYND15" by DKFZp434N127
# "FAM20A" by DKFZp434F2322
# "LY86" by MD-1
# "CCDC88B" ...
# "NOD2" by CARD15
# "NLRC4" ...
# "KIF21B" by KIAA0449
# "CORO7" by FLJ22021

"DKFZp547K1113"%in%rownames(exprs)
"FLJ22021"%in%featureData(eset)$Gene_Symbol

plaque.sig <- c("ZNF710", "KIAA1608", "KIAA1607", "FLJ00024", "SYK", "IL4I1", "DKFZp434N127", "DKFZp434F2322",
                "MD-1", "CCDC88B", "DISC1", "LILRB3", "CARD15", "RUNX1", "NLRC4", "KIAA0449", "ADRBK2", "CPVL", "FLJ22021")
int <- intersect(plaque.sig, rownames(exprs))

aheatmap(exprs[int, ])

# singscore macrophage plaque
# rank data
data.sig <- exprs
colnames(data.sig) <- pData(eset)$geo_accession
data.sig[is.na(data.sig)] <- 0
rankData <- rankGenes(data.sig)

# Calculate SingScore
scoredf <- simpleScore(rankData, upSet = int)
head(scoredf)
head(rankData[,2,drop = FALSE])
plotRankDensity(rankData[,1,drop = FALSE], upSet = int, isInteractive = F)
plotDispersion(scoredf,annot = NULL,isInteractive = FALSE)

# calculate score
hist(scoredf$TotalScore)
summary(scoredf$TotalScore)


# singscore modules

turquoise <- geneModuleMembership[order(geneModuleMembership$MMturquoise, decreasing = T), ]
turquoise <- turquoise[1:120, ]
brown <- geneModuleMembership[order(geneModuleMembership$MMbrown, decreasing = T), ]
brown <- brown[1:120, ]
blue <- geneModuleMembership[order(geneModuleMembership$MMblue, decreasing = T), ]
blue <- blue[1:120, ]

# Calculate SingScore
int.turquoise <- intersect(turquoise$Symbol, rownames(exprs))
scoredf.turquoise <- simpleScore(rankData, upSet = int.turquoise)
plotRankDensity(rankData[,1,drop = FALSE], upSet = int.turquoise, isInteractive = F)
plotDispersion(scoredf.turquoise,annot = NULL,isInteractive = FALSE)

int.brown <- intersect(brown$Symbol, rownames(exprs))
scoredf.brown <- simpleScore(rankData, upSet = int.brown)
plotRankDensity(rankData[,1,drop = FALSE], upSet = int.brown, isInteractive = F)
plotDispersion(scoredf.brown,annot = NULL,isInteractive = FALSE)

int.blue <- intersect(blue$Symbol, rownames(exprs))
scoredf.blue <- simpleScore(rankData, upSet = int.blue)
plotRankDensity(rankData[,1,drop = FALSE], upSet = int.blue, isInteractive = F)
plotDispersion(scoredf.blue,annot = NULL,isInteractive = FALSE)


# heatmap

score <- data.frame(scoredf$TotalScore, scoredf.turquoise$TotalScore, scoredf.brown$TotalScore, scoredf.blue$TotalScore)
colnames(score) <- c("Macrophage plaque score", "Turquoise hub score", "Brown hub score", "Blue hub score")
rownames(score) <- colnames(data.sig)
head(score)

write.csv(score, file = "data/all.scores.csv")

par(mfrow=c(1,1))

aheatmap(data.sig[int, ], annCol = score, 
         annColors = list("gray","turquoise","brown","blue"),
         main = "Macrophage plaque gene signature score", sub = "Dataset GSE11138")

hm.turquoise <- t(data.sig[int.turquoise, ])
names.turquoise <- score[order(score$`Turquoise hub score`, decreasing = T), ]
hm.turquoise <- hm.turquoise[rownames(names.turquoise), ]
head(hm.turquoise)
aheatmap(hm.turquoise, annRow = names.turquoise, Rowv = NA, 
         annColors = list("gray","turquoise","brown","blue"),
         main = "Turquoise Hub Gene Score", sub = "Dataset GSE11138")

hm.brown <- t(data.sig[int.brown, ])
names.brown <- score[order(score$`Brown hub score`, decreasing = T), ]
hm.brown <- hm.brown[rownames(names.brown), ]
head(hm.brown)
aheatmap(hm.brown, annRow = names.brown, Rowv = NA, 
         annColors = list("gray","turquoise","brown","blue"),
         main = "Brown Hub Gene Score", sub = "Dataset GSE11138")

hm.blue <- t(data.sig[int.blue, ])
names.blue <- score[order(score$`Blue hub score`, decreasing = T), ]
hm.blue <- hm.blue[rownames(names.blue), ]
head(hm.blue)
aheatmap(hm.blue, annRow = names.blue, Rowv = NA, 
         annColors = list("gray","turquoise","brown","blue"),
         main = "Blue Hub Gene Score", sub = "Dataset GSE11138")

##

head(score)

p1 <- ggscatter(score, x = "Macrophage plaque score", y = "Turquoise hub score", 
          add = "reg.line", conf.int = TRUE, size = 1, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Macrophage plaque score", ylab = "Turquoise hub score")
p2 <- ggscatter(score, x = "Macrophage plaque score", y = "Brown hub score", 
                add = "reg.line", conf.int = TRUE, size = 1, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "Macrophage plaque score", ylab = "Brown hub score")
p3 <- ggscatter(score, x = "Macrophage plaque score", y = "Blue hub score", 
                add = "reg.line", conf.int = TRUE, size = 1, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "Macrophage plaque score", ylab = "Blue hub score")
p4 <- ggscatter(score, x = "Turquoise hub score", y = "Brown hub score", 
                add = "reg.line", conf.int = TRUE, size = 1, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "Turquoise hub score", ylab = "Brown hub score")
p5 <- ggscatter(score, x = "Turquoise hub score", y = "Blue hub score", 
                add = "reg.line", conf.int = TRUE, size = 1, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "Turquoise hub score", ylab = "Blue hub score")
p6 <- ggscatter(score, x = "Blue hub score", y = "Brown hub score", 
                add = "reg.line", conf.int = TRUE, size = 1, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "Blue hub score", ylab = "Brown hub score")

p.list <- list(p1,p2,p3,p4,p5,p6)
do.call("grid.arrange", c(p.list, ncol=3))


##

df.turquoise <- turquoise[rownames(turquoise)%in%int.turquoise, ]
write.csv(df.turquoise, file="data/df.turquoise.csv")

df.brown <- brown[rownames(brown)%in%int.brown, ]
write.csv(df.brown, file="data/df.brown.csv")

df.blue <- blue[rownames(blue)%in%int.blue, ]
write.csv(df.blue, file="data/df.blue.csv")



# Step VIII: PPI for the final targets ----

rm(list=ls())

# Genes selected as final targets met the following criteria:
# "blue" module membership
# negative correlation with AIP
# top of blue membership score (on the 120 genes with highest membership score)
# present in the blue score associated with macrophage plaque formation in GSE11138
# 47 genes fullfill these criteria
# 2 of them are excluded manually: DNAJB12, ACTR2
# 45 genes left:

load(file="data/eset_GSE11138.RData")
exprs <- exprs(eset)
rownames(exprs) <- featureData(eset)$Gene_Symbol
head(exprs)

geneModuleMembership <- read.csv("data/geneModuleMembership.csv", row.names = 1)
head(geneModuleMembership)

blue <- geneModuleMembership[order(geneModuleMembership$MMblue, decreasing = T), ]
blue <- blue[1:120, ]

int.blue <- intersect(blue$Symbol, rownames(exprs))
int.blue <- int.blue[!int.blue %in% c("DNAJB12","ACTR2","MAP3K8")]
write.csv(int.blue, file="data/final_blue_genes.csv")

###


int.blue <- read.csv("data/final_blue_genes.csv", row.names = 1)
int.blue <- int.blue$x
head(int.blue)

get_interaction_resources()


## Strategy 1: all interactors ----

#interactions <- import_omnipath_interactions(resources=c("BioGRID","IntAct", "HPRD", "HuRI"))
interactions <- import_all_interactions(resources=c("BioGRID","IntAct", "HPRD", "HuRI"))
#interactions <- import_pathwayextra_interactions(resources=c("BioGRID","STRING"),organism = 9606)
#interactions = import_omnipath_interactions() %>% as_tibble()
print_interactions(head(interactions))

# test AIP interactions
interactions_AIP <- dplyr::filter(interactions, source_genesymbol == "AIP" |
                                     target_genesymbol == "AIP")
print_interactions(interactions_AIP)

list.interactions <- list()
for(i in 1:length(int.blue)){
  list.interactions[[i]] <- dplyr::filter(interactions, source_genesymbol == int.blue[i] | target_genesymbol == int.blue[i])
}
names(list.interactions) <- int.blue
lapply(list.interactions, dim)
genes.interactions <- unique(c(unlist(lapply(list.interactions, function(x) x$source_genesymbol)),
                          unlist(lapply(list.interactions, function(x) x$target_genesymbol))))
genes.pathfindr <- as.data.frame(genes.interactions)
genes.pathfindr$pval <- 0.05
colnames(genes.pathfindr) <- c("symbol","pval")
head(genes.pathfindr)

genesets <- c("KEGG","Reactome","BioCarta","GO-BP")

for(i in 1:length(genesets)){
    tmp1 <- run_pathfindR(genes.pathfindr,
                          gene_sets = genesets[i],
                          min_gset_size = 10, 
                          max_gset_size = 100,
                          p_val_threshold = 0.05,
                          convert2alias = T,
                          plot_enrichment_chart = F,
                          output_dir = paste0("blue_module_interactors_", genesets[i], "_pathfindR_Results"))
    jpeg(filename = paste0("blue_module_interactors_", genesets[i], "_enrichment_chart.jpeg") , width = 960, height = 480, quality = 100)
    p1 <- enrichment_chart(tmp1)
    plot(p1)
    dev.off()
    jpeg(filename = paste0("blue_module_interactors_", genesets[i],"_term_gene_graph.jpeg") , width = 960, height = 480, quality = 100)
    p2 <- term_gene_graph(tmp1, use_description = T, num_terms = 3)
    plot(p2)
    dev.off()
    rm(tmp1, p1, p2)
}  


string_db <- STRINGdb$new( version="11", species=9606,
                           score_threshold=400, input_directory="")
par(mfrow=c(1,1))
int1 <- as.data.frame(genes.interactions)
colnames(int1) <- "gene_id"
string1 <- string_db$map(int1, "gene_id", removeUnmappedRows = TRUE )
hits <- string1$STRING_id[]
string_db$plot_network(hits)

write.csv(int1, file="results/functional_analyses/PPI/all.interactors.csv")



## Strategy 2: common interactors ----

interactions.biogrid <- import_omnipath_interactions(resources=c("BioGRID"))
interactions.intact <- import_omnipath_interactions(resources=c("IntAct"))
interactions.hprd <- import_omnipath_interactions(resources=c("HPRD"))
interactions.huri <- import_omnipath_interactions(resources=c("HuRI"))

list.biogrid <- list()
for(i in 1:length(int.blue)){
  list.biogrid[[i]] <- dplyr::filter(interactions.biogrid, source_genesymbol == int.blue[i] | target_genesymbol == int.blue[i])
}
names(list.biogrid) <- int.blue
lapply(list.biogrid, dim)
genes.biogrid <- unique(c(unlist(lapply(list.biogrid, function(x) x$source_genesymbol)),
                   unlist(lapply(list.biogrid, function(x) x$target_genesymbol))))

list.intact <- list()
for(i in 1:length(int.blue)){
  list.intact[[i]] <- dplyr::filter(interactions.intact, source_genesymbol == int.blue[i] | target_genesymbol == int.blue[i])
}
names(list.intact) <- int.blue
lapply(list.intact, dim)
genes.intact <- unique(c(unlist(lapply(list.intact, function(x) x$source_genesymbol)),
                         unlist(lapply(list.intact, function(x) x$target_genesymbol))))

list.hprd <- list()
for(i in 1:length(int.blue)){
  list.hprd[[i]] <- dplyr::filter(interactions.hprd, source_genesymbol == int.blue[i] | target_genesymbol == int.blue[i])
}
names(list.hprd) <- int.blue
lapply(list.hprd, dim)
genes.hprd <- unique(c(unlist(lapply(list.hprd, function(x) x$source_genesymbol)),
                       unlist(lapply(list.hprd, function(x) x$target_genesymbol))))

list.huri <- list()
for(i in 1:length(int.blue)){
  list.huri[[i]] <- dplyr::filter(interactions.huri, source_genesymbol == int.blue[i] | target_genesymbol == int.blue[i])
}
names(list.huri) <- int.blue
lapply(list.huri, dim)
genes.huri <- unique(c(unlist(lapply(list.huri, function(x) x$source_genesymbol)), 
                       unlist(lapply(list.huri, function(x) x$target_genesymbol))))

genes.all <- list(genes.biogrid, genes.intact, genes.hprd, genes.huri)
names(genes.all) <- c("Biogrid","Intact","Hprd","Huri")

#par(mfrow=c(1,2))
colors <- c("Steelblue3","Darkgoldenrod1","darkolivegreen2","pink")
plot(eulerr::venn(genes.all),quantities = list(TRUE, cex=1.5), 
     fills = colors, edges = colors,
     labels = names(genes.all), cex = 2)
#plot(euler(genes.all),quantities = list(TRUE, cex=1.5), 
#     fills = colors, edges = colors,
#     labels = names(genes.all), cex = 2)


int2 <- intersect(genes.all$Biogrid, genes.all$Intact)
int2 <- as.data.frame(int2)
colnames(int2) <- "gene_id"
write.csv(int2, file="results/functional_analyses/PPI/common_interactors.csv")

genes.pathfindr <- int2
genes.pathfindr$pval <- 0.05
colnames(genes.pathfindr) <- c("symbol","pval")
head(genes.pathfindr)


for(i in 1:length(genesets)){
  tmp1 <- run_pathfindR(genes.pathfindr,
                        gene_sets = genesets[i],
                        min_gset_size = 10, 
                        max_gset_size = 100,
                        p_val_threshold = 0.05,
                        convert2alias = T,
                        plot_enrichment_chart = F,
                        output_dir = paste0("blue_module_interactors_", genesets[i], "_pathfindR_Results"))
  jpeg(filename = paste0("blue_module_interactors_", genesets[i], "_enrichment_chart.jpeg") , width = 960, height = 480, quality = 100)
  p1 <- enrichment_chart(tmp1)
  plot(p1)
  dev.off()
  jpeg(filename = paste0("blue_module_interactors_", genesets[i],"_term_gene_graph.jpeg") , width = 960, height = 480, quality = 100)
  p2 <- term_gene_graph(tmp1, use_description = T, num_terms = 3)
  plot(p2)
  dev.off()
  rm(tmp1, p1, p2)
}  


string_db <- STRINGdb$new( version="11", species=9606,
                           score_threshold=400, input_directory="")
string2 <- string_db$map(int2, "gene_id", removeUnmappedRows = TRUE )
hits2 <- string2$STRING_id[]
string_db$plot_network(hits2)

# common between Intact and Biogrid
# "LATS1"  "CDK1"   "FBXO11" "YAP1"   "WWTR1"  "BCL6" 



# end ---------------------------------------------------------------------
sessionInfo()

