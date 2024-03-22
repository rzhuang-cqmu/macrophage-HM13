
# Network Preservation
# 
# https://cran.r-project.org/web/packages/NetRep/vignettes/NetRep.html


# settings ----

rm(list=ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(viridis)
  library(RColorBrewer)
  library(stringr)
  library("NetRep")
  library("WGCNA")
  library(GEOquery)
})


setwd()
list.files()

set.seed(444)


# STAGE ----

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSe40231
# Multi-organ expression profiling uncovers a gene module in coronary artery disease 
# involving transendothelial migration of leukocytes and LIM domain binding 2: 
# the Stockholm Atherosclerosis Gene Expression (STAGE) study
# https://pubmed.ncbi.nlm.nih.gov/19997623/
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2780352/

# associated with other publications:
# https://www.nature.com/articles/s44161-021-00009-1
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4855300/
# https://pubmed.ncbi.nlm.nih.gov/27540175/
# https://pubmed.ncbi.nlm.nih.gov/31196451/
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9583458/
# 

## set up the data ----

module_df <- read.csv("data/gene_modules_identifiers.csv", row.names = 1)
head(module_df)
dim(module_df)
length(unique(module_df$gene_id))
module_df <- module_df[!duplicated(module_df$gene_id), ]

module_data <- as.vector(module_df$colors)
names(module_data) <- module_df$gene_id
#write.csv(module_data, file = "data/module_data.csv")

load(file="data/input_mat.RData")
input_mat[1:5,1:10]
dim(input_mat)
dim(na.omit(input_mat))

load("data/eset_GSE40231_all.tissues.RData")
load("data/input_all.RData")

dim(input_all)
dim(na.omit(input_all))

length(unique(colnames(input_mat)))
length(unique(colnames(input_all)))
input_mat <- input_mat[, !duplicated(colnames(input_mat))]
input_all <- input_all[, !duplicated(colnames(input_all))]

all(names(module_data)==colnames(input_mat))
all(names(module_data)==colnames(input_all))

# select categories of interest
pdata <- pData(eset)[, c(1, 2, 6, 8, 10, 32)]
head(pdata)
table(pdata$`tissue:ch1`)
sel.categories <- c("Liver")
pdata2 <- pdata[pdata$`tissue:ch1` %in% sel.categories, ]
data.liver <- input_all[rownames(pdata2), ]
dim(data.liver)
sel.categories <- c("Skeletal muscle")
pdata2 <- pdata[pdata$`tissue:ch1` %in% sel.categories, ]
data.skeletal <- input_all[rownames(pdata2), ]
dim(data.skeletal)
sel.categories <- c("Visceral fat")
pdata2 <- pdata[pdata$`tissue:ch1` %in% sel.categories, ]
data.visceral <- input_all[rownames(pdata2), ]
dim(data.visceral)


## obtain networks ----

adj <- adjacency(input_mat, type = "signed", power = 8, corFnc = "bicor")
adj[1:5,1:10]
cor.data <- WGCNA::cor1(input_mat)
cor.data[1:5,1:10]

data.liver[1:5,1:10]
adj.liver <- adjacency(data.liver, type = "signed", power = 8, corFnc = "bicor")
adj.liver[1:5,1:10]
cor.liver <- WGCNA::cor1(data.liver, use = "all.obs")
cor.liver[1:5,1:10]

data.skeletal[1:5,1:10]
adj.skeletal <- adjacency(data.skeletal, type = "signed", power = 8, corFnc = "bicor")
adj.skeletal[1:5,1:10]
cor.skeletal <- WGCNA::cor1(data.skeletal, use = "all.obs")
cor.skeletal[1:5,1:10]

data.visceral[1:5,1:10]
adj.visceral <- adjacency(data.visceral, type = "signed", power = 8, corFnc = "bicor")
adj.visceral[1:5,1:10]
cor.visceral <- WGCNA::cor1(data.visceral, use = "all.obs")
cor.visceral[1:5,1:10]



#module_data <- as.disk.matrix(x=module_data, file="data/module_data.rds", serialize=TRUE)
input_mat <- as.disk.matrix(x=input_mat, file="data/input_mat.rds", serialize=TRUE)
adj <- as.disk.matrix(x=adj, file="data/adj.rds", serialize=TRUE)
cor.data <- as.disk.matrix(x=cor.data, file="data/cor.data.rds", serialize=TRUE)
data.liver <- as.disk.matrix(x=data.liver, file="data/data.liver.rds", serialize=TRUE)
adj.liver <- as.disk.matrix(x=adj.liver, file="data/adj.liver.rds", serialize=TRUE)
cor.liver <- as.disk.matrix(x=cor.liver, file="data/cor.liver.rds", serialize=TRUE)
data.skeletal <- as.disk.matrix(x=data.skeletal, file="data/data.skeletal.rds", serialize=TRUE)
adj.skeletal <- as.disk.matrix(x=adj.skeletal, file="data/adj.skeletal.rds", serialize=TRUE)
cor.skeletal <- as.disk.matrix(x=cor.skeletal, file="data/cor.skeletal.rds", serialize=TRUE)
data.visceral <- as.disk.matrix(x=data.visceral, file="data/data.visceral.rds", serialize=TRUE)
adj.visceral <- as.disk.matrix(x=adj.visceral, file="data/adj.visceral.rds", serialize=TRUE)
cor.visceral <- as.disk.matrix(x=cor.visceral, file="data/cor.visceral.rds", serialize=TRUE)



## assess preservation ----

rm(list=ls())

module_df <- read.csv("data/gene_modules_identifiers.csv", row.names = 1)
head(module_df)
dim(module_df)
length(unique(module_df$gene_id))
module_df <- module_df[!duplicated(module_df$gene_id), ]

module_data <- as.vector(module_df$colors)
names(module_data) <- module_df$gene_id

#module_data <- attach.disk.matrix("data/module_data.rds")
#module_data <- read.csv(file = "data/module_data.csv")

input_mat <- attach.disk.matrix("data/input_mat.rds")
adj <- attach.disk.matrix("data/adj.rds")
cor.data <- attach.disk.matrix("data/cor.data.rds")

data.liver <- attach.disk.matrix("data/data.liver.rds")
adj.liver <- attach.disk.matrix("data/adj.liver.rds")
cor.liver <- attach.disk.matrix("data/cor.liver.rds")

data.skeletal <- attach.disk.matrix("data/data.skeletal.rds")
adj.skeletal <- attach.disk.matrix("data/adj.skeletal.rds")
cor.skeletal <- attach.disk.matrix("data/cor.skeletal.rds")

data.visceral <- attach.disk.matrix("data/data.visceral.rds")
adj.visceral <- attach.disk.matrix("data/adj.visceral.rds")
cor.visceral <- attach.disk.matrix("data/cor.visceral.rds")


# liver

data_list <- list(arterial=input_mat, liver=data.liver)
correlation_list <- list(arterial=cor.data, liver=cor.liver)
network_list <- list(arterial=adj, liver=adj.liver)

preservation.liver <- NetRep::modulePreservation(
  network=network_list, data=data_list, correlation= correlation_list, 
  selfPreservation = FALSE,
  moduleAssignments=module_data, discovery="arterial", test="liver", 
  nPerm=1000, nThreads=12
)

preservation.liver$observed
preservation.liver$p.values
# Get the maximum permutation test p-value
max_pval <- apply(preservation.liver$p.value, 1, max)
max_pval
#    black        blue       brown       green        grey     magenta        pink         red   turquoise      yellow 
# 0.000999001 0.000999001 0.000999001 1.000000000 1.000000000 1.000000000 0.000999001 0.000999001 0.000999001 1.000000000 

write.csv(as.data.frame(preservation.liver$observed), file = "results/module_preservation/liver_observed.csv")
write.csv(as.data.frame(preservation.liver$p.values), file = "results/module_preservation/liver_pvalues.csv")
save(preservation.liver, file = "data/module_preservation_liver.RData")


# skeletal

data_list <- list(arterial=input_mat, skeletal=data.skeletal)
correlation_list <- list(arterial=cor.data, skeletal=cor.skeletal)
network_list <- list(arterial=adj, skeletal=adj.skeletal)

preservation.skeletal <- NetRep::modulePreservation(
  network=network_list, data=data_list, correlation= correlation_list, 
  selfPreservation = FALSE,
  moduleAssignments=module_data, discovery="arterial", test="skeletal", 
  nPerm=1000, nThreads=12
)

preservation.skeletal$observed
preservation.skeletal$p.values
# Get the maximum permutation test p-value
max_pval <- apply(preservation.skeletal$p.value, 1, max)
max_pval
#   black        blue       brown       green        grey     magenta        pink         red   turquoise      yellow 
# 0.000999001 0.000999001 0.000999001 1.000000000 1.000000000 0.673326673 0.000999001 0.000999001 0.000999001 0.576423576 

write.csv(as.data.frame(preservation.skeletal$observed), file = "results/module_preservation/skeletal_observed.csv")
write.csv(as.data.frame(preservation.skeletal$p.values), file = "results/module_preservation/skeletal_pvalues.csv")
save(preservation.skeletal, file = "data/module_preservation_skeletal.RData")


# visceral

data_list <- list(arterial=input_mat, visceral=data.visceral)
correlation_list <- list(arterial=cor.data, visceral=cor.visceral)
network_list <- list(arterial=adj, visceral=adj.visceral)

preservation.visceral <- NetRep::modulePreservation(
  network=network_list, data=data_list, correlation= correlation_list, 
  selfPreservation = FALSE,
  moduleAssignments=module_data, discovery="arterial", test="visceral", 
  nPerm=1000, nThreads=12
)

preservation.visceral$observed
preservation.visceral$p.values
# Get the maximum permutation test p-value
max_pval <- apply(preservation.visceral$p.value, 1, max)
max_pval
#    black        blue       brown       green        grey     magenta        pink         red   turquoise      yellow 
# 0.000999001 0.000999001 0.000999001 1.000000000 1.000000000 0.702297702 0.031968032 0.000999001 0.000999001 0.001998002 

write.csv(as.data.frame(preservation.visceral$observed), file = "results/module_preservation/visceral_observed.csv")
write.csv(as.data.frame(preservation.visceral$p.values), file = "results/module_preservation/visceral_pvalues.csv")
save(preservation.visceral, file = "data/module_preservation_visceral.RData")


## generate plots ----

rm(list=ls())

module_df <- read.csv("data/gene_modules_identifiers.csv", row.names = 1)
head(module_df)
dim(module_df)
length(unique(module_df$gene_id))
module_df <- module_df[!duplicated(module_df$gene_id), ]

module_data <- as.vector(module_df$colors)
names(module_data) <- module_df$gene_id

input_mat <- attach.disk.matrix("data/input_mat.rds")
adj <- attach.disk.matrix("data/adj.rds")
cor.data <- attach.disk.matrix("data/cor.data.rds")

data.liver <- attach.disk.matrix("data/data.liver.rds")
adj.liver <- attach.disk.matrix("data/adj.liver.rds")
cor.liver <- attach.disk.matrix("data/cor.liver.rds")

data.skeletal <- attach.disk.matrix("data/data.skeletal.rds")
adj.skeletal <- attach.disk.matrix("data/adj.skeletal.rds")
cor.skeletal <- attach.disk.matrix("data/cor.skeletal.rds")

data.visceral <- attach.disk.matrix("data/data.visceral.rds")
adj.visceral <- attach.disk.matrix("data/adj.visceral.rds")
cor.visceral <- attach.disk.matrix("data/cor.visceral.rds")


# liver

data_list <- list(arterial=input_mat, liver=data.liver)
correlation_list <- list(arterial=cor.data, liver=cor.liver)
network_list <- list(arterial=adj, liver=adj.liver)

par(mar=c(10,10,5,10)) 

jpeg(filename = "results/module_preservation/topology_arterial.vs.liver.jpeg" , width = 900, height = 1200)
plotModule(
  data=data_list, correlation=correlation_list, network=network_list, 
  main = "Preserved modules Arterial vs Liver", 
  moduleAssignments=module_data, modules= c("blue","turquoise","brown"),
  discovery="arterial", test="liver"
)
dev.off()

par(mar=c(20,15,5,30)) 

jpeg(filename = "results/module_preservation/correlation_arterial.vs.liver.jpeg" , width = 900, height = 900)
plotCorrelation(
  data=data_list, correlation=correlation_list, network=network_list, 
  main = "Arterial vs Liver", 
  moduleAssignments=module_data, modules=c("blue","turquoise","brown"), discovery="arterial",
  test="liver", symmetric=TRUE, orderModules=FALSE
)
dev.off()


# skeletal

data_list <- list(arterial=input_mat, skeletal=data.skeletal)
correlation_list <- list(arterial=cor.data, skeletal=cor.skeletal)
network_list <- list(arterial=adj, skeletal=adj.skeletal)

par(mar=c(10,10,5,10)) 

jpeg(filename = "results/module_preservation/topology_arterial.vs.skeletal.jpeg" , width = 900, height = 1200)
plotModule(
  data=data_list, correlation=correlation_list, network=network_list, 
  main = "Preserved modules Arterial vs Skeletal",
  moduleAssignments=module_data, modules= c("blue","turquoise","brown"),
  discovery="arterial", test="skeletal"
)
dev.off()

par(mar=c(20,15,5,30)) 

jpeg(filename = "results/module_preservation/correlation_arterial.vs.skeletal.jpeg" , width = 900, height = 900)
plotCorrelation(
  data=data_list, correlation=correlation_list, network=network_list, 
  main = "Arterial vs Skeletal",
  moduleAssignments=module_data, modules=c("blue","turquoise","brown"), discovery="arterial",
  test="skeletal", symmetric=TRUE, orderModules=FALSE
)
dev.off()


# visceral

data_list <- list(arterial=input_mat, visceral=data.visceral)
correlation_list <- list(arterial=cor.data, visceral=cor.visceral)
network_list <- list(arterial=adj, visceral=adj.visceral)

par(mar=c(10,10,5,10)) 

jpeg(filename = "results/module_preservation/topology_arterial.vs.visceral.jpeg" , width = 900, height = 1200)
plotModule(
  data=data_list, correlation=correlation_list, network=network_list, 
  main = "Preserved modules Arterial vs Visceral",
  moduleAssignments=module_data, modules= c("blue","turquoise","brown"),
  discovery="arterial", test="visceral"
)
dev.off()

par(mar=c(20,15,5,30)) 

jpeg(filename = "results/module_preservation/correlation_arterial.vs.visceral.jpeg" , width = 900, height = 900)
plotCorrelation(
  data=data_list, correlation=correlation_list, network=network_list, 
  main = "Arterial vs Visceral",
  moduleAssignments=module_data, modules=c("blue","turquoise","brown"), discovery="arterial",
  test="visceral", symmetric=TRUE, orderModules=FALSE
)
dev.off()



# GSE225650 ----

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225650
# 	Multi-ancestry genetic analysis of gene regulation in coronary artery prioritizes disease risk loci
# 	Integrative single-cell meta-analysis reveals disease-relevant vascular cell states and markers in human atherosclerosis
# 	https://pubmed.ncbi.nlm.nih.gov/37950869/

gse <- getGEO("GSE225650", GSEMatrix = TRUE)
show(gse)
# 0 features, 138 samples
show(pData(phenoData(gse[[1]]))[1:5, c(1, 6, 8)])



# GSE120521 ----

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120521
# RNA-seq of stable and unstable section of human atherosclerotic plaques
# The Human-Specific and Smooth Muscle Cell-Enriched LncRNA SMILR Promotes Proliferation by Regulating Mitotic CENPF mRNA and Drives Cell-Cycle Progression Which Can Be Targeted to Limit Vascular Remodeling
# https://pubmed.ncbi.nlm.nih.gov/31339449/

## load data ----

rm(list=ls())

load(file="data/input_mat.RData")
input_mat[1:5,1:10]
dim(input_mat)
input_mat <- input_mat[, !duplicated(colnames(input_mat))]

gse <- getGEO("GSE120521", GSEMatrix = TRUE)
show(gse)
# 0 features, 8 samples
show(pData(phenoData(gse[[1]])))

pdata <- pData(phenoData(gse[[1]]))[, c(1, 2, 6, 8, 10, 42, 43)]
head(pdata)
table(pdata$source_name_ch1)
table(pdata$characteristics_ch1)
table(pdata$`stability of region:ch1`)

test.data <- read.csv("data/GSE120521_Athero_RNAseq_FPKM.csv")
head(test.data)
test.data <- test.data[!duplicated(test.data$name), ]
rownames(test.data) <- test.data$name
test.data <- test.data[, -1]
test.data <- as.matrix(test.data)
dim(test.data)
test.data <- t(test.data)
table(goodGenes(test.data))
test.data <- test.data[, goodGenes(test.data)]

int1 <- intersect(colnames(test.data), colnames(input_mat))
test.data <- test.data[, int1]
test.data[1:5,1:8]
any(!is.finite(test.data))

dim(input_mat)
input_mat <- input_mat[, int1]
any(!is.finite(input_mat))

all(colnames(input_mat)==colnames(test.data))

module_df <- read.csv("data/gene_modules_identifiers.csv", row.names = 1)
head(module_df)
dim(module_df)
length(unique(module_df$gene_id))
module_df <- module_df[!duplicated(module_df$gene_id), ]
module_data <- as.vector(module_df$colors)
names(module_data) <- module_df$gene_id
module_data <- module_data[int1]

save(module_data, file = "data/module_data_GSE120521.RData")


## obtain networks ----

adj <- adjacency(input_mat, type = "signed", power = 8, corFnc = "bicor")
adj[1:5,1:10]
cor.data <- WGCNA::cor1(input_mat)
cor.data[1:5,1:10]

test.adj <- adjacency(test.data, type = "signed hybrid", power = 8, corFnc = "bicor")
test.adj[1:5,1:10]
any(!is.finite(test.adj))
test.cor <- WGCNA::cor1(test.data)
test.cor[1:5,1:10]
any(!is.finite(test.adj))

disco.data <- as.disk.matrix(x=input_mat, file="data/disco.data_.GSE120521.rds", serialize=TRUE)
disco.adj <- as.disk.matrix(x=adj, file="data/disco.adj_.GSE120521.rds", serialize=TRUE)
disco.cor <- as.disk.matrix(x=cor.data, file="data/disco.cor_.GSE120521.rds", serialize=TRUE)

test.data <- as.disk.matrix(x=test.data, file="data/test.data_.GSE120521.rds", serialize=TRUE)
test.adj <- as.disk.matrix(x=test.adj, file="data/test.adj_.GSE120521.rds", serialize=TRUE)
test.cor <- as.disk.matrix(x=test.cor, file="data/test.cor_.GSE120521.rds", serialize=TRUE)


## assess preservation ----

rm(list=ls())

load(file = "data/module_data_GSE120521.RData")

disco.data <- attach.disk.matrix(file="data/disco.data_.GSE120521.rds")
disco.adj <- attach.disk.matrix(file="data/disco.adj_.GSE120521.rds")
disco.cor <- attach.disk.matrix(file="data/disco.cor_.GSE120521.rds")

test.data <- attach.disk.matrix(file="data/test.data_.GSE120521.rds")
test.adj <- attach.disk.matrix(file="data/test.adj_.GSE120521.rds")
test.cor <- attach.disk.matrix(file="data/test.cor_.GSE120521.rds")

data_list <- list(discovery=disco.data, test=test.data)
correlation_list <- list(discovery=disco.adj, test=test.adj)
network_list <- list(discovery=disco.cor, test=test.cor)

preservation.liver <- NetRep::modulePreservation(
  network=network_list, data=data_list, correlation= correlation_list, 
  selfPreservation = FALSE,
  moduleAssignments=module_data, discovery="discovery", test="test", 
  nPerm=1000, nThreads=12
)

preservation.liver$observed
preservation.liver$p.values
# Get the maximum permutation test p-value
max_pval <- apply(preservation.liver$p.value, 1, max)
max_pval
#   black      blue     brown     green      grey   magenta      pink       red turquoise    yellow 
# 0.7092907 0.9350649 0.7472527 0.8201798 0.9920080 0.9440559 0.4985015 0.9390609 0.7052947 0.8941059 

write.csv(as.data.frame(preservation.liver$observed), file = "results/module_preservation/GSE120521_observed.csv")
write.csv(as.data.frame(preservation.liver$p.values), file = "results/module_preservation/GSE120521_pvalues.csv")
save(preservation.liver, file = "data/module_preservation_GSE120521.RData")


# GSE236610 ----

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE236610
# Coronary Plaque Sampling Reveals Molecular Insights Into Coronary Artery Disease
# https://pubmed.ncbi.nlm.nih.gov/37539553/

rm(list=ls())

load(file="data/input_mat.RData")
input_mat[1:5,1:10]
dim(input_mat)
input_mat <- input_mat[, !duplicated(colnames(input_mat))]

gse <- getGEO("GSE236610", GSEMatrix = TRUE)
show(gse)
# 0 features, 27 samples
show(pData(phenoData(gse[[1]])))

pdata <- pData(phenoData(gse[[1]]))[, c(8, 10, 11)]
head(pdata)
table(pdata$source_name_ch1)
table(pdata$characteristics_ch1)
table(pdata$characteristics_ch1.1)

test.data <- read.csv("data/GSE236610_SMART_seq_27samples_rawcounts_070523.csv")
head(test.data)
test.data <- test.data[!duplicated(test.data$ID), ]
rownames(test.data) <- test.data$ID
test.data <- test.data[, -1]
test.data <- as.matrix(test.data)
dim(test.data)

library(edgeR)
y = DGEList(counts = test.data, samples = pdata)
y$genes <- as.data.frame(rownames(test.data))
library(AnnotationDbi)
library(org.Hs.eg.db)
y$genes$Symbol <- mapIds(org.Hs.eg.db, keys = rownames(test.data), keytype = "ENSEMBL", column = "SYMBOL")
head(y$genes)
logcpm <- cpm(y, normalized.lib.sizes = T, log = T)
rownames(logcpm) <- y$genes$Symbol
head(logcpm)

test.data <- t(logcpm)
table(goodGenes(test.data))
test.data <- test.data[, goodGenes(test.data)]

int1 <- intersect(colnames(test.data), colnames(input_mat))
test.data <- test.data[, int1]
test.data[1:5,1:8]
any(!is.finite(test.data))

dim(input_mat)
input_mat <- input_mat[, int1]
any(!is.finite(input_mat))

all(colnames(input_mat)==colnames(test.data))

module_df <- read.csv("data/gene_modules_identifiers.csv", row.names = 1)
head(module_df)
dim(module_df)
length(unique(module_df$gene_id))
module_df <- module_df[!duplicated(module_df$gene_id), ]
module_data <- as.vector(module_df$colors)
names(module_data) <- module_df$gene_id
module_data <- module_data[int1]

save(module_data, file = "data/module_data_GSE236610.RData")


## obtain networks ----

adj <- adjacency(input_mat, type = "signed", power = 8, corFnc = "bicor")
#adj <- adjacency(input_mat, type = "signed", power = 8, corOptions = list(use = "p", method = "spearman"))
adj[1:5,1:10]
cor.data <- WGCNA::cor1(input_mat)
cor.data[1:5,1:10]

test.adj <- adjacency(test.data, type = "signed", power = 8, corFnc = "bicor")
#test.adj <- adjacency(test.data, type = "signed", power = 8, corOptions = list(use = "p", method = "spearman"))
test.adj[1:5,1:10]
any(!is.finite(test.adj))
test.cor <- WGCNA::cor1(test.data)
test.cor[1:5,1:10]
any(!is.finite(test.adj))

disco.data <- as.disk.matrix(x=input_mat, file="data/disco.data_.GSE236610.rds", serialize=TRUE)
disco.adj <- as.disk.matrix(x=adj, file="data/disco.adj_.GSE236610.rds", serialize=TRUE)
disco.cor <- as.disk.matrix(x=cor.data, file="data/disco.cor_.GSE236610.rds", serialize=TRUE)

test.data <- as.disk.matrix(x=test.data, file="data/test.data_.GSE236610.rds", serialize=TRUE)
test.adj <- as.disk.matrix(x=test.adj, file="data/test.adj_.GSE236610.rds", serialize=TRUE)
test.cor <- as.disk.matrix(x=test.cor, file="data/test.cor_.GSE236610.rds", serialize=TRUE)


## assess preservation ----

rm(list=ls())

load(file = "data/module_data_GSE236610.RData")

disco.data <- attach.disk.matrix(file="data/disco.data_.GSE236610.rds")
disco.adj <- attach.disk.matrix(file="data/disco.adj_.GSE236610.rds")
disco.cor <- attach.disk.matrix(file="data/disco.cor_.GSE236610.rds")

test.data <- attach.disk.matrix(file="data/test.data_.GSE236610.rds")
test.adj <- attach.disk.matrix(file="data/test.adj_.GSE236610.rds")
test.cor <- attach.disk.matrix(file="data/test.cor_.GSE236610.rds")

data_list <- list(discovery=disco.data, test=test.data)
correlation_list <- list(discovery=disco.adj, test=test.adj)
network_list <- list(discovery=disco.cor, test=test.cor)

preservation.liver <- NetRep::modulePreservation(
  network=network_list, data=data_list, correlation= correlation_list, 
  selfPreservation = FALSE,
  moduleAssignments=module_data, discovery="discovery", test="test", 
  nPerm=1000, nThreads=12
)

preservation.liver$observed
preservation.liver$p.values
# Get the maximum permutation test p-value
max_pval <- apply(preservation.liver$p.value, 1, max)
max_pval
# black      blue     brown     green      grey   magenta      pink       red turquoise    yellow 
# 0.3166833 0.9990010 0.7392607 0.9100899 0.5434565 0.8131868 0.7782218 0.5154845 0.9400599 0.9980020 

write.csv(as.data.frame(preservation.liver$observed), file = "results/module_preservation/GSE236610_observed.csv")
write.csv(as.data.frame(preservation.liver$p.values), file = "results/module_preservation/GSE236610_pvalues.csv")
save(preservation.liver, file = "data/module_preservation_GSE236610.RData")



## generate plots ----

rm(list=ls())


load(file = "data/module_data_GSE236610.RData")

disco.data <- attach.disk.matrix(file="data/disco.data_.GSE236610.rds")
disco.adj <- attach.disk.matrix(file="data/disco.adj_.GSE236610.rds")
disco.cor <- attach.disk.matrix(file="data/disco.cor_.GSE236610.rds")

test.data <- attach.disk.matrix(file="data/test.data_.GSE236610.rds")
test.adj <- attach.disk.matrix(file="data/test.adj_.GSE236610.rds")
test.cor <- attach.disk.matrix(file="data/test.cor_.GSE236610.rds")

data_list <- list(discovery=disco.data, test=test.data)
correlation_list <- list(discovery=disco.adj, test=test.adj)
network_list <- list(discovery=disco.cor, test=test.cor)

par(mar=c(10,10,5,10)) 

jpeg(filename = "results/module_preservation/topology_GSE236610.jpeg" , width = 900, height = 1200)
plotModule(
  data=data_list, correlation=correlation_list, network=network_list, 
  main = "Preserved modules GSE236610", 
  moduleAssignments=module_data, modules= c("blue","turquoise","brown"),
  discovery="discovery", test="test"
)
dev.off()

par(mar=c(20,15,5,30)) 

jpeg(filename = "results/module_preservation/correlation_GSE236610.jpeg" , width = 900, height = 900)
plotCorrelation(
  data=data_list, correlation=correlation_list, network=network_list, 
  main = "GSE236610", 
  moduleAssignments=module_data, modules=c("blue","turquoise","brown"), discovery="discovery",
  test="test", symmetric=TRUE, orderModules=FALSE
)
dev.off()



# end ---------------------------------------------------------------------

sessionInfo()
