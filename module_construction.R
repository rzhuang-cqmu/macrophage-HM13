# description -------------------------------------------------------------

# Identification of myeloid cell gene drivers contributing to early atherosclerosis through promoting foam cell expansion

# In this study, we implemented a targeted weighted gene co-expression network analysis
# (WGCNA), a robust algorithm to uncover molecular rewiring correlated with phenotypic/genotypic
# traits of interest, to identify biologically meaningful clusters (modules) of interconnected genes highly
# correlated with the four key drivers in the CAD-causal RGNs (AIP, DRAP1, POLR2I, and PQBP1) in
# CAD patients. WGCNA was performed in transcriptomic datasets obtained from the CAD cohort.

# Series GSE40231:
# Affymetrix transcriptomic microarray data from atherosclerotic aortic root wall samples of 40 CAD patients
# (marked “STANS”) and non-atherosclerotic internal mammary arterial wall (marked “MAMMARY”) from 40 CAD patients

# GPL570	[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array

# Gene-expression values were pre-processed with the robust multichip average (RMA) procedure
# in three steps (background adjustment, quantile normalization, summarization).

# WGCNA to identify the gene expression profiling associated with AIP, DRAP1, POLR2I, or PQBP1
# Functional enrichment and STRING PPI analyses
# Correlations of the Module Hub Gene Signature Scores and Macrophage Plaque Gene Signature Score with the Final Target Gene(s)
# PPI validation for the final targets




# settings ----------------------------------------------------------------

rm(list = ls())

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
})

setwd()
list.files()

set.seed(444)


# Step I: Inspection & Preprocessing GSE40231 ----

# starting from cell-gene readcount matrix; 
# sample identity, corresponding to the 'sample number in processed data' numbers, is appended to each cell barcode column name

gse <- getGEO('GSE40231',GSEMatrix=TRUE)
show(gse)
# 54675 features, 278 samples
show(pData(phenoData(gse[[1]]))[1:5,c(1,6,8)])

pdata <- pData(phenoData(gse[[1]]))[, c(1,2,6,8,10,32)]
head(pdata)
table(pdata$`tissue:ch1`)
sel.categories <- c("Atherosclerotic aortic wall", "Internal mammary artery")
pdata <- pdata[pdata$`tissue:ch1`%in%sel.categories, ]

eset <- gse$GSE40231_series_matrix.txt.gz
eset <- eset[, rownames(pdata)]
# 54675 features, 80 samples

# add gene annotations

featureData(eset)

rownames(eset) %in% keys(hgu133plus2.db) %>%  summary

# function to collapse duplicate entries
collapser <- function(x){
  x %>% unique %>% sort %>% paste(collapse = "|")
}

ae.annots <- AnnotationDbi::select(
  x       = hgu133plus2.db,
  keys    = rownames(eset),
  columns = c("PROBEID", "ENSEMBL", "ENTREZID", "SYMBOL"),
  keytype = "PROBEID"
) %>%
  group_by(PROBEID) %>%
  summarise_each(funs(collapser)) %>%
  ungroup

# check
all(ae.annots$PROBEID == rownames(eset))
# FALSE

setdiff(ae.annots$PROBEID, rownames(eset))
setdiff(rownames(eset), ae.annots$PROBEID)

head(ae.annots)
head(rownames(eset))

tail(ae.annots)
tail(rownames(eset))

table(ae.annots$PROBEID == rownames(eset))
# FALSE  TRUE 
# 1118 53557

# re-order
fd <- data.frame(ae.annots, stringsAsFactors = F)
rownames(fd) <- fd$PROBEID
fd <- fd[rownames(eset), ]
head(fd)
table(fd$PROBEID == rownames(eset))
# all TRUE

# add to eset
fd <- new("AnnotatedDataFrame",
          data = fd
)

featureData(eset) <- fd


# check expression distribution

boxplot(log2(exprs(eset)[1:1000, ]), las = 2)
# count data was RMA-normalized

save(eset, file="data/eset_GSE40231.RData")


# Step II: WGCNA ----

rm(list = ls())

load("data/eset_GSE40231.RData")

pdata <- pData(eset)
datTraits <- pdata$`tissue:ch1`
table(datTraits)
# Atherosclerotic aortic wall     Internal mammary artery
#            40                          40

exprs <- exprs(eset)
rownames(exprs) <- featureData(eset)$SYMBOL
# alternatively run WGCNA with array identifiers and change to symbol with:
# #colnames(input_mat) <- mapIds(hgu133plus2.db, keys = colnames(input_mat), keytype = "PROBEID", column = "SYMBOL")
head(exprs, 20)
# remove probes without assigned gene symbol
exprs <- exprs[!rownames(exprs) == "", ]
dim(exprs) # 44662

input_mat <- log2(t(exprs))
input_mat[1:5, 1:10]

# Identification of outlying samples

# check by clustering
plotClusterTreeSamples(datExpr = input_mat)

# determine the mean expression per array and the number of missing values per array
meanExpressionByArray <- apply(input_mat, 1, mean, na.rm = T)
NumberMissingByArray <- apply(is.na(data.frame(input_mat)), 1, sum)
barplot(meanExpressionByArray,
        xlab = "Sample", ylab = "Mean expression",
        main = "Mean expression across samples",
        names.arg = c(1:80), cex.names = 0.7
)

# no obvious outliers

# Keep only arrays containing less than 500 missing entries
KeepArray <- NumberMissingByArray < 500
table(KeepArray)

# all arrays are OK

# Handling missing data and zero variance in probe profiles
NumberMissingByGene <- apply(is.na(data.frame(input_mat)), 2, sum)
summary(NumberMissingByGene)

# Calculate the variances of the probes and the number of present entries
variancedatExpr <- as.vector(apply(as.matrix(input_mat), 2, var, na.rm = T))
no.presentdatExpr <- as.vector(apply(!is.na(as.matrix(input_mat)), 2, sum))
table(no.presentdatExpr)

# no missings in this dataset

# Keep only genes whose variance is non-zero and have at least 4 present entries
KeepGenes <- variancedatExpr > 0 & no.presentdatExpr >= 4

par(mar = c(5, 5, 5, 5))
hist(variancedatExpr, breaks = 100)
summary(variancedatExpr)
quantile(variancedatExpr, 0.75)

input_mat <- input_mat[, variancedatExpr > quantile(variancedatExpr, 0.80)]
# input_mat <- input_mat[, variancedatExpr > quantile(variancedatExpr, 0.95)] # shorter genelist for parameters tuning
dim(input_mat)
# 80 8933 # 80%
# duplicate symbols are kept for now as different probes may have different behaviour


# select the same features for the full dataset
load("data/eset_GSE40231_all.tissues.RData")
pdata <- pData(eset)
datTraits <- pdata$`tissue:ch1`
table(datTraits)
exprs <- exprs(eset)
rownames(exprs) <- featureData(eset)$SYMBOL
head(exprs, 20)
exprs <- exprs[!rownames(exprs) == "", ]
dim(exprs) # 44662
input_all <- log2(t(exprs))
input_all <- input_all[, variancedatExpr > quantile(variancedatExpr, 0.80)]
input_all[1:5, 1:10]
dim(input_all)
# 278 8933 # 80%
save(input_all, file="data/input_all.RData")


# check for genes of interest

colnames(input_mat)[grep("AIP", colnames(input_mat))]
colnames(input_mat)[grep("DRAP1", colnames(input_mat))]
colnames(input_mat)[grep("POLR2I", colnames(input_mat))]
colnames(input_mat)[grep("PQBP1", colnames(input_mat))]

# save(input_mat, file="data/input_mat.RData")


# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(
  input_mat, # <= Input data
  # blockSize = 30,
  powerVector = powers,
  verbose = 5
)

par(mfrow = c(1, 2), mar = c(5, 5, 5, 5))
cex1 <- 0.9
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence") # , ylim = c(0,1)
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.80, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red"
)

# Soft threshold power was chosen to be the smallest value such that approximate scale-free topology of R2 > 0.8 was achieved
picked_power <- 8
# A soft-threshold power of 8 was used to achieve approximate scale-free topology (R2>0.8).

allowWGCNAThreads() # allow multi-threading (optional)
enableWGCNAThreads()

temp_cor <- cor
cor <- WGCNA::cor # Force it to use WGCNA cor function (fix a namespace conflict issue)


## loop through alternative parameters ----

# corType.options <- c("pearson", "bicor")
# power.options <- c(10, 12, 14)
deepSplit.options <- c(0, 2, 4)
minModuleSize.options <- c(50, 100, 200)
pamStage.options <- c(TRUE, FALSE)
mergeCutHeight.options <- c(0.15, 0.25)

number.of.combinations <- c(length(deepSplit.options) * length(minModuleSize.options) * length(pamStage.options) * length(mergeCutHeight.options))
# 36

par(mfrow = c(number.of.combinations, 1), mar = c(0, 20, 0, 0))

for (s in deepSplit.options) {
  for (m in minModuleSize.options) {
    for (p in pamStage.options) {
      for (h in mergeCutHeight.options) {
        netwk <- blockwiseModules(
          input_mat,
          maxBlockSize = 15000,
          randomSeed = 54321,
          loadTOM = FALSE,
          corType = "bicor", # compare 'pearson' to 'bicor' (robust)
          power = picked_power,
          networkType = "signed", # consider sign of correlations in downstream analyses!
          saveTOMs = F,
          deepSplit = s, # compare from less (0) to most sensitive (4)
          minModuleSize = m, # compare with 50, 100 and 200
          pamStage = p, # compare T vs F
          mergeCutHeight = h, # compare with 0.15
          
          numericLabels = T,
          nThreads = 12,
          verbose = 3, indent = 0
        )
        
        cor <- temp_cor # Return cor function to original namespace
        
        # Convert labels to colors for plotting
        mergedColors <- labels2colors(netwk$colors)
        plotColorUnderTree(netwk$dendrograms[[1]],
                           colors = mergedColors[netwk$blockGenes[[1]]],
                           paste0("Split.", s, "_", "Size.", m, "_", "PAM.", p, "_", "Height", h)
        )
        
        rm(netwk)
      }
    }
  }
}


## re-sampling ----

# function sampledBlockwiseModules does not account for balanced sampling
# between the two conditions

rm(list = ls())

load(file = "data/input_mat.RData")
input_mat[1:5, 1:10]
dim(input_mat) # 80 8933

load("data/eset_GSE40231.RData")
pdata <- pData(eset)
colnames(pdata)
table(pdata$`tissue:ch1`)

all(rownames(input_mat) == rownames(pdata))

group1 <- pdata[pdata$`tissue:ch1` == "Atherosclerotic aortic wall", ]
group2 <- pdata[pdata$`tissue:ch1` == "Internal mammary artery", ]

picked_power <- 8
allowWGCNAThreads() # allow multi-threading (optional)
enableWGCNAThreads()
temp_cor <- cor
cor <- WGCNA::cor # Force it to use WGCNA cor function (fix a namespace conflict issue)

par(mfrow = c(20, 1), mar = c(0, 8, 0, 0))

for (i in 1:20) {
  samples.sel <- c(
    sample(rownames(group1), 35),
    sample(rownames(group2), 35)
  )
  
  input.sel <- input_mat[samples.sel, ]
  
  netwk <- blockwiseModules(
    input.sel,
    maxBlockSize = 15000,
    randomSeed = 543,
    loadTOM = FALSE,
    corType = "bicor", 
    power = picked_power,
    networkType = "signed", 
    saveTOMs = F,
    deepSplit = 2, # 
    minModuleSize = 100,
    pamStage = F, #
    mergeCutHeight = 0.25, 
    numericLabels = T,
    nThreads = 12,
    verbose = 3, indent = 0
  )
  
  mergedColors <- labels2colors(netwk$colors)
  plotColorUnderTree(netwk$dendrograms[[1]],
                     colors = mergedColors[netwk$blockGenes[[1]]],
                     paste0("Resampling.", i)
  )
  
  rm(netwk)
}





## final network ----

rm(list = ls())

load(file = "data/input_mat.RData")
dim(input_mat) # 80 8933

picked_power <- 8
allowWGCNAThreads() # allow multi-threading (optional)
enableWGCNAThreads()
temp_cor <- cor
cor <- WGCNA::cor # Force it to use WGCNA cor function (fix a namespace conflict issue)

netwk <- blockwiseModules(
  input_mat,
  maxBlockSize = 15000,
  randomSeed = 54321,
  loadTOM = FALSE,
  corType = "bicor",
  power = picked_power,
  networkType = "signed",
  saveTOMs = F,
  deepSplit = 2,
  minModuleSize = 100,
  pamStage = F,
  mergeCutHeight = 0.25,
  numericLabels = T,
  nThreads = 12,
  verbose = 3, indent = 0
)

cor <- temp_cor # Return cor function to original namespace
save(netwk, file = "data/netwk.RData")

table(netwk$colors)
# 0    1    2    3    4    5    6    7    8    9
# 2404 1870 1500 1108  662  532  279  232  191  155

# Convert labels to colors for plotting
mergedColors <- labels2colors(netwk$colors)
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)

# netwk$colors[netwk$blockGenes[[1]]]
# table(netwk$colors)

moduleLabels <- netwk$colors
moduleColors <- labels2colors(netwk$colors)
MEs <- netwk$MEs
geneTree <- netwk$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "data/networkConstruction.RData"
)


# Relate Module (cluster) Assignments to Groups

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

module_df[1:5, ]
table(module_df$colors)
# black      blue     brown     green      grey   magenta      pink       red turquoise    yellow
# 232      1500      1108       532      2404       155       191       279      1870       662

write.csv(module_df, file = "data/gene_modules_identifiers.csv")

module_df <- read.csv("data/gene_modules_identifiers.csv", row.names = 1)
head(module_df)
module_df[grep("AIP", module_df$gene_id), ] # turquoise
module_df[grep("DRAP1", module_df$gene_id), ] # DRAP1 is not on the 20% most variable genes used for WGCNA
module_df[grep("POLR2I", module_df$gene_id), ] # brown
module_df[grep("PQBP1", module_df$gene_id), ] # turquoise



# Step III: Relate genes to traits ----

rm(list=ls())

load("data/eset_GSE40231.RData")
load("data/networkConstruction.RData")

exprs <- exprs(eset)
rownames(exprs) <- featureData(eset)$SYMBOL
head(exprs)
input_mat.long <- log2(t(exprs))
dim(input_mat.long)
input_mat.long[1:5, 1:10]

load("data/input_mat.RData")
input_mat[1:5, 1:10]
dim(input_mat)
# input_mat.long <- input_mat

pdata <- pData(eset)
pdata$AIP <- input_mat.long[, "AIP"] 
pdata$DRAP1 <- input_mat.long[, "DRAP1"] # DRAP1 is not on the 20% most variable genes used for WGCNA
pdata$POLR2I <- input_mat.long[, "POLR2I"] 
pdata$PQBP1 <- input_mat.long[, "PQBP1"]
# head(input_mat.long[, grep("DRAP1", colnames(input_mat.long))])
head(pdata)

colnames(pdata)
allTraits <- pdata[, c(32:36)]
# allTraits <- pdata[, c(32:35)]
head(allTraits)
allTraits$`tissue:ch1` <- gsub("Internal mammary artery", 0, allTraits$`tissue:ch1`)
allTraits$`tissue:ch1` <- gsub("Atherosclerotic aortic wall", 1, allTraits$`tissue:ch1`)
colnames(allTraits)[1] <- "group"
allTraits$group <- as.numeric(allTraits$group)

# Re-cluster samples
sampleTree1 <- hclust(dist(input_mat), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors <- numbers2colors(allTraits, signed = FALSE)
traitColors[, 1] <- gsub("#FFFFFF", "lightblue", traitColors[, 1]) # internal mammary artery
traitColors[, 1] <- gsub("#FF3300", "orange", traitColors[, 1]) # atherosclerotic aortic wall
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree1, traitColors,
                    groupLabels = names(allTraits), abCol = "blue",
                    main = "Sample dendrogram and trait heatmap"
)


####

datTraits <- allTraits$group

# Define numbers of genes and samples
nGenes <- ncol(input_mat)
nSamples <- nrow(input_mat)
# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(input_mat, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, allTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
# Will display correlations and their p-values
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")",
                    sep = ""
)
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(5, 15, 5, 5), mfrow = c(1, 1))
# Display the correlation values within a heatmap plot
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(allTraits),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 1,
  zlim = c(-1, 1),
  main = paste("Module-trait relationships")
)

# Gene relationship to trait and important modules: Gene Significance and Module Membership
###

input_mat[1:5, 1:10]

group <- as.data.frame(allTraits$group)
names(group) <- "group"
head(group)

modNames <- substring(names(MEs), 3)

geneModuleMembership <- as.data.frame(cor(input_mat, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")
#write.csv(geneModuleMembership, file = "data/geneModuleMembership.csv")
#write.csv(MMPvalue, file = "data/MMPvalue.csv")

geneTraitSignificance <- as.data.frame(cor(input_mat, group, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(group), sep = "")
names(GSPvalue) <- paste("p.GS.", names(group), sep = "")
#write.csv(geneTraitSignificance, file = "geneTraitSignificance.csv")
#write.csv(GSPvalue, file = "GSPvalue.csv")

# correlation and p value for genes of interest:
# AIP: 0.19278427526075	0.086659447793442
# DRAP1: not on the modules
# POLR2I: -0.483531893692713	5.53426529393598E-06
# PQBP1: -0.096943375239457	0.392297306612912


#

module <- "turquoise"
module <- "brown"
module <- "blue"
module <- "red"
module <- "yellow"
module <- "tan"

column <- match(module, modNames)
moduleGenes <- moduleColors == module
par(mfrow = c(1, 3), mar = c(10, 5, 5, 5))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for main trait",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module
)

# par(mfrow = c(3,4), mar=c(5,5,5,5))
par(mfrow = c(4, 3), mar = c(5, 5, 5, 5))

for (i in 1:length(modNames)) {
  module <- modNames[i]
  column <- match(module, modNames)
  moduleGenes <- moduleColors == module
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for main trait",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module
  )
}








# end ---------------------------------------------------------------------
sessionInfo()
