
# Modular Differential Connectivity (MDC)

# Differential coexpression analysis (as in Zhang et al.)
# https://pubmed.ncbi.nlm.nih.gov/23622250/


# settings ----

rm(list=ls())

library(DGCA)
library(SummarizedExperiment)
library(ggplot2)

setwd()
list.files()

set.seed(444)


# module-based differential correlation analysis ----

modules <- read.csv("data/gene_modules_identifiers.csv", row.names = 1)
table(modules$colors)
head(modules)
dim(modules)

load(file="data/input_mat.RData")
input_mat <- t(input_mat)
input_mat[1:5,1:10]
dim(input_mat)

load("data/eset_GSE40231.RData")
eset
sumexp <- makeSummarizedExperimentFromExpressionSet(eset)
colData(sumexp) <- colData(sumexp)[, "tissue.ch1", drop = F]

table(sumexp$tissue.ch1)
design_mat <- model.matrix(~0+sumexp$tissue.ch1)
head(design_mat)
dim(design_mat)


## module-based differential correlation ----

# This analysis finds the average (median or mean) change in correlation between gene symbols in the two conditions, 
# the significance of that change in correlation, as well as the top genes with a gain and/or 
# loss in correlation with the other genes in the module between the conditions, if any of them are significant.

moduleDC_res = moduleDC(inputMat = input_mat, design = design_mat,
                        compare = c("sumexp$tissue.ch1Atherosclerotic aortic wall", "sumexp$tissue.ch1Internal mammary artery"), 
                        genes = modules$gene_id,
                        labels = modules$colors, nPerm = 100, number_DC_genes = 3,  # 20
                        dCorAvgMethod = "median") # mean
moduleDC_res

#write.csv(moduleDC_res, file = "results/module_characterization/differential_correlation.csv")
write.csv(moduleDC_res, file = "results/module_characterization/differential_correlation_3genes.csv")




# It is also possible to take one module and measure differential correlation strength for each of its genes compared to all of the others in the module:

mod1_genes = modules[modules$ind == "mod1", "values"]
darmanis_mod1 = darmanis[mod1_genes, ]
moduleDC_res = ddcorAll(inputMat = darmanis_mod1, design = design_mat,
                        compare = c("oligodendrocyte", "neuron"), nPerm = 50,
                        getDCorAvg = TRUE, dCorAvgType = "gene_average",
                        dCorAvgMethod = "median")
head(moduleDC_res[["avg_dcor"]])
tail(moduleDC_res[["avg_dcor"]])


## GO enrichment ----

library(GOstats, quietly = T)
library(HGNChelper, quietly = T)
library(org.Hs.eg.db)
moduleGO_res = moduleGO(genes = modules$gene_id, labels = modules$colors,
                        universe = rownames(input_mat), pval_GO_cutoff = 1)
#save(moduleGO_res, file = "data/moduleGO.RData")

moduleGO_df = extractModuleGO(moduleGO_res)
head(moduleGO_df)
write.csv(moduleGO_df, file = "results/module_characterization/moduleGO.csv")

plotModuleGO(moduleGO_df, nTerms = 3, text_size = 12, coord_flip = F)



# end ----
sessionInfo()


