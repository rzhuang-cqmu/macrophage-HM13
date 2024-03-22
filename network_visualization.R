
# settings ----

rm(list=ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(viridis)
  library(RColorBrewer)
  library(stringr)
  library(BioNERO)
  library(SummarizedExperiment)
  library(WGCNA)
  })

setwd()
list.files()

set.seed(123)


# inspection ----

load("data/eset_GSE40231.RData")
eset

load(file="data/input_mat.RData")
dim(input_mat)
input_mat[1:5,1:10]

sumexp <- SummarizedExperiment(assays = t(input_mat),
                               colData = pData(eset)[, "tissue:ch1"]
)
colnames(colData(sumexp)) <- "tissue"

table(sumexp$tissue)
# Atherosclerotic aortic wall     Internal mammary artery 
#            40                          40 

final_exp <- exp_preprocess(sumexp,  Zk_filtering = F,
                            min_exp = 1, 
                            remove_confounders = F,
                            variance_filter = F, percentile = 0.2)
final_exp
dim(final_exp)
#  8933   80

# Heatmap of sample correlations
p <- plot_heatmap(final_exp, type = "samplecor", show_rownames = FALSE)
p

# Heatmap of gene expression (here, only the first 50 genes)
p <- plot_heatmap(
  final_exp[1:50, ], type = "expr", show_rownames = FALSE, show_colnames = FALSE
)
p

# PCA

plot_PCA(final_exp)


save(final_exp, file = "data/final_exp_BioNERO.RData")


# build exp2gcn object from WGCNA ----

load(file="data/input_mat.RData")
dim(input_mat)
colnames(input_mat) <- make.names(colnames(input_mat), unique = T)
input_mat[1:5,1:5]

modules <- read.csv("data/gene_modules_identifiers.csv", row.names = 1)
head(modules)
colnames(modules) <- c("Genes","Modules")
#modules$Modules <- paste0("ME", modules$Modules)
length(unique(modules$Genes))
modules$Genes <- colnames(input_mat)

load("data/networkConstruction.RData")
colnames(MEs) <- paste0("ME", levels(as.factor(modules$Modules)))
MEs[1:5,1:5]

cor_matrix <- WGCNA::bicor(input_mat, maxPOutliers = 0.1)
cor_matrix[1:5,1:5]
adj_matrix <- WGCNA::adjacency.fromSimilarity(
  cor_matrix, power = 8, type = "signed"
)
adj_matrix[1:5,1:5]

#TOM <- WGCNA::TOMsimilarity(adj_matrix, TOMType = "signed")
#Hierarchically cluster genes
#dissTOM <- 1-TOM #hclust takes a distance structure
#dissTOM[1:5,1:5]
#geneTree <- hclust(as.dist(dissTOM), method="average")

kwithin <- WGCNA::intramodularConnectivity(adj_matrix, modules$Modules)
head(kwithin)

params = list(net_type = "signed",
              module_merging_threshold = 0.7,
              SFTpower = 8,
              cor_method = "biweight")

net2 <- list(adjacency_matrix = adj_matrix,
             MEs = MEs, 
             genes_and_modules = modules,
             kIN = kwithin,
             correlation_matrix = cor_matrix,
             params = params,
             dendro_plot_objects = list(
               tree = geneTree,
               Unmerged = modules$Modules,
               Merged = modules$Modules)
)

save(net2, file = "data/BioNERO_net2.RData")


## plot modules ----

rm(list=ls())

load(file = "data/BioNERO_net2.RData")
load(file = "data/final_exp_BioNERO.RData")


# Dendro and colors
plot_dendro_and_colors(net2)
# Eigengene networks
plot_eigengene_network(net2)
plot_ngenes_per_module(net2)


## module stability ----
module_stability(final_exp, net2, nRuns = 10)

## module-trait associations ----
MEtrait <- module_trait_cor(exp = final_exp, MEs = net2$MEs)
head(MEtrait)

plot_module_trait_cor(MEtrait)

## module expression profile ----
plot_expression_profile(
  exp = final_exp, 
  net = net2, 
  plot_module = TRUE, 
  modulename = "turquoise"
)

## Enrichment analysis for conserved protein domains (Interpro) ----
data(zma.interpro)
interpro_enrichment <- module_enrichment(
  net = net2, 
  background_genes = rownames(final_exp),
  annotation = zma.interpro
)
# Print results without geneIDs for better visualization
interpro_enrichment[, -6]

## Hub gene identification ----
hubs <- get_hubs_gcn(final_exp, net2)
head(hubs)
table(hubs$Module)

## extracting subgraphs ----
#edges <- get_edge_list(net2, module="brown")

# Remove edges based on optimal scale-free topology fit
edges_filtered <- get_edge_list(net2, module = "blue", filter = TRUE)
## The correlation threshold that best fits the scale-free topology is 0.9
dim(edges_filtered)
edges_filtered[grep("HM13", edges_filtered$Gene1), ]
#Gene1    Gene2    Weight
#1874644 HM13.1 GPRIN3.1 0.9134916
#1876144 HM13.1  RNASET2 0.9015441
#1948144 HM13.1    TTC28 0.9148915
#2029144 HM13.1 VCPIP1.1 0.9352651
edges_filtered[grep("HM13", edges_filtered$Gene2), ]
#Gene1  Gene2    Weight
#1714622      CCIN HM13.1 0.9155604
#1714628 LINC00260 HM13.1 0.9005821
#1714638     EPB42 HM13.1 0.9368580
#1714819    NUP205 HM13.1 0.9104776
#1714866      NME4 HM13.1 0.9034149
#1714949    KAZN.1 HM13.1 0.9100190
#1714994   DNAJB12 HM13.1 0.9040172
#1715400      GPT2 HM13.1 0.9082824

# similar results with:
# Remove edges based on minimum correlation
edges_filtered <- get_edge_list(
  net2, module = "turquoise", 
  filter = TRUE, method = "min_cor", rcutoff = 0.95
)
dim(edges_filtered)

## Network visualization ----
plot_gcn(
  edgelist_gcn = edges_filtered, 
  net = net2, 
  show_labels = "all",
#  top_n_hubs = 10,
  color_by = "module", 
  hubs = hubs
)

edges_filtered <- get_edge_list(net2, module = "blue", filter = TRUE)
## The correlation threshold that best fits the scale-free topology is 0.9
dim(edges_filtered)

edges_filtered <- get_edge_list(net2, module = c("blue","brown","turquoise"), filter = TRUE)
dim(edges_filtered)

plot_gcn(
  edgelist_gcn = edges_filtered, 
  net = net2, 
  color_by = "module", 
  hubs = hubs
)


#
plot_gcn(
  edgelist_gcn = edges_filtered, 
  net = net2,
  color_by = "module",
  hubs = hubs,
  interactive = TRUE,
  dim_interactive = c(500, 500)
)






hubs <- get_hubs_gcn(t(input_mat), net2)
head(hubs)

edges_filtered <- get_edge_list(net2, module = c("blue","brown","turquoise"), filter = TRUE)
dim(edges_filtered)

edges_filtered <- get_edge_list(net2, module = c("blue"), filter = TRUE)
dim(edges_filtered)

# Remove edges based on minimum correlation
edges_filtered <- get_edge_list(
  net2, module = "turquoise", 
  filter = TRUE, method = "min_cor", rcutoff = 0.95 # method = "optimalSFT" # 
)
dim(edges_filtered)
## [1] 588   3

plot_gcn(
  edgelist_gcn = edges_filtered,
  net2,
  show_labels = "all",
  top_n_hubs = 10,
  color_by = "module",
  hubs = hubs)




# coexpression network with BioNERO ----

rm(list=ls())

load(file = "data/final_exp_BioNERO.RData")

sft <- SFT_fit(final_exp, net_type = "signed", cor_method = "biweight")
sft$power
power <- sft$power
sft$plot

## infer the GCN ----

net <- exp2gcn(
  final_exp, net_type = "signed", SFTpower = 8, 
  module_merging_threshold = 0.7,
  cor_method = "biweight"
)


# Dendro and colors
plot_dendro_and_colors(net)
# Eigengene networks
plot_eigengene_network(net)
plot_ngenes_per_module(net)


## module stability ----
module_stability(final_exp, net, nRuns = 2)

## module-trait associations ----
MEtrait <- module_trait_cor(exp = final_exp, MEs = net$MEs)
head(MEtrait)

plot_module_trait_cor(MEtrait)

## module expression profile ----
plot_expression_profile(
  exp = final_exp, 
  net = net, 
  plot_module = TRUE, 
  modulename = "magenta"
)

## Enrichment analysis for conserved protein domains (Interpro) ----
data(zma.interpro)
interpro_enrichment <- module_enrichment(
  net = net, 
  background_genes = rownames(final_exp),
  annotation = zma.interpro
)
# Print results without geneIDs for better visualization
interpro_enrichment[, -6]

## Hub gene identification ----
hubs <- get_hubs_gcn(final_exp, net)
head(hubs)

## extracting subgraphs ----
edges <- get_edge_list(net, module="magenta")
head(edges)

# Remove edges based on optimal scale-free topology fit
edges_filtered <- get_edge_list(net, module = "magenta", filter = TRUE)
## The correlation threshold that best fits the scale-free topology is 0.7
dim(edges_filtered)
## [1] 588   3

# Remove edges based on p-value
edges_filtered <- get_edge_list(
  net, module = "midnightblue",
  filter = TRUE, method = "pvalue", 
  nSamples = ncol(final_exp)
)
dim(edges_filtered)
## [1] 921   3

# Remove edges based on minimum correlation
edges_filtered <- get_edge_list(
  net, module = "midnightblue", 
  filter = TRUE, method = "min_cor", rcutoff = 0.7
)
dim(edges_filtered)
## [1] 588   3

## Network visualization ----
plot_gcn(
  edgelist_gcn = edges_filtered, 
  net = net, 
  color_by = "module", 
  hubs = hubs
)

plot_gcn(
  edgelist_gcn = edges_filtered, 
  net = net,
  color_by = "module",
  hubs = hubs,
  interactive = TRUE,
  dim_interactive = c(500, 500)
)


## Network statistics ----

net_stats()



# end ----
sessionInfo()
