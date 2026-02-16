#### Weighted Gene Co-expression Network Analysis (hdWGCNA) Pipeline ####


## Load core single-cell analysis libraries
library(Seurat)
library(SeuratDisk)


## Load plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)


## Load co-expression network analysis packages
library(WGCNA)
library(hdWGCNA)


## Using the cowplot theme for ggplot
theme_set(theme_cowplot())


## Set random seed for reproducibility
set.seed(12345)


## Load the Seurat object
bcg_combined_filtered <- LoadH5Seurat("I:/GSE248730_RAW/New folder/bcg_combined_filtered.h5Seurat")


## Enable multi threading for WGCNA
enableWGCNAThreads(nThreads = 16)


## Prepare Seurat object for hdWGCNA
seurat_obj <- SetupForWGCNA(
  seurat_obj =bcg_combined_filtered,
  gene_select = "fraction", # the gene selection approach
 fraction = 0.1, # minimum fraction of cells expressing a gene
  wgcna_name = "wgcn1") # hdWGCNA experiment name


## Construct Meta cells by sample
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
 group.by = c('sample'), # specify the columns in seurat_obj@meta.data to group by
 reduction = 'harmony', # select the dimensionality reduction to perform KNN on
 ident.group = 'capture') # set the Idents of the metacell seurat object


## Normalize and scale Meta cell expression matrix
seurat_obj <- NormalizeMetacells(seurat_obj)
seurat_obj <- ScaleData(object = seurat_obj)


## Define expression matrix for network construction
seurat_obj <- SetDatExpr(
  seurat_obj,
 assay = 'RNA', # using RNA assay
 slot = 'scale.data') # using normalized data


## Soft-threshold power selection
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed') # network type


## Visualize soft power results
plot_list <- PlotSoftPowers(seurat_obj)
wrap_plots(plot_list, ncol=2) # sub plotting

power_table <- GetPowerTable(seurat_obj)


## Construct co-expression network
seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = 'hspc_mt') # name of the topological overlap matrix written to disk


## Plot gene dendrogram
PlotDendrogram(seurat_obj, main='hspc hdWGCNA Dendrogram')


## Compute Module Eigengenes 
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="capture")


## Extract module Eigengenes
hMEs <- GetMEs(seurat_obj)
MEs <- GetMEs(seurat_obj, harmonized=FALSE)


## Compute module connectivity (kME)
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'capture'
)


## Rename modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "BCG-M")

modules <- GetModules(seurat_obj) %>% subset(module != 'grey')


## Visualize hub genes by module 
p <- PlotKMEs(seurat_obj, ncol = 5, n_hubs = 5)


### Define traits for module–trait analysis 
cur_traits <- c( 'gender','bcg_status','bcg_td')


## Module–trait correlation analysis
seurat_obj <- ModuleTraitCorrelation(
  seurat_obj,
  traits = cur_traits,
  group.by='capture')

mt_cor <- GetModuleTraitCorrelation(seurat_obj)

names(mt_cor)
head(mt_cor$cor$all_cells[,1:10])


## Visualize module–trait correlations
PlotModuleTraitCorrelation(
  seurat_obj,
  label = 'fdr',
  label_symbol = 'stars',
  text_size = 2,
  text_digits = 2,
 text_color = 'white',
  high_color = 'yellow',
  mid_color = 'black',
 low_color = 'purple',
  plot_max = 0.2,
  combine=FALSE)


## Get the module assignment table
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')


## print the first 6 columns
head(modules[,1:6])


#### Identify module-specific genes
genebyepi<-modules$gene_name[modules$module=='BCG-M1']

