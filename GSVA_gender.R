#### Gene Set Variation Analysis (GSVA): Gender-Specific Pathway Activity ####


## Load gene set & pathway analysis libraries
library(GSVA)
library(msigdbr)


## Load single-cell analysis & visualization libraries
library(Seurat)
library(SeuratDisk)
library(EnhancedVolcano)
library(ggplot2)


## Load Seurat object
bcg_combined_filtered <- LoadH5Seurat("I:/GSE248730_RAW/New folder/bcg_combined_filtered.h5Seurat")


## Extract normalized expression matrix
expr <- as.matrix(GetAssayData(bcg_combined_filtered, slot = "data"))


## Retrieve MSigDB hallmark gene sets
msig <- msigdbr(species = "Homo sapiens", category = "H")
gene_sets <- split(msig$gene_symbol, msig$gs_name)


## Run GSVA on single-cell expression data
gsva_scores <- gsva(expr,
  gene_sets,
  method = "gsva",
  kcdf = "Gaussian",
  verbose = TRUE)


## Add GSVA scores back to Seurat object
seurat_obj[["GSVA"]] <- CreateAssayObject(gsva_scores)
DefaultAssay(seurat_obj) <- "GSVA"


## Visualize pathway activity by gender
FeaturePlot(seurat_obj,
  features = "", # pathway name(s) from previous analysis
  split.by = "gender")


## Differential pathway activity analysis
Idents(seurat_obj) <- "gender"

FindMarkers(seurat_obj,
  assay = "GSVA",
  slot = "data",
  ident.1 = "male",
  ident.2 = "female")


