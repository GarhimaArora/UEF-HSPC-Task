#### Identify differentially expressed (DE) genes ####


## Load libraries
library(Seurat)
library(SeuratDisk)
library(EnhancedVolcano)
library(ggplot2)


## Load the Seurat object
bcg_combined_filtered <- LoadH5Seurat("I:/GSE248730_RAW/New folder/bcg_combined_filtered.h5Seurat")


## Set cell identities (Gender)
b1seurat <- SetIdent(b1seurat, value = b1seurat@meta.data$gender)


## Identify differentially expressed genes
newcells.markers <- FindMarkers(b1seurat, ident.1  = "male", ident.2 = "female")


## Inspect top DE genes
head(newcells.markers)


## Visualize DE genes using volcano plot
EnhancedVolcano(newcells.markers,
                lab = rownames(newcells.markers),
                x = 'avg_log2FC',
                y = 'p_val',pCutoff = 1e-2,
                FCcutoff = 2,
                pointSize = 3.0,
                labSize = 6.0)


## Save plot
ggsave(filename="volcano",plot = last_plot(), path="I:/GSE248730_RAW/New folder/", device = "png", dpi = 500)



