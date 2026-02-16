#### Script to integrate scRNA-Seq datasets to correct for batch effects ####


## Load core single-cell analysis libraries
library(Seurat)
library(SeuratDisk)

library(ggplot2)
library(gridExtra)


## Load plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)


## Load co-expression network analysis packages
library(WGCNA)
library(hdWGCNA)


## Using the cowplot theme for ggplot
theme_set(theme_cowplot())

setwd("I:/GSE248730_RAW/SCRNA/")


## Get data location
dirs <- list.dirs(path = 'I:/GSE248730_RAW/SCRNA', recursive = F, full.names = F)

for(x in dirs[1:14]){
  name <- gsub('_pbmc_','', x)
  
  cts <- ReadMtx(mtx = paste0('I:/GSE248730_RAW/SCRNA/',x,'/matrix.mtx.gz'),
                 features = paste0('I:/GSE248730_RAW/SCRNA/',x,'/features.tsv.gz'),
                 cells = paste0('I:/GSE248730_RAW/SCRNA/',x,'/barcodes.tsv.gz'),feature.column = 2)
  
  
## Create Seurat object
  assign(name, CreateSeuratObject(counts = cts))}


## Data integration through Seurat merge 
bcg_combined <- merge(capture1,y = c(capture10,capture11, capture12, capture13, capture14,
                                     capture2,capture3,capture4,capture5,capture6,capture7,capture8 ,capture9 ) , add.cell.ids = c(dirs), project = "bcg")

rm(capture1,capture10,capture11, capture12, capture13, capture14,
   capture2,capture3,capture4,capture5,capture6,capture7,capture8 ,capture9,cts)

bcg_combined<-JoinLayers(bcg_combined, assay = "RNA")


## Data integration through Seurat integration anchors RPCA (or CCA)
bcglist <- SplitObject(bcg_combined, split.by = "capture")

features <- SelectIntegrationFeatures(object.list = bcg_combined)


bcglist <- lapply(X = bcglist, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

bcglist.anchors <- FindIntegrationAnchors(object.list = bcglist, anchor.features = features, reduction = "rpca")
bcg_combined <- IntegrateData(anchorset = bcglist.anchors)


cellAnnotations <- read.csv("I:/GSE248730_RAW/SCRNA/all_cells_unfiltered_metadata.csv")


bcg_combined$orig.ident <- cellAnnotations$orig.ident
bcg_combined$nCount_RNA=cellAnnotations$nCount_RNA
bcg_combined$nFeature_RNA=cellAnnotations$nFeature_RNA
bcg_combined$batchID=cellAnnotations$batchID
bcg_combined$percent.mt=cellAnnotations$percent.mt
bcg_combined$demuxlet_status=cellAnnotations$demuxlet_status

##
bcg_combined$demuxlet_status=cellAnnotations$gender
bcg_combined$demuxlet_status=cellAnnotations$bcg_status # (placebo/bcg)
bcg_combined$demuxlet_status=cellAnnotations$bcg_td # (t0/t90)


##
bcg_combined$mitoPercent <- PercentageFeatureSet(bcg_combined, pattern='^MT-')
bcg_combined$percent_ribo  <- PercentageFeatureSet(bcg_combined, pattern='^RP[SL]')
bcg_combined$percent_hb <- PercentageFeatureSet(bcg_combined, pattern='^HB[^(P)]')


## Explore QC
bcg_combined$log10GenesPerUMI <- log10(bcg_combined$nFeature_RNA) / log10(bcg_combined$nCount_RNA)

bcg_combined_filtered <- subset(bcg_combined, subset = nCount_RNA > 500 &
                                  nFeature_RNA > 250 &   nFeature_RNA <6000 &
                                  mitoPercent < 20& log10GenesPerUMI > 0.80&
                                  percent_ribo > 1 )# can use demuxlet status 

 
## Normalization using CLR and visualization
bcg_combined_filtered <- NormalizeData(object = bcg_combined_filtered,normalization.method ='CLR')
bcg_combined_filtered <- FindVariableFeatures(object = bcg_combined_filtered)
bcg_combined_filtered <- ScaleData(object = bcg_combined_filtered)


## Normalization using SC transform
bcg_combined_filtered <- SCTransform(bcg_combined_filtered, vars.to.regress = "percent.mt", verbose = FALSE)


## Integration using Harmony
bcg_combined_filtered <- RunHarmony(bcg_combined_filtered,
                                       group.by.vars = c("sample", "batchID"),
                                       reduction = "pca", reduction.save = "harmony")


## Cell type identification
library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()

library(scRNAseq)
counts_matrix=counts_matrix = as.matrix(bcg_combined_filtered@assays$RNA@counts[,colnames(bcg_combined_filtered)]) 

pred.hesc <- SingleR(test = counts_matrix, ref = hpca.se, assay.type.test=1,labels = hpca.se$label.main)
bcg_combined_filtered[["SingleR.labels"]] <- pred.hesc$labels


## Dimension reduction
bcg_combined_filtered <- RunPCA(object = bcg_combined_filtered)
bcg_combined_filtered <- FindNeighbors(object = bcg_combined_filtered, dims = 1:20)
bcg_combined_filtered <- FindClusters(object = bcg_combined_filtered,resolution = 0.1)
bcg_combined_filtered <- RunUMAP(object = bcg_combined_filtered, dims = 1:20)


## Visualization
p1 <- DimPlot(bcg_combined_filtered, reduction = 'umap', group.by = 'seurat_clusters',raster=FALSE)
p2 <- DimPlot(bcg_combined_filtered, reduction = 'umap', group.by = 'orig.ident',raster=FALSE)

p3 <- DimPlot(bcg_combined_filtered, reduction = 'umap', group.by = 'gender',raster=FALSE)
p4 <- DimPlot(bcg_combined_filtered, reduction = 'umap', group.by = 'bcg_status',raster=FALSE)
p5 <- DimPlot(bcg_combined_filtered, reduction = 'umap', group.by = 'bcg_td',raster=FALSE)

p6 <- DimPlot(bcg_combined_filtered, reduction = 'umap', group.by = 'SingleR.labels',raster=FALSE)


## Cell marker visualization
features <- c("LYZ", "CCL5", "IL32", "PTPRCAP", "FCGR3A", "PF4") ###(example markers)
p7<-DotPlot(bcg_combined_filtered, features = features) + RotatedAxis()


## Save and load
SaveH5Seurat(bcg_combined_filtered,'I:/GSE248730_RAW/New folder/bcg_combined_filtered.h5Seurat' ,overwrite = TRUE)

bcg_combined_filtered <- LoadH5Seurat("I:/GSE248730_RAW/New folder/bcg_combined_filtered.h5Seurat")



######################## This completes the initial pre processing and sanity checks for the scRNA-seq data #######################


