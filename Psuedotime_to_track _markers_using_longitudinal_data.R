#### Trajectory Analysis Using Monocle3 (t0 vs td90 Data) ####


## Set seed for reproducibility
set.seed(1234)


## Load required libraries
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(monocle3)

library(EnhancedVolcano)
library(ggplot2)
library(tidyverse)


## Load Seurat object
bcg_combined_filtered_pt <- LoadH5Seurat("I:/GSE248730_RAW/New folder/bcg_combined_filtered.h5Seurat")


#### MONOCLE3 WORKFLOW


## Ensure correct UMAP slot naming for Monocle3 compatibility
bcg_combined_filtered_pt[["UMAP"]] <- bcg_combined_filtered_pt[["umap"]]
bcg_combined_filtered_pt[["umap"]] <- NULL


## Convert Seurat object to cell_data_set for Monocle3
cds <- SeuratWrappers::as.cell_data_set(bcg_combined_filtered_pt) #change to cds here


## Cluster cells in Monocle3
cds <- cluster_cells(cds)


## Inspect cell and gene metadata 

# Cell metadata
colData(cds)

# Gene metadata
fData(cds)
rownames(fData(cds))[1:10]

# Add gene short names (required by Monocle3)
fData(cds)$gene_short_name <- rownames(fData(cds))

# Inspect raw counts
counts(cds)


## Assign partitions and clusters from Seurat

# Recreate a single partition
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)


cds@clusters$UMAP$partitions <- reacreate.partition

# Assign Seurat cluster identities
list_cluster <- bcg_combined_filtered_pt@active.ident
cds@clusters$UMAP$clusters <- list_cluster


## Assign UMAP Embeddings
cds@int_colData@listData$reducedDims$UMAP <- bcg_combined_filtered_pt@reductions$UMAP@cell.embeddings


## Learn trajectory graph
cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)


## Visualize clusters before trajectory inference
cluster.before.trajectory <- plot_cells(cds,
                                        color_cells_by = 'cluster',
                                        label_groups_by_cluster = FALSE,
                                        group_label_size = 5) +
  theme(legend.position = "right")

cluster.names <- plot_cells(cds,
                            color_cells_by = "redefined_cluster",
                            label_groups_by_cluster = FALSE,
                            group_label_size = 5) +
  scale_color_manual(values = c('red', 'blue', 'green', 'maroon', 'yellow', 'grey', 'cyan')) +
  theme(legend.position = "right")


## Define root node selection function
get_earliest_principal_node <- function(cds, time_bin="1-10"){
  cell_ids <- which(colData(cds)[, "monocle3_pseudotime"] == time_bin)
  
closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
root_pr_nodes}


## Order cells along Pseudotime -------------------
cds <- order_cells(cds, reduction_method = 'UMAP', root_pr_nodes=get_earliest_principal_node(cds))


## Visualize Pseudotime trajectory
plot_cells(cds,
           color_cells_by = 'pseudotime',
          label_groups_by_cluster = TRUE,
           label_branch_points = TRUE,
          label_roots = TRUE,
          label_leaves = FALSE,
          graph_label_size = 4)


## Extract and store Pseudotime values
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))


## Pseudotime distribution across clusters
ggplot(data.pseudo, aes(monocle3_pseudotime, 
                        reorder(redefined_cluster, 
                                monocle3_pseudotime, median), 
                        fill = redefined_cluster)) +
geom_boxplot()


## Identify genes varying along Pseudotime
deg_bcells <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)

deg_bcells %>% 
  arrange(q_value) %>% 
  filter(status == 'OK') %>% 
  head()


## Visualize example Pseudotime-dependent genes
FeaturePlot(bcg_combined_filtered_pt, features = c('E2F2', 'STMN1', 'CD52'))##exxample features can be used from the prexious analysis


## Transfer Pseudotime back to Seurat object
bcg_combined_filtered_pt$pseudotime <- pseudotime(cds)

Idents(bcg_combined_filtered_pt) <- bcg_combined_filtered_pt$bcg_td


FeaturePlot(bcg_combined_filtered_pt, 
            features = "pseudotime", 
            label = T, raster=FALSE)

