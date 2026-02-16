  ### Scripts

  The folder contains all the code required to reproduce the proposed single-cell RNA-seq analysis pipeline for investigating sex-specific differences in human bone marrow cells.  

- Downloading_RAW_data_count_SRR.bash  
  Bash script for downloading raw sequencing data from the SRA using the SRA Toolkit. This script retrieves FASTQ files corresponding to the SRR accessions associated with the study.  

- cellranger_data.sh  
  Shell script for running Cell Ranger on the downloaded FASTQ files to perform read alignment, barcode processing, and UMI counting, generating gene–cell expression matrices for downstream single-cell analysis.

- Preprocessing_data.R  
  R script for preprocessing single-cell data using Seurat, including data loading, quality control, normalization, batch correction, dimensionality reduction, clustering, and cell-type annotation.  

- Deseq_gene_identification.R  
  R script for identifying differentially expressed genes between male and female bone marrow cells using Seurat’s differential expression framework and visualizing results with volcano plots.  

- hdWGCNA_for_gender_specific_eigene_gene_module_identification.R  
  R script implementing high-dimensional weighted gene co-expression network analysis (hdWGCNA) to identify sex-associated gene modules and hub genes in bone marrow cells.  

- GSVA_gender.R  
  R script for performing Gene Set Variation Analysis (GSVA) to assess sex-specific pathway activity using curated Hallmark gene sets.  

- Psuedotime_to_track_markers_using_longitudinal_data.R  
  R script using Monocle3 to perform trajectory and pseudotime analysis, modeling transcriptional changes across longitudinal time points (t0 vs td90) and relating them to sex-specific effects.

  ---

  ### Analysis_Toolkit

  The text file contains information on all required software tools.

  ---

  ### R_Package_List

  The text file lists all required R packages for the analysis.

  ### Documentation_Final

  The PDF file contains the proposed computational pipeline for investigating sex-specific transcriptional differences in bone marrow scRNA-seq data across vaccination status and longitudinal time points.

  
