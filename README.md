# scRNAseq-SSCs
scRNAseq analysis for Paula Cohen lab SSCs project

## INTRODUCTION

## METHOD DETAILS

### single-cell RNA sequencing
Immediately post-sorting, ?? SSCs were run on the 10X Chromium (10X Genomics) and then through library preparation following Chromium Single Cell 3' Reagent Kit (v2 Chemistry). 

## Data analysis
Unsupervised cell clustering analysis was carried out using the Seurat 2.1 R package. Cells with <200 UMIs and genes detected in <3 cells were excluded from the analysis. Gene expression raw counts were normalized following a global-scaling normalization method with a scale factor of 10,000 and a log transformation, using the Seurat NormalizeData function. The top 2000 highly variable genes from C57BL/6J and B6129PF1/J datasets were selected, followed by a canonical correlation analysis (CCA) to identify common sources of variation between the two datasets and minimize batch effect. The first 11 CCA results were chosen for principal component analysis (PCA). Cells where the variance explained by CCA was more than 2-fold (ratio > 0.5) compared to PCA were used for 2-dimensional t-Distributed Stochastic Neighbor Embedding (tSNE) (ref van der maaten and hinton 2008) with 0.8 resolution. Cells contained in cluster 11 (hematopoietic cells) were further subjected to a second round of unsupervised analysis following the same approach, resulting in a tSNE analysis with ~0.1 resolution. The modified Seurat function FindAllMarkers was used to calculate average differential expression among cell clusters. The p-value was calculated using likelihood-ratio test and adjusted by Benjamini-Hochberg method.

## How to use this repository

#### Software Setup
R version 3.4.3 http://cran.us.r-project.org/bin/macosx/R-3.4.3.pkg <br />
dplyr_0.7.4 <br />
Seurat_2.1.0 https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_2.1.0.tar.gz (other Seurat versions will generate slightly different results )<br />

After pulling this repository, create folders **_data_** and **_output_** in the top working folder.
Move Cell Ranger analysis results into **_data_** folder.

### 1. Seurat_setup.R
Unsupervised cell clustering analysis was carried out using the Seurat 2.2 R package. Cells with <500 genes and genes detected within <3 cells were excluded from the analysis. Gene expression raw counts were normalized following a global-scaling normalization method with a scale factor of 10,000 and a log transformation, using the Seurat NormalizeData function. The top 1000 highly variable genes from young C57BL/6J and aged C57BL/6J datasets were selected, followed by canonical correlation analysis (CCA) to identify common sources of variation between the two datasets and minimize the batch effect. The first 20 CCA results were chosen for principal component analysis (PCA). Cells were used for 2-dimensional t-Distributed Stochastic Neighbor Embedding (tSNE) (ref van der maaten and hinton 2008) with 0.8 resolution.

 After running this script, a `mouse_eyes_alignment.Rda` file will be generated inside **_data_** folder.
 Do not modify any files in **_data_** folder.
 
 
### 2. Identify_Cell_Types_Manually.R
All clusters are examed against 122(number may change) CD marker genes.
All cell types are predicted by at least two marker genes with the adjusted p-value(FDR) smaller than 10^-30.

Endothelial cells were identified by Cdh5, Flt1, Kdr, Pecam1, Plvap, Ptprb, and Vwf.<br />
Pericytes were identified by Dcn, Des, Mylk, Pdgfrb, and Rgs5.<br />
Hematopoietic cells were identified by Laptm5 and Ptprc.<br />
Melanocytes were identified by Mlana and Pmel.<br />
Myelinating Schwann cells were identified by Mbp and Mpz.<br />
Retinal pigment epitheliums were identified by Rlbp1 and Rpe65.<br />

`DotPlot` is implemented for visualizing differential expressed marker genes.
Multiple plots and table will be generated, save them if you want. I prefer to keep the original identity of `mouse_eyes_alignment.Rda` intact for further downstream analysis.

### 3. Differential_analysis.R
Modified FindAllMarkers() `FindAllMarkers.UMI()` will generate similar dataframe plus two extra columns `pct.1_UMI` and `pct.2_UMI` to record nUMI. `pct.1_UMI` is nUMI of current cluster, `pct.2_UMI` is average nUMI of rest of clusters.<br />
`FindAllMarkers(object, test.use = "bimod")` : Likelihood-ratio test for single cell gene expression, (McDavid et al., Bioinformatics, 2013)<br />
`p.adjust(p, method = "BH")`:Benjamini & Hochberg (1995) ("BH" or its alias "fdr").<br />

#### 3.1 Compare differential exression across all clusters, generate CSV files
Generate CSV file in **_output_** folder.

#### 3.2~3.3 Zoom in Retinal pigment epithelium and Hematopoietic cells
Further, subset the cell types into small clusters.
Re-run `RunPCA()`, `FindClusters()`,`RunTSNE()`
Generate CSV file in **_output_** folder.

Below is a example of `./output/129_B6.csv` file with first 6 rows.

| row.name |   p_val | avg_logFC |  pct.1 |  pct.2 | p_val_adj |  pct.1_UMI |  pct.2_UMI |  cluster | gene
| -----    | ------  | -------- | ----  | ----- | ------- | ------- | ------| --- | --- |
| Lum    | 0.000 | 1.596 | 0.983 | 0.236 | 0.000 | 3.048 | 0.567 | 0)Pericytes | Lum
| Cygb    | 0.000 | 1.531 | 0.981 | 0.311 | 0.000 | 2.721 | 0.651 | 0)Pericytes | Cygb
| Igfbp4 | 0.000 | 1.507 | 1.000 | 0.720 | 0.000 | 3.897 | 1.773 | 0)Pericytes | Igfbp4
| Serpine2 | 0.000 | 1.462 | 0.999 | 0.506 | 0.000 | 3.488 | 1.191 | 0)Pericytes | Serpine2
| Dcn | 0.000 | 0.964 | 0.968 | 0.347 | 0.000 | 3.400 | 0.983 | 0)Pericytes | Dcn
| Cxcl12 | 0.000 | 1.222 | 0.981 | 0.531 | 0.000 | 2.957 | 1.148 | 0)Pericytes | Cxcl12


The results data frame has the following columns :

p_val: p_val (unadjusted) is calculated using likelihood-ratio test for single-cell gene expression, (McDavid et al., Bioinformatics, 2013) <br />
avg_logFC: log fold-change of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group.<br />
pct.1: The percentage of cells where the gene is detected in the first group.<br />
pct.2: The percentage of cells where the gene is detected in the second group.<br />
p_val_adj: Adjusted p-value, based on Benjamini & Hochberg (1995) ("BH" or its alias "fdr")<br />
pct.1_UMI is nUMI of the current cluster.
pct.2_UMI is average nUMI of rest of clusters.<br />
cluster : either cell types or original clusters in `./data/mouse_eyes_alignment.Rda`.
row.name and gene column are identical.<br />

##  
sessionInfo()
R version 3.4.4 (2018-03-15)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS High Sierra 10.13.4

Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] dplyr_0.7.4         Seurat_2.1.0        Biobase_2.38.0      BiocGenerics_0.24.0 Matrix_1.2-14       cowplot_0.9.2       ggplot2_2.2.1      

loaded via a namespace (and not attached):
  [1] diffusionMap_1.1-0   Rtsne_0.14           VGAM_1.0-5           colorspace_1.3-2    
  ggridges_0.5.0       class_7.3-14         modeltools_0.2-21    mclust_5.4          
  htmlTable_1.11.2     base64enc_0.1-3      proxy_0.4-22         rstudioapi_0.7      
  DRR_0.0.3            flexmix_2.3-14       lubridate_1.7.4      prodlim_2018.04.18  
  mvtnorm_1.0-7        ranger_0.9.0         codetools_0.2-15     splines_3.4.4       
  R.methodsS3_1.7.1    mnormt_1.5-5         doParallel_1.0.11    robustbase_0.93-0   
  knitr_1.20           tclust_1.3-1         RcppRoll_0.2.2       Formula_1.2-2       
  caret_6.0-79         ica_1.0-1            broom_0.4.4          gridBase_0.4-7      
  ddalpha_1.3.3        cluster_2.0.7-1      kernlab_0.9-26       R.oo_1.22.0         
  sfsmisc_1.1-2        compiler_3.4.4       backports_1.1.2      assertthat_0.2.0    
  lazyeval_0.2.1       lars_1.2             acepack_1.4.1        htmltools_0.3.6     
  tools_3.4.4          bindrcpp_0.2.2       igraph_1.1.0         gtable_0.2.0        
  glue_1.2.0           reshape2_1.4.3       Rcpp_0.12.16         NMF_0.21.0          
  trimcluster_0.1-2    gdata_2.18.0         ape_5.1              nlme_3.1-137        
  iterators_1.0.9      fpc_2.1-11           psych_1.8.3.3        timeDate_3043.102   
  gower_0.1.2          stringr_1.3.0        irlba_2.3.2          rngtools_1.2.4      
  gtools_3.5.0         DEoptimR_1.0-8       MASS_7.3-50          scales_0.5.0        
  ipred_0.9-6          RColorBrewer_1.1-2   yaml_2.1.19          pbapply_1.3-4       
  gridExtra_2.3        pkgmaker_0.22        segmented_0.5-3.0    rpart_4.1-13        
  latticeExtra_0.6-28  stringi_1.1.7        foreach_1.4.4        checkmate_1.8.5     
  caTools_1.17.1       ggjoy_0.4.0          lava_1.6.1           geometry_0.3-6      
  dtw_1.18-1           SDMTools_1.1-221     rlang_0.2.0          pkgconfig_2.0.1     
  prabclus_2.2-6       bitops_1.0-6         lattice_0.20-35      ROCR_1.0-7          
  purrr_0.2.4          bindr_0.1.1          recipes_0.1.2        htmlwidgets_1.2     
 tidyselect_0.2.4     CVST_0.2-1           plyr_1.8.4           magrittr_1.5        
 R6_2.2.2             gplots_3.0.1         Hmisc_4.1-1          dimRed_0.1.0        
 sn_1.5-2             withr_2.1.2          pillar_1.2.2         foreign_0.8-70      
 mixtools_1.1.0       survival_2.42-3      scatterplot3d_0.3-41 abind_1.4-5         
 nnet_7.3-12          tsne_0.1-3           tibble_1.4.2         KernSmooth_2.23-15  
 grid_3.4.4           data.table_1.10.4-3  FNN_1.1              ModelMetrics_1.1.0  
 digest_0.6.15        diptest_0.75-7       xtable_1.8-2         numDeriv_2016.8-1
 tidyr_0.8.0          R.utils_2.6.0        stats4_3.4.4         munsell_0.4.3    
 registry_0.5         magic_1.5-8