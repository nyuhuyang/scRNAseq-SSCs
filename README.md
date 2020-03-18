# Dynamic Transcriptome Profiles within Spermatogonial and Spermatocyte Populations During Postnatal Testis Maturation Revealed by Single-Cell Sequencing

The scripts developed for the single-cell RNA-seq study in [Grive, K. J., **Y. Hu**, E. Shu, A. Grimson, O. Elemento, J. K. Grenier, and P. E. Cohen. Dynamic Transcriptome Profiles within Spermatogonial and Spermatocyte Populations During Postnatal Testis Maturation Revealed by Single-Cell Sequencing. PLoS Genet (2019)](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007810)

## INTRODUCTION

## METHOD DETAILS

### Single-cell transcriptome analysis 
Data normalization, unsupervised cell clustering, and differential expression were carried out using the Seurat R package. Batch effect and cell-cycle effect were removed by Combat and Seurat together.

The Cells with less than 500 genes or 2000 UMIs or more than 15% of mitochondria genes were excluded from the analysis. Gene expression raw counts were normalized following a global-scaling normalization method with a scale factor of 10,000 and a log 2 transformation, using the Seurat NormalizeData function. The top 4000 highly variable genes were selected using the expression and dispersion (variance/mean) of genes. Combat removed the batch effect. Seurat regressed the difference between the G2M and S phase, then followed by principal component analysis (PCA). The most significant principal components (1-30) were used for unsupervised clustering and t-Distributed Stochastic Neighbor Embedding (tSNE). 

Cell types were manually identified by marker genes, and confirmed by SingleR (Single-cell Recognition) package. Differential expression analysis was performed based on the MAST (Model-based Analysis of Single Cell Transcriptomics). Gene Set Enrichment Analysis was performed based on GSEA software and Molecular Signature Database(MSigDB). Pathways were visualized by EnrichmentMap in Cytoscape.

## How to use this repository

#### Software Setup
R version 3.5.1<br />
dplyr_0.7.6 <br />
Seurat_2.3.4 <br />
SingleR_0.2.0 <br />

After pulling this repository, create folders **_data_** and **_output_** in the top working folder.
Move Cell Ranger analysis results into **_data_** folder.

### 1. Seurat_setup.R
<a href="https://github.com/nyuhuyang/scRNAseq-SSCs/blob/master/R/Seurat_setup.R">Seurat_setup.R</a><br />
The Cells with less than 500 genes or 2000 UMIs or more than 15% of mitochondria genes were excluded from the analysis. Gene expression raw counts were normalized following a global-scaling normalization method with a scale factor of 10,000 and a log 2 transformation, using the Seurat NormalizeData function. The top 4000 highly variable genes were selected using the expression and dispersion (variance/mean) of genes. Combat removed the batch effect. Seurat regressed the difference between the G2M and S phase, then followed by principal component analysis (PCA). The most significant principal components (1-30) were used for unsupervised clustering with 0.8 resolution and t-Distributed Stochastic Neighbor Embedding (tSNE). 

![](https://github.com/nyuhuyang/scRNAseq-SSCs/blob/master/Figs/TSNEPlot_PND06_PND14_PND18_PND18pre_PND25_PND30_Ad-depleteSp_Ad-Thy1_res.0.8.jpeg)
After running this script, a `SSCs_(date).Rda` file will be generated inside **_data_** folder.

### 2. Identify_Cell_Types_Manually.R
<a href="https://github.com/nyuhuyang/scRNAseq-SSCs/blob/master/R/Identify_Cell_Types_Manually.R">Identify_Cell_Types_Manually.R</a><br />
All clusters are tested against marker genes.
All cell types are predicted by at least two marker genes.

Smooth muscle were identified by Acta2 and Col1a2.<br />
Endothelial & Hematopoietic cells were identified by Laptm5, Ptprc, Vcam1, Insl3, Lyz2 and Hbb-bt.<br />
Sertoli cells were identified by Ptgds, Adgra3 and Wt1.<br />
Spermatogonia were identified by Gfra1, Zbtb16, Sall4, Dmrt1 and Dazl.<br />
Early Spermatocytes were identified by Kit and Dazl.<br />
Spermatocytes were identified by Id4, Sycp3, Mtl5, Nxt1 and Shcbp1l.<br />
Round Spermatids were identified by Lyzl1, Acrv1 and Hemgn.<br />
Spermatids were identified by Txndc8, Tssk6, Oaz3, Prm1, and Prm2.<br />

![](https://github.com/nyuhuyang/scRNAseq-SSCs/blob/master/Figs/Featureplot~.jpeg)
![](https://github.com/nyuhuyang/scRNAseq-SSCs/blob/master/Figs/dotplot.jpeg)

### 3. SingleR.R (optional)
<a href="https://github.com/nyuhuyang/scRNAseq-SSCs/blob/master/R/SingleR.R">SingleR.R</a><br />
Cell types were identified by SingleR (Single-cell Recognition) package. SingleR is a novel computational method for unbiased cell type recognition of scRNA-seq. SingleR leverages reference transcriptomic datasets of pure cell types to infer the cell of origin of each of the single cells independently.

We processed and annotated three reference mouse datasets:

Immunological Genome Project (ImmGen): a collection of 830 microarray samples, which we classified to 20 main cell types and further annotated to 253 subtypes.

A dataset of 358 mouse RNA-seq samples annotated to 28 cell types. This dataset was collected, processed and shared, courtesy of Bérénice Benayoun. This data set is especially useful for brain-related samples.

GSE43717: Cellular source and mechanisms of high transcriptome complexity in the mammalian testis (RNA-Seq cells). This dataset contains five mouse testis cell types including Sertoli cells, Spermatogonia, Spermatocytes, Spermatids, and Spermatozoa.

![](https://github.com/nyuhuyang/scRNAseq-SSCs/blob/master/Figs/cell_tpes.PNG)
![](https://github.com/nyuhuyang/scRNAseq-SSCs/blob/master/Figs/journal.pgen.1007810.g001.PNG)
After running this script, a `singler_labels.RData` file will be generated inside **_output_** folder.

### 4. Differential_analysis.R
<a href="https://github.com/nyuhuyang/scRNAseq-SSCs/blob/master/R/Differential_analysis.R">Differential_analysis.R</a><br />
Modified FindAllMarkers() `FindAllMarkers.UMI()` will generate similar dataframe plus two extra columns `UMI.1` and `UMI.2` to record nUMI. `UMI.1` is average nUMI of current cluster, `UMI.2` is average nUMI of all rest of clusters.<br />
`FindAllMarkers(object, test.use = "MAST")` : MAST (Model-based Analysis of Single Cell Transcriptomics), a GLM-framework that treates cellular detection rate as a covariate (Finak et al, Genome Biology, 2015)<br />

Below is an example of a Differential analysis output file.

| gene |   p_val | avg_logFC |  pct.1 |  pct.2 | p_val_adj |  UMI.1 |  UMI.2 |  cluster
| -----    | ------  | -------- | ----  | ----- | ------- | ------- | ------| --- |
| Acta2    | 0 | 3.9340 | 0.939 | 0.055 | 0 | 3.5565 | 0.0339 | Smooth muscle
| Myl9    | 0 | 2.9163 | 0.99 | 0.161 | 0 | 3.0622 | 0.1834 | Smooth muscle
| H19    | 0 | 2.9105 | 0.959 | 0.042 | 0 | 2.6070 | 0.0365 | Smooth muscle
| Meg3    | 0 | 2.7729 | 0.965 | 0.072 | 0 | 2.6652 | 0.0931 | Smooth muscle
| Col1a2 | 0 | 2.7221 | 0.995 | 0.035 | 0 | 2.7625 | 0.0496 | Smooth muscle

The results data frame has the following columns :

gene: gene name.<br />
p_val: p_val (unadjusted) is calculated using MAST (Model-based Analysis of Single Cell Transcriptomics, Finak et al., Genome Biology, 2015) <br />
avg_logFC: log fold-change of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group.<br />
pct.1: The percentage of cells where the gene is detected in the first group.<br />
pct.2: The percentage of cells where the gene is detected in the second group.<br />
p_val_adj: Adjusted p-value, based on bonferroni<br />
UMI.1 is average nUMI of the current cluster.<br />
UMI.2 is average nUMI of rest of clusters.<br />
cluster : either cell type or corresponding cluster.

![](https://github.com/nyuhuyang/scRNAseq-SSCs/blob/master/Figs/Doheatmap_top10_object_Cell.Types.jpeg)
### 5. GSEA.R
<a href="https://github.com/nyuhuyang/scRNAseq-SSCs/blob/master/R/GSEA.R">GSEA.R</a></li>
After running this script, a expression txt and label cls file will be generated inside **_output_** folder, for Gene Set Enrichment Analysis.

### 6. Main_figures.R
<a href="https://github.com/nyuhuyang/scRNAseq-SSCs/blob/master/R/Main_figures.R">GSEA.R</a></li>
Source code to produce all figures.
