# scRNAseq-SSCs
Dynamic transcriptome profiles within spermatogonial and spermatocyte populations during postnatal testis maturation revealed by single-cell sequencing
Kathryn J Grive, Yang Hu, Eileen Shu, Andrew Grimson, Olivier Elemento, Jennifer K Grenier, Paula E. Cohen
https://www.biorxiv.org/content/early/2018/11/06/464149

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
Tree structure of directory:
<pre>
|── LICENSE.md<br>
|── R<br>
|  |── Differential_analysis.R<br>
|  |── GSEA.R<br>
|  |── Identify_Cell_Types_Manually.R<br>
|  |── Seurat_functions.R<br>
|  |── Seurat_setup.R<br>
|  └── SingleR.R<br>
|── README.md<br>
|── data<br>
|  |── Ad-Thy1<br>
|  |  └── outs<br>
|  |      |── filtered_gene_bc_matrices<br>
|  |      |  └── mm10<br>
|  |      |      |── barcodes.tsv<br>
|  |      |      |── genes.tsv<br>
|  |      |      └── matrix.mtx<br>
...
|  |── PND30<br>
|  |  └── outs<br>
|  |      |── filtered_gene_bc_matrices<br>
|  |      |  └── mm10<br>
|  |      |      |── barcodes.tsv<br>
|  |      |      |── genes.tsv<br>
|  |      |      └── matrix.mtx<br>
|── doc<br>
|── output<br>
└── scRNAseq-SSCs.Rproj<br>
</pre>
</li>

### 1. Seurat_setup.R
<a href="https://github.com/nyuhuyang/scRNAseq-SSCs/blob/master/R/Seurat_setup.R">Seurat_setup.R</a><br />
The Cells with less than 500 genes or 2000 UMIs or more than 15% of mitochondria genes were excluded from the analysis. Gene expression raw counts were normalized following a global-scaling normalization method with a scale factor of 10,000 and a log 2 transformation, using the Seurat NormalizeData function. The top 4000 highly variable genes were selected using the expression and dispersion (variance/mean) of genes. Combat removed the batch effect. Seurat regressed the difference between the G2M and S phase, then followed by principal component analysis (PCA). The most significant principal components (1-30) were used for unsupervised clustering with 0.8 resolution and t-Distributed Stochastic Neighbor Embedding (tSNE). 

 After running this script, a `SSCs_(date).Rda` file will be generated inside **_data_** folder.
 Do not modify any files in **_data_** folder.
 
### 2. SingleR.R (optional)
<a href="https://github.com/nyuhuyang/scRNAseq-SSCs/blob/master/R/SingleR.R">SingleR.R</a><br />
Cell types were identified by SingleR (Single-cell Recognition) package. SingleR is a novel computational method for unbiased cell type recognition of scRNA-seq. SingleR leverages reference transcriptomic datasets of pure cell types to infer the cell of origin of each of the single cells independently.

We processed and annotated three reference mouse datasets:

Immunological Genome Project (ImmGen): a collection of 830 microarray samples, which we classified to 20 main cell types and further annotated to 253 subtypes.

A dataset of 358 mouse RNA-seq samples annotated to 28 cell types. This dataset was collected, processed and shared, courtesy of Bérénice Benayoun. This data set is especially useful for brain-related samples.

GSE43717: Cellular source and mechanisms of high transcriptome complexity in the mammalian testis (RNA-Seq cells). This dataset contains five mouse testis cell types including Sertoli cells, Spermatogonia, Spermatocytes, Spermatids, and Spermatozoa.

After running this script, a `singler_labels.RData` file will be generated inside **_output_** folder.

### 3. Identify_Cell_Types_Manually.R
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

`DotPlot` is implemented for visualizing differential expressed marker genes.
Multiple plots and table will be generated, save them when necessary. I prefer to keep the original identity of `SSCs_(date).Rda` intact for further downstream analysis.

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

### 5. GSEA.R
<a href="https://github.com/nyuhuyang/scRNAseq-SSCs/blob/master/R/GSEA.R">GSEA.R</a></li>
After running this script, a expression txt and label cls file will be generated inside **_output_** folder, for Gene Set Enrichment Analysis.