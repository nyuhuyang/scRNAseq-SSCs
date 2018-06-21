########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
source("./R/Seurat_functions.R")
########################################################################
#
#  1 Seurat Alignment 
# 
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# Load the dataset

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here
SSCs_raw <- list()
SSCs_Seurat <- list()
samples <- c("Ad_depleteSp","Ad_Thy1","PND14","PND18","PND18pre",
             "PND25","PND30","PND6" )
conditions <- c("Adault-SSCs","Adault-SSCs","first-wave","first-wave",
                "first-wave","first-wave","first-wave","first-wave" )

for(i in 1:length(samples)){
    SSCs_raw[[i]] <- Read10X(data.dir = paste0("./data/",
                                                     samples[i],"/outs/filtered_gene_bc_matrices/mm10/"))
    colnames(SSCs_raw[[i]]) <- paste0(samples[i],"-",colnames(SSCs_raw[[i]]))
    SSCs_Seurat[[i]] <- CreateSeuratObject(SSCs_raw[[i]],
                                           min.cells = 3,
                                           min.genes = 200,
                                           names.delim = "-")
    SSCs_Seurat[[i]]@meta.data$conditions <- conditions[i]
}
SSCs_Seurat <- lapply(SSCs_Seurat, FilterCells, 
                            subset.names = "nGene", 
                            low.thresholds = 200, 
                            high.thresholds = Inf)
SSCs_Seurat <- lapply(SSCs_Seurat, NormalizeData)
SSCs_Seurat <- lapply(SSCs_Seurat, ScaleData)
SSCs_Seurat <- lapply(SSCs_Seurat, FindVariableGenes, do.plot = FALSE)

# we will take the union of the top 1k variable genes in each dataset for
# alignment note that we use 1k genes in the manuscript examples, you can
# try this here with negligible changes to the overall results
genes.use <- lapply(SSCs_Seurat, function(x) head(rownames(x@hvg.info), 600))
genes.use <- unique(unlist(genes.use))
for(i in 1:length(samples)){
        genes.use <- intersect(genes.use, rownames(SSCs_Seurat[[i]]@scale.data))
}
length(genes.use) # 1/10 of total sample size


#======1.2 Perform a canonical correlation analysis (CCA) =========================
# run a canonical correlation analysis to identify common sources
# of variation between the two datasets.
remove(SSCs_raw)
GC()
SSCs <- RunMultiCCA(object.list = SSCs_Seurat, 
                    genes.use = genes.use,
                    niter = 25, num.ccs = 30,
                    standardize =TRUE)
save(SSCs, file = "./data/SSCs_alignment.Rda")

# CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = SSCs, reduction.use = "cca", group.by = "conditions", 
              pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = SSCs, features.plot = "CC1", group.by = "conditions", 
              do.return = TRUE)
plot_grid(p1, p2)

p3 <- MetageneBicorPlot(SSCs, grouping.var = "conditions", dims.eval = 1:30, 
                        display.progress = FALSE) # run on cluster
p3 + geom_smooth(method = 'loess')

PrintDim(object = SSCs, reduction.type = "cca", dims.print = 1:2, genes.print = 10)

DimHeatmap(object = SSCs, reduction.type = "cca", cells.use = 500, dim.use = c(1:3,11:20), 
           do.balanced = TRUE)

DimHeatmap(object = SSCs, reduction.type = "cca", cells.use = 500, dim.use = 10:18, 
           do.balanced = TRUE)

#======1.3 QC =========================
SSCs <- CalcVarExpRatio(object = SSCs, reduction.type = "pca",
                              grouping.var = "conditions", dims.use = 1:20)
SSCs <- SubsetData(SSCs, subset.name = "var.ratio.pca",accept.low = 0.5)

mito.genes <- grep(pattern = "^mt-", x = rownames(x = SSCs@data), value = TRUE)
percent.mito <- Matrix::colSums(SSCs@raw.data[mito.genes, ])/Matrix::colSums(SSCs@raw.data)
SSCs <- AddMetaData(object = SSCs, metadata = percent.mito, col.name = "percent.mito")
#SSCs <- ScaleData(object = SSCs, genes.use = genes.use, display.progress = FALSE, 
#                         vars.to.regress = "percent.mito")
#Now we can run a single integrated analysis on all cells!
VlnPlot(object = SSCs, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

SSCs <- FilterCells(object = SSCs, subset.names = c("nGene", "percent.mito"), 
                          low.thresholds = c(500, -Inf), high.thresholds = c(8000, 0.15))

par(mfrow = c(1, 2))
GenePlot(object = SSCs, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = SSCs, gene1 = "nUMI", gene2 = "nGene")

#======1.4 align seurat objects =========================
#Now we align the CCA subspaces, which returns a new dimensional reduction called cca.aligned
set.seed(42)
SSCs <- AlignSubspace(object = SSCs, reduction.type = "cca", grouping.var = "conditions", 
                            dims.align = 1:10)
#Now we can run a single integrated analysis on all cells!

SSCs <- FindClusters(object = SSCs, reduction.type = "cca.aligned", dims.use = 1:10, 
                           resolution = 0.6, force.recalc = T, save.SNN = TRUE)

SSCs <- RunTSNE(object = SSCs, reduction.use = "cca.aligned", dims.use = 1:10, 
                      do.fast = TRUE)

p1 <- TSNEPlot(SSCs, do.return = T, pt.size = 1, group.by = "conditions")
p2 <- TSNEPlot(SSCs, do.label = F, do.return = T, pt.size = 1)
#png('./output/TSNESplot_alignment.png')
plot_grid(p1, p2)
#dev.off()

TSNEPlot(object = SSCs,do.label = TRUE, group.by = "ident", 
         do.return = TRUE, no.legend = TRUE,
         pt.size = 1,label.size = 8 )+
    ggtitle("TSNEplot for all cell clusters")+
    theme(text = element_text(size=20),     #larger text including legend title							
          plot.title = element_text(hjust = 0.5)) #title in middle

SplitTSNEPlot(object = SSCs)
#dev.off()
save(SSCs, file = "./data/SSCs_alignment.Rda")
