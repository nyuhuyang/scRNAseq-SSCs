########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(ggrepel)
source("../R/Seurat_functions.R")
########################################################################
#
#  1 Seurat Alignment 
# 
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# Load the dataset

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here

# rename all "_" into "-" in sample names 
SSCs_raw <- list()
SSCs_Seurat <- list()
samples <- c("Ad-depleteSp","Ad-Thy1","PND14","PND18","PND18pre",
             "PND25","PND30","PND6" )
conditions <- c("Adault-SSCs","Adault-SSCs","first-wave","first-wave",
                "first-wave","first-wave","first-wave","first-wave" )

for(i in 1:length(samples)){
    SSCs_raw[[i]] <- Read10X(data.dir = paste0("./data/",
                                                     samples[i],"/outs/filtered_gene_bc_matrices/mm10/"))
    colnames(SSCs_raw[[i]]) <- paste0(samples[i],"_",colnames(SSCs_raw[[i]]))
    SSCs_Seurat[[i]] <- CreateSeuratObject(SSCs_raw[[i]],
                                           min.cells = 3,
                                           min.genes = 200,
                                           names.delim = "_")
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
                    niter = 25, num.ccs = 50,
                    standardize =TRUE)
save(SSCs, file = "./data/SSCs_alignment.Rda")

# Calculate median UMI per cell
SSCs_raw_data <- as.matrix(x = SSCs@raw.data)
mean(colSums(SSCs_raw_data))
median(colSums(SSCs_raw_data))
boxplot(colSums(SSCs_raw_data))
# CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = SSCs, reduction.use = "cca", group.by = "orig.ident", 
              pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = SSCs, features.plot = "CC1", group.by = "orig.ident", 
              do.return = TRUE)
plot_grid(p1, p2)

AverageExpression()
p3 <- MetageneBicorPlot(SSCs, grouping.var = "orig.ident", dims.eval = 1:50, 
                        display.progress = TRUE, display.progress = FALSE,
                        smooth = TRUE) # run on cluster
p3 + geom_smooth(method = 'loess')

PrintDim(object = SSCs, reduction.type = "cca", dims.print = 1:2, genes.print = 10)

DimHeatmap(object = SSCs, reduction.type = "cca", cells.use = 500, dim.use = c(1:9), 
           do.balanced = TRUE)
DimHeatmap(object = SSCs, reduction.type = "cca", cells.use = 500, dim.use = c(10:18), 
           do.balanced = TRUE)
#======1.3 QC =========================
SSCs <- CalcVarExpRatio(object = SSCs, reduction.type = "pca",
                              grouping.var = "orig.ident", dims.use = 1:10)
SSCs <- SubsetData(SSCs, subset.name = "var.ratio.pca",accept.low = 0.5) #15555 out of 15650

mito.genes <- grep(pattern = "^mt-", x = rownames(x = SSCs@data), value = TRUE)
percent.mito <- Matrix::colSums(SSCs@raw.data[mito.genes, ])/Matrix::colSums(SSCs@raw.data)
SSCs <- AddMetaData(object = SSCs, metadata = percent.mito, col.name = "percent.mito")
#SSCs <- ScaleData(object = SSCs, genes.use = genes.use, display.progress = FALSE, 
#                         vars.to.regress = "percent.mito")
#Now we can run a single integrated analysis on all cells!
VlnPlot(object = SSCs, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 1)
SSCs <- FilterCells(object = SSCs, subset.names = c("nGene", "percent.mito"), 
                          low.thresholds = c(500, -Inf), high.thresholds = c(7000, 0.15))

par(mfrow = c(1, 2))
GenePlot(object = SSCs, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = SSCs, gene1 = "nUMI", gene2 = "nGene")

#======1.4 align seurat objects =========================
#Now we align the CCA subspaces, which returns a new dimensional reduction called cca.aligned
set.seed(42)
SSCs <- AlignSubspace(object = SSCs, reduction.type = "cca", grouping.var = "orig.ident", 
                            dims.align = 1:10)
SSCs@project.name <- "Paula"
#Now we can run a single integrated analysis on all cells!

SSCs <- FindClusters(object = SSCs, reduction.type = "cca.aligned", dims.use = 1:10, 
                           resolution = 1.0, force.recalc = T, save.SNN = TRUE)

SSCs <- RunTSNE(object = SSCs, reduction.use = "cca.aligned", dims.use = 1:10, 
                      do.fast = TRUE)

p1 <- TSNEPlot(SSCs, do.return = T, pt.size = 1, group.by = "orig.ident")
p2 <- TSNEPlot(SSCs, do.return = T, pt.size = 1, group.by = "ident")
#png('./output/TSNESplot_alignment.png')
plot_grid(p1, p2)
#dev.off()
TSNEPlot(object = SSCs,do.label = T, group.by = "ident", 
         do.return = TRUE, no.legend = T,
         pt.size = 1,label.size = 8 )+
        ggtitle("All samples")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle

#======1.5 RunPCA==================================
SSCs <- FindVariableGenes(object = SSCs, mean.function = ExpMean, dispersion.function = LogVMR, 
                         do.plot = FALSE)
hv.genes <- head(rownames(SSCs@hvg.info), 1000)
SSCs <- RunPCA(object = SSCs, pc.genes = hv.genes, pcs.compute = 100, do.print = TRUE, 
              pcs.print = 1:5, genes.print = 5)
PCElbowPlot(object = SSCs, num.pc = 100)
PCHeatmap(SSCs, pc.use = c(1:3, 25:30), cells.use = 500, do.balanced = TRUE)
SSCs <- FindClusters(object = SSCs, reduction.type = "pca", dims.use = 1:30, resolution = 1, 
                    save.SNN = TRUE, n.start = 10, nn.eps = 0.5, print.output = FALSE)
SSCs <- RunTSNE(object = SSCs, reduction.use = "pca", dims.use = 1:30, 
                do.fast = TRUE)
TSNEPlot(object = SSCs,do.label = T, group.by = "orig.ident", 
         do.return = TRUE, no.legend = T,
         pt.size = 1,label.size = 8 )+
        ggtitle("All samples")+
        theme(text = element_text(size=15),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) 

save(SSCs, file = "./data/SSCs_alignment.Rda")
