########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(Matrix)
library(sva)
library(SingleR)
source("../R/Seurat_functions.R")
########################################################################
#
#  1 Data preprocessing
# 
# ######################################################################
#======1.1 Load the data files and Set up Seurat object =========================
# Load the dataset

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here

# rename all "_" into "-" in sample names 
SSCs_raw <- list()
SSCs_Seurat <- list()
samples <- c("PND06","PND14","PND18","PND18pre",
             "PND25","PND30","Ad-depleteSp","Ad-Thy1")
conditions <- c("first-wave","first-wave","first-wave","first-wave",
                "first-wave","first-wave","Adault-SSCs","Adault-SSCs")
for(i in 1:length(samples)){
        SSCs_raw[[i]] <- Read10X(data.dir = paste0("./data/",
                                                   samples[i],"/outs/filtered_gene_bc_matrices/mm10/"))
        colnames(SSCs_raw[[i]]) <- paste0(samples[i],"_",colnames(SSCs_raw[[i]]))
        SSCs_Seurat[[i]] <- CreateSeuratObject(SSCs_raw[[i]],
                                               min.cells = 3,
                                               min.genes = 200,
                                               names.delim = "_",
                                               project = "paula")
        SSCs_Seurat[[i]]@meta.data$conditions <- conditions[i]
}
SSCs <- Reduce(function(x, y) MergeSeurat(x, y, do.normalize = F), SSCs_Seurat)
remove(SSCs_raw,SSCs_Seurat);GC()
SSCs <- FilterCells(SSCs, subset.names = "nGene",
                    low.thresholds = 200,
                    high.thresholds = Inf) %>%
        NormalizeData() %>%
        ScaleData(display.progress = FALSE) %>%
        FindVariableGenes(do.plot = FALSE, display.progress = FALSE)
save(SSCs, file = "./data/SSCs_20180825.Rda")

#======1.2 QC, pre-processing and normalizing the data=========================
# 1.2.1 Calculate median UMI per cell
SSCs_raw_data <- as.matrix(x = SSCs@raw.data)
mean(colSums(SSCs_raw_data))
median(colSums(SSCs_raw_data))
min(colSums(SSCs_raw_data))
remove(SSCs_raw_data);GC()

# 1.2.2 Identifying Outlier Cells
# Plot genes per cell
# How many genes expressed per cells
complexity.per.cell <- apply(SSCs@data, 2, function(x) sum(x>0))
# Mean count per cell.
mean.count.per.cell <- apply(SSCs@data, 2, function(x) mean(x))
# Gene prevalence
gene.prevalence <- apply(SSCs@raw.data, 1, function(x) sum(x>0))
# Complexity by mean expression
par(mfrow = c(2,1))
plot(complexity.per.cell, mean.count.per.cell)

# 1.2.3 calculate mitochondria percentage
mito.genes <- grep(pattern = "^mt-", x = rownames(x = SSCs@data), value = TRUE)
percent.mito <- Matrix::colSums(SSCs@raw.data[mito.genes, ])/Matrix::colSums(SSCs@raw.data)
SSCs <- AddMetaData(object = SSCs, metadata = percent.mito, col.name = "percent.mito")

g1 <- VlnPlot(object = SSCs, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 1,
             x.lab.rot = T, do.return = T)
g1
SSCs <- FilterCells(object = SSCs, subset.names = c("nGene","nUMI","percent.mito"),
                    low.thresholds = c(500,2000, -Inf), 
                    high.thresholds = c(8000,125000, 0.15))

g2 <- VlnPlot(object = SSCs, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 1,
              x.lab.rot = T, do.return = T)

par(mfrow = c(2, 1))
GenePlot(object = SSCs, gene1 = "nUMI", gene2 = "percent.mito",use.raw = T)
GenePlot(object = SSCs, gene1 = "nUMI", gene2 = "nUMI",use.raw = T)

VlnPlot(object = SSCs, features.plot = c("nUMI"), nCol = 1,
              group.by = "ident", x.lab.rot = T, do.return = T)
plot_grid(g1,g2)
# After removing unwanted cells from the dataset, the next step is to normalize the data.
SSCs <- NormalizeData(object = SSCs, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
SSCs <- FindVariableGenes(object = SSCs, mean.function = ExpMean, 
                          dispersion.function = LogVMR, do.plot = FALSE, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(SSCs@var.genes)
#======1.3 1st run of pca-tsne  =========================
SSCs <- ScaleData(object = SSCs) %>%
        RunPCA() %>%
        FindClusters(dims.use = 1:20, force.recalc = T, print.output = FALSE) %>%
        RunTSNE()
#SSCs@meta.data$orig.ident <- gsub("PND18pre","PND18",SSCs@meta.data$orig.ident)
TSNEPlot(object = SSCs, do.label = F, group.by = "orig.ident", 
         do.return = TRUE, no.legend = F,# colors.use = singler.colors,
         pt.size = 1,label.size = 8 )+
        ggtitle("TSNEPlot of all samples")+
        theme(text = element_text(size=15),							
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

save(SSCs, file = "./data/SSCs_20180825.Rda")
Iname = load("./data/SSCs_20180825.Rda")

#======1.4 Add Cell-cycle score =========================
# Read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "./data/seurat_resources/regev_lab_cell_cycle_genes.txt")
# We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- MouseGenes(SSCs,cc.genes[1:43])
g2m.genes <- MouseGenes(SSCs,cc.genes[44:97])
# Assign Cell-Cycle Scores
SSCs <- CellCycleScoring(object = SSCs, s.genes = s.genes, g2m.genes = g2m.genes, 
                         set.ident = TRUE)
# Visualize the distribution of cell cycle markers across
RidgePlot(object = SSCs, features.plot = MouseGenes(SSCs,c("PCNA", "TOP2A", "MCM6", "MKI67")), 
          nCol = 2)
# regressing out the difference between the G2M and S phase scores
SSCs@meta.data$CC.Difference <- SSCs@meta.data$S.Score - SSCs@meta.data$G2M.Score
# view cell cycle scores and phase assignments
head(x = SSCs@meta.data)

#======1.5 Add batch id =========================
batchname = SSCs@meta.data$orig.ident
batch.effect = rep(NA,length(batchname))
batch.effect[batchname %in% c("PND14","PND18","PND18pre","PND30")] = 1
batch.effect[batchname %in% c("PND06","PND25","Ad-depleteSp","Ad-Thy1")] = 2
names(batch.effect) = rownames(SSCs@meta.data)
SSCs <- AddMetaData(object = SSCs, metadata = batch.effect, col.name = "batch.effect")
table(SSCs@meta.data$batch.effect)
head(x = SSCs@meta.data)
#======1.6 batch-correct using ComBat =========================
SingleFeaturePlot.1(SSCs,"nUMI",threshold=15000)
SingleFeaturePlot.1(SSCs,"batch.effect",threshold=1.0)
SingleFeaturePlot.1(SSCs,"percent.mito",threshold=0.05)
SingleFeaturePlot.1(SSCs,"CC.Difference",threshold=0.05)
m = as.matrix(SSCs@data)
m = m[rowSums(m)>0,]
com = ComBat(m, batch.effect, prior.plots=FALSE, par.prior=TRUE)
#----save files just in case------
saveRDS(Matrix(as.matrix(com)), file = "./data/Combat_data.Rda")
saveRDS(SSCs@data, file = "./data/SSCs_data.Rda")
remove(m,com);GC()
#---------------------
SSCs@data = readRDS("./data/Combat_data.Rda")
SSCs@scale.data = NULL
SSCs <- ScaleData(object = SSCs,#genes.use = SSCs@var.genes,
                  model.use = "negbinom", do.par=T, do.center = T, do.scale = T,
                  vars.to.regress = c("CC.Difference"),#"CC.Difference","percent.mito"--nogood,"nUMI"--nogood
                  display.progress = T)
SSCs@data =  readRDS("./data/SSCs_data.Rda")

gene.use <- rownames(SSCs@scale.data);length(gene.use)
SSCs@data = SSCs@data[gene.use,]
#save(SSCs, file = "./data/SSCs_20181001.Rda") #do.center = F, do.scale = T
#======1.7 unsupervised clustering =========================
SSCs <- RunPCA(object = SSCs, pc.genes = SSCs@var.genes, pcs.compute = 100, 
               do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCAPlot(object = SSCs)
PCElbowPlot(object = SSCs, num.pc = 100)
PCHeatmap(SSCs, pc.use = c(1:3, 25:30), cells.use = 500, do.balanced = TRUE)
SSCs <- RunTSNE(object = SSCs, reduction.use = "pca", dims.use = 1:30, 
                do.fast = TRUE, perplexity= 30)

SSCs <- FindClusters(object = SSCs, reduction.type = "pca", dims.use = 1:30, 
                     resolution = 0.8, 
                     k.param = 30,force.recalc = T,
                     save.SNN = TRUE, n.start = 100, nn.eps = 0, print.output = FALSE)
#SSCs@meta.data$orig.ident <- gsub("PND18pre","PND18",SSCs@meta.data$orig.ident)
TSNEPlot.1(object = SSCs, do.label = T, group.by = "ident", 
         do.return = TRUE, no.legend = T, 
         text.repel = T, label.repel = F,
        #colors.use = singler.colors,
         pt.size = 1,label.size = 6 )+
        ggtitle("TSNEPlot, resolution = 0.8, k.param = 30")+
        theme(text = element_text(size=15),							
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 
SSCs <- StashIdent(object = SSCs, save.name = "res.0.8")
save(SSCs, file = "./data/SSCs_20180825.Rda")
