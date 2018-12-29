########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
library(DropletUtils)
library(dplyr)
library(scater)
library(kableExtra)
#install_github("MarioniLab/scran") #BiocManager::install("BiocNeighbors", version = "devel")
library(scran)
library(EnsDb.Mmusculus.v79)
library(devtools)
library(Matrix)
library(devtools)
#library(scRNAseq)#BiocInstaller::biocLite("scRNAseq")
source("../R/Seurat_functions.R")
source("../R/scatter_utils.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, rLynchursive = T)
if(!dir.exists("./data/")) dir.create("data")
########################################################################
#
#  0.1~0.5 scater
# 
# ######################################################################
# 0.1. Setting up the data
# 0.1.1 Reading in a sparse matrix
df_samples <- readxl::read_excel("doc/20180530_sample_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
(samples <- df_samples$sample)
(conditions <- df_samples$conditions)

sce_list <- list()
for(i in 1:length(samples)){
        fname <- paste0("data/",samples[i],
                        "/outs/filtered_gene_bc_matrices/mm10")
        sce_list[[i]] <- read10xCounts.1(fname, col.names=TRUE,
                                         add.colnames = samples[i])
}

names(sce_list) <- samples
# 0.1.2 Annotating the rows
for(i in 1:length(samples)){
        rownames(sce_list[[i]]) <- uniquifyFeatureNames(rowData(sce_list[[i]])$ID,
                                                        rowData(sce_list[[i]])$Symbol)
        print(head(rownames(sce_list[[i]]),3))
        print(length(rownames(sce_list[[i]])))
}

# We also identify the chromosomal location for each gene. 
# The mitochondrial percentage is particularly useful for later quality control.
for(i in 1:length(samples)){
        location <- mapIds(EnsDb.Mmusculus.v79, keys=rownames(sce_list[[i]]), 
                           column="SEQNAME", keytype="GENENAME") #GENEID
        rowData(sce_list[[i]])$CHR <- location
        print(summary(location=="MT"))
}

# 0.2 Calling cells from empty droplets (skip)

########################################################################
# 0.2.2 Examining cell-calling diagnostics
# 0.4 Quality control on the cells#########################
# It is entirely possible for droplets to contain damaged or dying cells,
# which need to be removed prior to downstream analysis. 
# We compute some QC metrics using  calculateQCMetrics() (McCarthy et al. 2017) 
# and examine their distributions in Figure 2.
sce_list <- lapply(sce_list, function(x) calculateQCMetrics(x,compact = FALSE,
                        feature_controls=list(Mito=which(location=="MT"))))
########################################################################
## ---------------------------------------------
## ----qchist, Histograms of QC metric distributions in the dataset."----
for(i in 1:length(samples)){
        jpeg(paste0(path,"/0_QC_",samples[i],".jpeg"), units="in", width=10, height=7,
             res=600)
        par(mfrow=c(1,3))
                hist(sce_list[[i]]$log10_total_counts, breaks=20, col="grey80",
             xlab="Log-total UMI count",
             main = paste("nUMI count of",samples[i]))
        hist(sce_list[[i]]$log10_total_features_by_counts, breaks=20, col="grey80",
             xlab="Log-total number of expressed features",
             main = paste("nGene count of",samples[i]))
        hist(sce_list[[i]]$pct_counts_Mito, breaks=20, col="grey80",
             xlab="Proportion of reads in mitochondrial genes",
             main = paste("mitochondrial % of",samples[i]))
        dev.off()
}
## -------------------------------------------------
########################################################################

# Ideally, we would remove cells with low library sizes or total number of expressed features as described previously.
# However, this would likely remove cell types with low RNA content,
# especially in a heterogeneous population with many different cell types.
# Thus, we use a more relaxed strategy and only remove cells with large mitochondrial proportions,
# using it as a proxy for cell damage. 
# (Keep in mind that droplet-based datasets usually do not have spike-in RNA.)
# Low-quality cells are defined as those with extreme values for these QC metrics and are removed.
for(i in 1:length(samples)){
        high.mito <- isOutlier(sce_list[[i]]$pct_counts_Mito, nmads=3, type="higher")
        low.lib <- isOutlier(sce_list[[i]]$log10_total_counts, type="lower", nmad=3)
        low.genes <- isOutlier(sce_list[[i]]$log10_total_features_by_counts, type="lower", nmad=3)
        discard <- high.mito | low.lib | low.genes
        data.frame(HighMito= sum(high.mito),LowLib=sum(low.lib), 
                   LowNgenes=sum(low.genes),Discard=sum(discard))
        sce_list[[i]] <- sce_list[[i]][,!discard]
        print(summary(!discard))
}

########################################################################
#
#  0.6~ scran
# 
# ######################################################################
# Use natural Log transform to fit Seurat

for(i in 1:length(sce_list)){
        logcounts(sce_list[[i]]) <- as(log1p(assay(sce_list[[i]], "counts")),"dgCMatrix")
}
# 0.6 Normalizing for cell-specific biases
clusters <- list()
for(i in 1:length(sce_list)){
        clusters[[i]] <- quickCluster(sce_list[[i]], method="igraph", min.mean=0.1,
                                      assay.type = "logcounts",
                                 irlba.args=list(maxit=1000)) # for convergence.
        print(table(clusters[[i]]))
        sce_list[[i]] <- computeSumFactors(sce_list[[i]], min.mean=0.1, 
                                           cluster=clusters[[i]])
        print(summary(sizeFactors(sce_list[[i]])))
        
}

save(sce_list, file = "data/sce_8_20181225.Rda")


####################
#--------------------
## ----sfplot, fig.cap="Size factors for all cells in the PBMC dataset, plotted against the library size."----
jpeg(paste0(path,"/0_Size_factors2.jpeg"), units="in", width=10, height=7,
     res=600)
par(mfrow=c(3,2))
for(i in 7:length(samples)) plot(sce_list[[i]]$total_counts,
                                 sizeFactors(sce_list[[i]]),main=samples[i], log="xy")
dev.off()
#--------------------
####################
#Filter out low sizeFactors cells
sce_list[[9]] <- sce_list[[9]][,which(sizeFactors(sce_list[[9]]) >0.02)]
sce_list[[11]] <- sce_list[[11]][,which(sizeFactors(sce_list[[11]]) >0.02)]

sce_list <- lapply(sce_list, function(x) normalize(x,exprs_values = "logcounts"))

# 0.7 Modelling the mean-variance trend
# The lack of spike-in transcripts complicates the modelling of the technical noise.
# One option is to assume that most genes do not exhibit strong biological variation,
# and to fit a trend to the variances of endogenous genes. 
# However, this assumption is generally unreasonable for a heterogeneous population. 
# Instead, we assume that the technical noise is Poisson and create a fitted trend on that basis 
# using the makeTechTrend() function.
#for(i in 1:length(samples))  metadata(sce_list[[i]])$log.exprs.offset =1

new.trend <- lapply(sce_list, function(y) makeTechTrend(x=y))

## ----trendplot, "Variance of normalized log-expression values for each gene in the dataset,
# plotted against the mean log-expression. 
# The blue line represents the mean-dependent trend fitted to the variances, 
# while the red line represents the Poisson noise."----
fit <- lapply(sce_list, function(x) trendVar(x, use.spikes=FALSE, loess.args=list(span=0.05)))
######################################
## ---------------------------------
jpeg(paste0(path,"/0_Variance_values2.jpeg"), units="in", width=10, height=7,
     res=600)
par(mfrow=c(3,2))
for(i in 7:length(samples)) {
        plot(fit[[i]]$mean, fit[[i]]$var, pch=16, main =samples[i])
        curve(fit[[i]]$trend(x), col="dodgerblue", add=TRUE)
        curve(new.trend[[i]](x), col="red", add=TRUE)
}
dev.off()
## ---------------------------------
######################################
# decompose the variance for each gene using the Poisson-based trend, 
# and examine the genes with the highest biological components.

fit0 <- fit
for(i in 1:length(samples)) fit[[i]]$trend <- new.trend[[i]]

dec <- mapply(decomposeVar, x = sce_list, fit = fit)
top.dec <- list()
for(i in 1:length(samples)) {
        top.dec[[i]] <- dec[[i]][order(dec[[i]]$bio, decreasing=TRUE),] 
        print(length(row.names(top.dec[[i]])))
}

## ----hvgplot, "Distributions of normalized log-expression values for the top 10 genes 
# with the largest biological components in the dataset. 
# Each point represents the log-expression value in a single cell."----
plotExpression(sce_list[[1]], features=rownames(top.dec[[1]])[1:10])

########################################################################
#
#  0.8 batch correction
# 
# ######################################################################
# 0.8.1 Feature selection across batches
rownames_top.dec <- lapply(top.dec, function(x) rownames(x))
universe <- Reduce(intersect, rownames_top.dec)
bio_top.dec <- lapply(dec, function(x) x[,"bio"])
mean.bio <- Reduce("+", bio_top.dec)/length(samples)
chosen <- universe[mean.bio > 0]
length(chosen)

#sce_list_universe <- lapply(sce_list,function(x) x[universe,])
#rescaled <- Reduce(multiBatchNorm, sce_list_universe)
save(sce_list,chosen, file = "./data/sce_list.Rda")

# 0.9 Performing MNN-based correction
#https://bioconductor.org/packages/3.8/workflows/vignettes/simpleSingleCell/inst/doc/work-5-mnn.html#4_performing_mnn-based_correction

set.seed(100)
original <- lapply(sce_list, function(x) logcounts(x)[chosen,])
# Slightly convoluted call to avoid re-writing code later.
# Equivalent to fastMNN(GSE81076, GSE85241, k=20, d=50, approximate=TRUE)
mnn.out <- do.call(fastMNN, c(original, list(k=20, d=50, approximate=TRUE)))
dim(mnn.out$corrected)
mnn.out$batch
mnn.out$pair

# 0.10 Examining the effect of correction
# 0.10.1 prepare integrated SingleCellExperiment object
counts_list <- lapply(sce_list, function(x) counts(x))
logcounts_list <- lapply(sce_list, function(x) logcounts(x))
colData_list <- lapply(sce_list, function(x) colData(x))
rowData_list <- lapply(sce_list, function(x) rowData(x))

all.counts <- do.call(cbind, counts_list)
all.logcounts <- do.call(cbind, logcounts_list)
all.colData <- do.call(rbind, colData_list)

sce <- SingleCellExperiment(list(counts=all.counts,
                                 logcounts=all.logcounts),
                            colData=all.colData,
                            rowData=rowData_list[[1]][,1:3])
reducedDim(sce, "MNN") <- mnn.out$corrected
sce$Batch <- as.character(mnn.out$batch)
sce

set.seed(100)
# Using irlba to set up the t-SNE, for speed.
sce <- runPCA(sce, ntop=Inf, method="irlba",feature_set = chosen)
osce <- runTSNE(sce, use_dimred="PCA")
ot <- plotTSNE(osce, colour_by="Batch") + ggtitle("Original")

set.seed(100)
sce_MNN <- runTSNE(sce, use_dimred="MNN",feature_set = chosen)
ct <- plotTSNE(sce_MNN, colour_by="Batch") + ggtitle("Corrected")
save(sce_MNN, file = "./data/MNN_batch_correct.Rda")
######################################
## ---------------------------------
jpeg(paste0(path,"/0_batech_correction.jpeg"), units="in", width=10, height=7,
     res=600)
multiplot(ot, ct, cols=2)
dev.off()
## ---------------------------------
######################################
lnames = load(file = "./data/MNN_batch_correct.Rda");lnames
MCL <- as.seurat(sce_MNN)
MCL <- FindClusters(object = MCL, reduction.type = "MNN", dims.use = 1:20, 
                    resolution = 0.6, force.recalc = T, save.SNN = TRUE)
DimPlot(object = MCL, reduction.use = "TSNE", dim.1 = 1, dim.2 = 2, 
        group.by = "ident", do.return = TRUE)
save(MCL, file = "./data/MCL_MNN_20181017.Rda")
