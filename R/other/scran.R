## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("scran")
library(scran)
browseVignettes("scran")
## ---- echo=FALSE, results="hide", message=FALSE--------------------------
require(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)

## ----style, echo=FALSE, results='asis'-----------------------------------
BiocStyle::markdown()

## ----setup, echo=FALSE, message=FALSE------------------------------------
set.seed(100)

## ------------------------------------------------------------------------
ngenes <- 10000
ncells <- 200
mu <- 2^runif(ngenes, -1, 5)
gene.counts <- matrix(rnbinom(ngenes*ncells, mu=mu, size=10), nrow=ngenes)

## ------------------------------------------------------------------------
library(org.Mm.eg.db)
all.ensembl <- unique(toTable(org.Mm.egENSEMBL)$ensembl_id)
rownames(gene.counts) <- sample(all.ensembl, ngenes)

## ------------------------------------------------------------------------
nspikes <- 100
ncells <- 200
mu <- 2^runif(nspikes, -1, 5)
spike.counts <- matrix(rnbinom(nspikes*ncells, mu=mu, size=10), nrow=nspikes)
rownames(spike.counts) <- paste0("ERCC-", seq_len(nspikes))
all.counts <- rbind(gene.counts, spike.counts)

## ------------------------------------------------------------------------
library(scran)
sce <- SingleCellExperiment(list(counts=all.counts))
isSpike(sce, "MySpike") <- grep("^ERCC", rownames(sce))

## ------------------------------------------------------------------------
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))

## ------------------------------------------------------------------------
assigned <- cyclone(sce, pairs=mm.pairs)
head(assigned$scores)

## ------------------------------------------------------------------------
table(assigned$phases)

## ------------------------------------------------------------------------
sce <- computeSumFactors(sce)
summary(sizeFactors(sce))

## ------------------------------------------------------------------------
larger.sce <- SingleCellExperiment(list(counts=cbind(all.counts, all.counts, all.counts)))
clusters <- quickCluster(larger.sce)
larger.sce <- computeSumFactors(larger.sce, cluster=clusters)

## ------------------------------------------------------------------------
sce2 <- computeSpikeFactors(sce)
summary(sizeFactors(sce2))

## ------------------------------------------------------------------------
sce <- computeSpikeFactors(sce, general.use=FALSE)

## ------------------------------------------------------------------------
sce <- normalize(sce)

## ------------------------------------------------------------------------
fit <- trendVar(sce, parametric=TRUE)

## ------------------------------------------------------------------------
decomp <- decomposeVar(sce, fit)
top.hvgs <- order(decomp$bio, decreasing=TRUE)
head(decomp[top.hvgs,])

## ---- fig.cap=""---------------------------------------------------------
plot(decomp$mean, decomp$total, xlab="Mean log-expression", ylab="Variance")
o <- order(decomp$mean)
lines(decomp$mean[o], decomp$tech[o], col="red", lwd=2)
points(fit$mean, fit$var, col="red", pch=16)

## ------------------------------------------------------------------------
alt.fit <- trendVar(sce, use.spikes=FALSE) 
alt.decomp <- decomposeVar(sce, alt.fit)

## ------------------------------------------------------------------------
batch <- rep(c("1", "2"), each=100)
design <- model.matrix(~batch)
alt.fit2 <- trendVar(sce, design=design)
alt.decomp2 <- decomposeVar(sce, alt.fit)

## ------------------------------------------------------------------------
null.dist <- correlateNull(ncol(sce))
# Only using the first 200 genes as a demonstration.
cor.pairs <- correlatePairs(sce, subset.row=top.hvgs[1:200], null.dist=null.dist)
head(cor.pairs)

## ------------------------------------------------------------------------
null.dist2 <- correlateNull(design=design, iter=1e5) # fewer iterations, to speed it up.
cor.pairs2 <- correlatePairs(sce, subset.row=top.hvgs[1:200], 
                             null.dist=null.dist2, design=design)

## ------------------------------------------------------------------------
cor.genes <- correlatePairs(sce, subset.row=top.hvgs[1:200], 
                            null.dist=null.dist, per.gene=TRUE)

## ------------------------------------------------------------------------
y <- convertTo(sce, type="edgeR")

## ------------------------------------------------------------------------
sessionInfo()

