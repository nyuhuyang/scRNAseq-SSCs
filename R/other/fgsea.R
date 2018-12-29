###################################################
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite(c("fgsea","org.Mm.eg.db","reactome.db"))

## ---generating exampleRanks---------------------------
library(GEOquery)
library(limma)
source("https://raw.githubusercontent.com/assaron/r-utils/master/R/exprs.R")
gse14308 <- GEOquery::getGEO("GSE14308")[[1]]
pData(gse14308)$condition <- sub("-.*$", "", gse14308$title)
es <- collapseBy(gse14308, fData(gse14308)$ENTREZ_GENE_ID, FUN=median)
es <- es[!grepl("///", rownames(es)), ]
es <- es[rownames(es) != "", ]
exprs(es) <- limma::normalizeBetweenArrays(log2(exprs(es)+1), method="quantile")
es.design <- model.matrix(~0+condition, data=pData(es))
fit <- lmFit(es, es.design)
fit2 <- contrasts.fit(fit, makeContrasts(conditionTh1-conditionNaive,
                                         levels=es.design))
fit2 <- eBayes(fit2)
de <- data.table(topTable(fit2, adjust.method="BH", number=12000, sort.by = "B"), keep.rownames = T)
ranks <- de[order(t), list(rn, t)]

# Visualise the six most significantly down- and up-regulated genes.
library(pheatmap)

my_group <- data.frame(group = pData(es)$condition)
row.names(my_group) <- colnames(exprs(es))

pheatmap(
        mat <- es[c(head(de[order(t), 1])$rn, tail(de[order(t), 1])$rn),],
        annotation_col = my_group,
        cluster_rows = FALSE,
        cellwidth=25,
        cellheight=15
)

## ----echo=F, message=F---------------------------------------------------
library(fgsea)
library(data.table)
library(ggplot2)

## ------------------------------------------------------------------------
data(examplePathways)
data(exampleRanks)
head(examplePathways)
head(exampleRanks)
## ------------------------------------------------------------------------
fgseaRes <- fgsea(pathways = examplePathways, 
                  stats = exampleRanks,
                  minSize=15,
                  maxSize=500,
                  nperm=10000)

## ------------------------------------------------------------------------
head(fgseaRes[order(pval), ])

## ------------------------------------------------------------------------
sum(fgseaRes[, padj < 0.01])

## ---- fig.width=7, fig.height=4------------------------------------------
plotEnrichment(examplePathways[["5991130_Programmed_Cell_Death"]],
               exampleRanks) + labs(title="Programmed Cell Death")

## ---- fig.width=7, fig.height=8, fig.retina=2----------------------------
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(examplePathways[topPathways], exampleRanks, fgseaRes, 
              gseaParam = 0.5)

## ---- fig.width=7, fig.height=8, fig.retina=2----------------------------
collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.01], 
                                      examplePathways, exampleRanks)
mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
        order(-NES), pathway]
plotGseaTable(examplePathways[mainPathways], exampleRanks, fgseaRes, 
              gseaParam = 0.5)

## ----message=FALSE-------------------------------------------------------
library(data.table)
fwrite(fgseaRes, file="./output/fgseaRes.txt", sep="\t", sep2=c("", " ", ""))

## ----message=FALSE-------------------------------------------------------
library(org.Mm.eg.db)
fgseaResMain <- fgseaRes[match(mainPathways, pathway)]
fgseaResMain[, leadingEdge := lapply(leadingEdge, mapIds, x=org.Mm.eg.db, keytype="ENTREZID", column="SYMBOL")]
fwrite(fgseaResMain, file="./output/fgseaResMain.txt", sep="\t", sep2=c("", " ", ""))

## ----message=F-----------------------------------------------------------
pathways <- reactomePathways(names(exampleRanks))
fgseaRes <- fgsea(pathways, exampleRanks, nperm=1000, maxSize=500)
head(fgseaRes)

## ------------------------------------------------------------------------
rnk.file <- system.file("extdata", "naive.vs.th1.rnk", package="fgsea")
gmt.file <- system.file("extdata", "mouse.reactome.gmt", package="fgsea")

## ------------------------------------------------------------------------
ranks <- read.table(rnk.file,
                    header=TRUE, colClasses = c("character", "numeric"))
ranks <- setNames(ranks$t, ranks$ID)
str(ranks)

## ------------------------------------------------------------------------
pathways <- gmtPathways(gmt.file)
str(head(pathways))

## ------------------------------------------------------------------------
fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize=500, nperm=1000)
head(fgseaRes)

#https://www.bioconductor.org/packages/devel/bioc/vignettes/EGSEA/inst/doc/EGSEA.pdf
