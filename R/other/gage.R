###############################################################################
### 1. Gene set and data preparation
###############################################################################

###################################################
### code chunk number 1-0: synopsis0 (eval = FALSE)
###################################################
## try http:// if https:// URLs are not supported
## source("https://bioconductor.org/biocLite.R")
## biocLite(c("gage","gageData","hgu133a.db","pathview"))

###################################################
### code chunk number 1-1: demo.data
###################################################
library(gage)
filename=system.file("extdata/gse16873.demo", package = "gage")
demo.data=readExpData(filename, row.names=1)
#check the data
head(demo.data)
str(demo.data)
#convert the data.frame into a matrix as to speed up the computing
demo.data=as.matrix(demo.data)
str(demo.data)


###################################################
### code chunk number 1-2: readList
###################################################
#an example GMT gene set data derived from MSigDB data
filename=system.file("extdata/c2.demo.gmt", package = "gage")
demo.gs=readList(filename)
demo.gs[1:3]
#to use these gene sets with gse16873, need to convert the gene symbols
#to Entrez IDs first
data(egSymb)
demo.gs.sym<-lapply(demo.gs, sym2eg)
demo.gs.sym[1:3]


###################################################
### code chunk number 1-3: gse16873.affyid
###################################################
library(gageData)
data(gse16873.affyid)
affyid=rownames(gse16873.affyid)

library(hgu133a.db)
egids2=hgu133aENTREZID[affyid]
annots=toTable(egids2)
str(annots)
gse16873.affyid=gse16873.affyid[annots$probe_id,]

#if multiple probe sets map to a gene, select the one with maximal IQR
iqrs=apply(gse16873.affyid, 1, IQR)
sel.rn=tapply(1:nrow(annots), annots$gene_id, function(x){
        x[which.max(iqrs[x])]
})
gse16873.egid=gse16873.affyid[sel.rn,]
rownames(gse16873.egid)=names(sel.rn)

cn=colnames(gse16873.egid)
hn=grep('HN',cn, ignore.case =T)
dcis=grep('DCIS',cn, ignore.case =T)
data(kegg.gs)
gse16873.kegg.p.affy <- gage(gse16873.egid, gsets = kegg.gs,
                             ref = hn, samp = dcis)
#result should be similar to that of using gse16873


###################################################
### code chunk number 1-4: pathview.conversion
###################################################
library(pathview)
data(bods)
print(bods)
#simulated human expression data with RefSeq ID
refseq.data <- sim.mol.data(mol.type = "gene", id.type = "REFSEQ",
                            nexp = 2, nmol = 1000)
#construct map between non-Entrez ID and Entrez Gene ID
id.map.refseq <- id2eg(ids = rownames(refseq.data), category =
                               "REFSEQ", org = "Hs")
#Map data onto Entrez Gene IDs, note different sum.method can be used
entrez.data <- mol.sum(mol.data = refseq.data, id.map = id.map.refseq,
                       sum.method = "mean")

###############################################################################
### 2. Generally Applicable Gene-set/Pathway Analysis
###############################################################################

###################################################
### code chunk number 2-0: synopsis0 (eval = FALSE)
###################################################
## ##step 0: setup (also need to map the reads outside R)
## source("http://bioconductor.org/biocLite.R")
## biocLite(c("pathview", "gage", "gageData", "GenomicAlignments",
##             "TxDb.Hsapiens.UCSC.hg19.knownGene"))


###################################################
### code chunk number 2: synopsis1 (eval = FALSE)
###################################################
## ##step 1: read counts
 library(TxDb.Hsapiens.UCSC.hg19.knownGene)
 exByGn <- exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene")
 library(GenomicAlignments)
 fls <- list.files("tophat_all/", pattern="bam$", full.names =T)
 bamfls <- BamFileList(fls)
 flag <- scanBamFlag(isSecondaryAlignment=FALSE, isProperPair=TRUE)
 param <- ScanBamParam(flag=flag)
 gnCnt <- summarizeOverlaps(exByGn, bamfls, mode="Union",
              ignore.strand=TRUE, singleEnd=FALSE, param=param)
 hnrnp.cnts=assay(gnCnt)


###################################################
### code chunk number 3: synopsis2 (eval = FALSE)
###################################################
## ##step 2: preprocessing
 require(gageData) #demo only
 data(hnrnp.cnts) #demo only
 cnts=hnrnp.cnts
 sel.rn=rowSums(cnts) != 0
 cnts=cnts[sel.rn,]
 ##joint workflow with DEseq/edgeR/limma/Cufflinks forks here
 libsizes=colSums(cnts)
 size.factor=libsizes/exp(mean(log(libsizes)))
 cnts.norm=t(t(cnts)/size.factor)
 cnts.norm=log2(cnts.norm+8)


###################################################
### code chunk number 4: synopsis3 (eval = FALSE)
###################################################
## ##step 3: gage
## ##joint workflow with DEseq/edgeR/limma/Cufflinks merges around here
 library(gage)
 ref.idx=5:8
 samp.idx=1:4
 data(kegg.gs)
 cnts.kegg.p <- gage(cnts.norm, gsets = kegg.gs, ref = ref.idx,
                     samp = samp.idx, compare ="unpaired")


###################################################
### code chunk number 5: synopsis4 (eval = FALSE)
###################################################
## ##step 4: pathview
 cnts.d= cnts.norm[, samp.idx]-rowMeans(cnts.norm[, ref.idx])
 sel <- cnts.kegg.p$greater[, "q.val"] < 0.1 &
          !is.na(cnts.kegg.p$greater[,"q.val"])
 path.ids <- rownames(cnts.kegg.p$greater)[sel]
 sel.l <- cnts.kegg.p$less[, "q.val"] < 0.1 &
            !is.na(cnts.kegg.p$less[,"q.val"])
 path.ids.l <- rownames(cnts.kegg.p$less)[sel.l]
 path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)
 library(pathview)
 pv.out.list <- sapply(path.ids2, function(pid) pathview(
                       gene.data = cnts.d, pathway.id = pid,
                       species = "hsa"))


###################################################
### code chunk number 6: start
###################################################
options(width=80)


###################################################
### code chunk number 7: install (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite(c("pathview", "gage", "gageData", "GenomicAlignments",
##             "TxDb.Hsapiens.UCSC.hg19.knownGene"))


###################################################
### code chunk number 8: readcount (eval = FALSE)
###################################################
## library(TxDb.Hsapiens.UCSC.hg19.knownGene)
## exByGn <- exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene")
## library(GenomicAlignments)
## fls <- list.files("tophat_all/", pattern="bam$", full.names =T)
## bamfls <- BamFileList(fls)
## flag <- scanBamFlag(isSecondaryAlignment=FALSE, isProperPair=TRUE)
## param <- ScanBamParam(flag=flag)
## #to run multiple core option: library(parallel); options("mc.cores"=4)
## gnCnt <- summarizeOverlaps(exByGn, bamfls, mode="Union",
##              ignore.strand=TRUE, singleEnd=FALSE, param=param)
## hnrnp.cnts=assay(gnCnt)


###################################################
### code chunk number 9: preprocessing
###################################################
require(gageData)
data(hnrnp.cnts)
cnts=hnrnp.cnts
dim(cnts)
sel.rn=rowSums(cnts) != 0
cnts=cnts[sel.rn,]
dim(cnts)
libsizes=colSums(cnts)
size.factor=libsizes/exp(mean(log(libsizes)))
cnts.norm=t(t(cnts)/size.factor)
range(cnts.norm)
cnts.norm=log2(cnts.norm+8)
range(cnts.norm)
#optional MA plot
pdf("hnrnp.cnts.maplots.pdf", width=8, height=10)
op=par(lwd=2, cex.axis=1.5, cex.lab=1.5, mfrow=c(2,1))
plot((cnts.norm[,6]+cnts.norm[,5])/2, (cnts.norm[,6]-cnts.norm[,5]), 
     main="(a) Control vs Control", xlab="mean", ylab="change",
     ylim=c(-5,5), xlim=c(0,20), lwd=1)
abline(h=0, lwd=2, col="red", lty="dashed")
plot((cnts.norm[,1]+cnts.norm[,5])/2, (cnts.norm[,1]-cnts.norm[,5]), 
     main="(b) Knockdown vs Control", xlab="mean", ylab="change",
     ylim=c(-5,5), xlim=c(0,20), lwd=1)
abline(h=0, lwd=2, col="red", lty="dashed")
dev.off()


###################################################
### code chunk number 10: gage
###################################################
library(gage)
ref.idx=5:8
samp.idx=1:4
data(kegg.gs)
#knockdown and control samples are unpaired
cnts.kegg.p <- gage(cnts.norm, gsets = kegg.gs, ref = ref.idx,
                    samp = samp.idx, compare ="unpaired")


###################################################
### code chunk number 11: pathview
###################################################
#differential expression: log2 ratio or fold change, uppaired samples
cnts.d= cnts.norm[, samp.idx]-rowMeans(cnts.norm[, ref.idx])

#up-regulated pathways (top 3) visualized by pathview
sel <- cnts.kegg.p$greater[, "q.val"] < 0.1 &
        !is.na(cnts.kegg.p$greater[,"q.val"])
path.ids <- rownames(cnts.kegg.p$greater)[sel]
path.ids2 <- substr(path.ids, 1, 8)
library(pathview)
pv.out.list <- sapply(path.ids2[1:3], function(pid) pathview(
        gene.data = cnts.d, pathway.id = pid,
        species = "hsa"))

#down-regulated pathways  (top 3) visualized by pathview
sel.l <- cnts.kegg.p$less[, "q.val"] < 0.1 &
        !is.na(cnts.kegg.p$less[,"q.val"])
path.ids.l <- rownames(cnts.kegg.p$less)[sel.l]
path.ids.l2 <- substr(path.ids.l, 1, 8)
pv.out.list.l <- sapply(path.ids.l2[1:3], function(pid) pathview(
        gene.data = cnts.d, pathway.id = pid,
        species = "hsa"))


###################################################
### code chunk number 12: goanalysis
###################################################
library(gageData)
data(go.sets.hs)
data(go.subs.hs)
lapply(go.subs.hs, head)
#Molecular Function analysis is quicker, hence run as demo
cnts.mf.p <- gage(cnts.norm, gsets = go.sets.hs[go.subs.hs$MF], 
                  ref = ref.idx, samp = samp.idx, compare ="unpaired")
#Biological Process analysis takes a few minutes if you try it
#cnts.bp.p <- gage(cnts.norm, gsets = go.sets.hs[go.subs.hs$BP], 
#    ref = ref.idx, samp = samp.idx, compare ="unpaired")


###################################################
### code chunk number 13: goresults
###################################################
for (gs in rownames(cnts.mf.p$less)[1:3]) {
        outname = gsub(" |:|/", "_", substr(gs, 12, 100))
        geneData(genes = go.sets.hs[[gs]], exprs = cnts.norm, ref = ref.idx, 
                 samp = samp.idx, outname = outname, txt = T, heatmap = T, 
                 limit = 3, scatterplot = T)
}


###################################################
### code chunk number 14: pergenescore
###################################################
cnts.t= apply(cnts.norm, 1, function(x) t.test(x[samp.idx], x[ref.idx],
                                               alternative = "two.sided", paired = F)$statistic)
cnts.meanfc= rowMeans(cnts.norm[, samp.idx]-cnts.norm[, ref.idx])
range(cnts.t)
range(cnts.meanfc)
cnts.t.kegg.p <- gage(cnts.t, gsets = kegg.gs, ref = NULL, samp = NULL)
cnts.meanfc.kegg.p <- gage(cnts.meanfc, gsets = kegg.gs, ref = NULL, samp = NULL)


###################################################
### code chunk number 15: deseq2
###################################################
library(DESeq2)
grp.idx <- rep(c("knockdown", "control"), each=4)
coldat=DataFrame(grp=factor(grp.idx))
dds <- DESeqDataSetFromMatrix(cnts, colData=coldat, design = ~ grp)
dds <- DESeq(dds)
deseq2.res <- results(dds)
#direction of fc, depends on levels(coldat$grp), the first level
#taken as reference (or control) and the second one as experiment.
deseq2.fc=deseq2.res$log2FoldChange
names(deseq2.fc)=rownames(deseq2.res)
exp.fc=deseq2.fc
out.suffix="deseq2"


###################################################
### code chunk number 16: deseq2
###################################################
require(gage)
data(kegg.gs)
fc.kegg.p <- gage(exp.fc, gsets = kegg.gs, ref = NULL, samp = NULL)
sel <- fc.kegg.p$greater[, "q.val"] < 0.1 &
        !is.na(fc.kegg.p$greater[, "q.val"])
path.ids <- rownames(fc.kegg.p$greater)[sel]
sel.l <- fc.kegg.p$less[, "q.val"] < 0.1 &
        !is.na(fc.kegg.p$less[,"q.val"])
path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)
require(pathview)
#view first 3 pathways as demo
pv.out.list <- sapply(path.ids2[1:3], function(pid) pathview(
        gene.data =  exp.fc, pathway.id = pid,
        species = "hsa", out.suffix=out.suffix))


###################################################
### code chunk number 17: deseq (eval = FALSE)
###################################################
## library(DESeq)
## grp.idx <- rep(c("knockdown", "control"), each=4)
## cds <- newCountDataSet(cnts, grp.idx)
## cds = estimateSizeFactors(cds)
## cds = estimateDispersions(cds)
## #this line takes several minutes
## system.time(
## deseq.res <- nbinomTest(cds, "knockdown", "control")
## )
## deseq.fc=deseq.res$log2FoldChange
## names(deseq.fc)=deseq.res$id
## sum(is.infinite(deseq.fc))
## deseq.fc[deseq.fc>10]=10
## deseq.fc[deseq.fc< -10]=-10
## exp.fc=deseq.fc
## out.suffix="deseq"


###################################################
### code chunk number 18: edger
###################################################
library(edgeR)
grp.idx <- rep(c("knockdown", "control"), each=4)
dgel <- DGEList(counts=cnts, group=factor(grp.idx))
dgel <- calcNormFactors(dgel)
dgel <- estimateCommonDisp(dgel)
dgel <- estimateTagwiseDisp(dgel)
et <- exactTest(dgel)
edger.fc=et$table$logFC
names(edger.fc)=rownames(et$table)
exp.fc=edger.fc
out.suffix="edger"


###################################################
### code chunk number 19: limma
###################################################
library(edgeR)
grp.idx <- rep(c("knockdown", "control"), each=4)
dgel2 <- DGEList(counts=cnts, group=factor(grp.idx))
dgel2 <- calcNormFactors(dgel2)
library(limma)
design <- model.matrix(~grp.idx)
log2.cpm <- voom(dgel2,design)
fit <- lmFit(log2.cpm,design)
fit <- eBayes(fit)
limma.res=topTable(fit,coef=2,n=Inf,sort="p")
limma.fc=limma.res$logFC
names(limma.fc)=limma.res$ID
exp.fc=limma.fc
out.suffix="limma"


###################################################
### code chunk number 20: cufflinks (eval = FALSE)
###################################################
## cuff.res=read.delim(file="gene_exp.diff", sep="\t")
## #notice the column name special character changes. The column used to be
## #cuff.res$log2.fold_change. for older versions of Cufflinks.
## cuff.fc=cuff.res$log2.FPKMy.FPKMx. 
## gnames=cuff.res$gene
## sel=gnames!="-"
## gnames=as.character(gnames[sel])
## cuff.fc=cuff.fc[sel]
## names(cuff.fc)=gnames
## gnames.eg=pathview::id2eg(gnames, category ="symbol")
## sel2=gnames.eg[,2]>""
## cuff.fc=cuff.fc[sel2]
## names(cuff.fc)=gnames.eg[sel2,2]
## range(cuff.fc)
## cuff.fc[cuff.fc>10]=10
## cuff.fc[cuff.fc< -10]=-10
## exp.fc=cuff.fc
## out.suffix="cuff"


