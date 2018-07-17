library(GSEABase)
library(hgu95av2.db)
library(GO.db)

data(sample.ExpressionSet) # from Biobase
experimentData(sample.ExpressionSet)
exprs(sample.ExpressionSet)[1:5,1:5]
head(pData(sample.ExpressionSet))
phenoData(sample.ExpressionSet)
head(featureData(sample.ExpressionSet))
dim(sample.ExpressionSet)
geneIds(sample.ExpressionSet)

## GeneSet from Broad XML; 'fl' could be a url
fl <- system.file("extdata", "Broad.xml", package="GSEABase")
gs2 <- getBroadSets(fl)[[2]] # actually, a list of two gene sets
geneIds <- geneIds(gs2) # any character vector would do
gs3 <- GeneSet(geneIds, setName="geneIds")
###################################################
### code chunk number 3: GeneSet
###################################################

egs <- GeneSet(sample.ExpressionSet[201:250,], setName="Sample")
egs


###################################################
### code chunk number 4: geneIds
###################################################
head(geneIds(egs))


###################################################
### code chunk number 5: details
###################################################
details(egs)


###################################################
### code chunk number 6: GSEABase.Rnw:96-99
###################################################
## FIXME: GeneSet(AnnotationIdentifier("hgu95av2")) --> non-empty
## FIXME: GeneSet(AnnotationIdentifier("hgu95av2"),
## collectionType=GOCollection()) filters on GOCollection (or KEGG)


###################################################
### code chunk number 7: GeneSet-methods
###################################################
showMethods("GeneSet", inherited=FALSE)


###################################################
### code chunk number 8: GeneIdentifierTypes
###################################################
names(slot(getClass("GeneIdentifierType"), "subclasses"))


###################################################
### code chunk number 9: mapIdentifiers
###################################################
mapIdentifiers(egs, EntrezIdentifier())


###################################################
### code chunk number 10: GeneSet_Identifiers
###################################################
library(annotate)                       # getEG
eids <- unique(getEG(geneIds(egs), "hgu95av2"))
eids <- eids[!is.na(eids)]
GeneSet(EntrezIdentifier(), geneIds=as.character(eids))


###################################################
### code chunk number 11: CollectionType
###################################################
names(slot(getClass("CollectionType"), "subclasses"))


###################################################
### code chunk number 12: GOCollection
###################################################
GeneSet(GOCollection(c("GO:0005488", "GO:0019825"),
                     evidenceCode="IDA"),
        geneIdType=EntrezIdentifier("org.Hs.eg.db"),
        setName="Sample GO Collection")


###################################################
### code chunk number 13: Broad
###################################################
fl <- system.file("extdata", "Broad1.xml", package="GSEABase")
bgs <- GeneSet(BroadCollection(), urls=fl)
bgs


###################################################
### code chunk number 14: Broad-to-annotation
###################################################
bgs1 <- mapIdentifiers(bgs, AnnotationIdentifier("hgu95av2"))
bgs1


###################################################
### code chunk number 15: subset
###################################################
bgs[1:5]
bgs[c("GALNS", "LOC646365")]


###################################################
### code chunk number 16: egs-bgs
###################################################
egs & bgs1


###################################################
### code chunk number 17: subset-ExpressionSet
###################################################
sample.ExpressionSet[bgs,]


###################################################
### code chunk number 18: GeneColorSet-setup
###################################################
conn <- textConnection("
                       Entrez ID, Gene Symbol, Expression level, Phenotype response
                       ##used to be MRP2
                       1244, ABCC2, Increase, Resistant
                       538, ATP7A, Increase, Resistant
                       540, ATP7B, Increase, Resistant
                       9961, MVP, Increase, Resistant
                       ##the LRP below must be MVP
                       ##LRP, Increase, Resistant - need to know which one
                       7507,XPA, Increase, Resistant
                       2067, ERCC1, Increase, Resistant
                       ##TOP, Increase, Resistant  - need to know which one, notes say II
                       672, BRCA1, Increase, Resistant
                       3725, JUN, Increase, Resistant
                       #GCS, Increase, Resistant  - my notes say alpha-GCS - so which one?
                       ##I only found gamma at PubMed as being related
                       2730, GCLM, Increase, Resistant")
tbl <- read.csv(conn, strip.white=TRUE,comment.char="#")
close(conn); unlink(conn)


###################################################
### code chunk number 19: GeneColorSet-phenotype
###################################################
tbl


###################################################
### code chunk number 20: GeneColorSet-constructor
###################################################
gcs <- GeneColorSet(EntrezIdentifier(),
                    setName="A color set",
                    geneIds=as.character(tbl$Entrez.ID),
                    phenotype="Cisplatin resistance",
                    geneColor=tbl$Expression.level,
                    phenotypeColor=tbl$Phenotype.response)
gcs


###################################################
### code chunk number 21: GeneSetCollection
###################################################
gsc <- GeneSetCollection(sample.ExpressionSet[201:250,], setType=GOCollection())
gsc
gsc[["GO:0005737"]]


###################################################
### code chunk number 22: GeneSetCollection-GOCollection
###################################################
GeneSetCollection(sample.ExpressionSet[201:300,],
                  setType=GOCollection(evidenceCode="IMP"))


###################################################
### code chunk number 23: GeneSetCollection-BroadCollection
###################################################
## FIXME: BroadCollection default to paste("c", 1:4, sep="")
## FIXME: GeneSetCollection(BroadCollection(), urls=fl); filters on bcCategory
fl <- system.file("extdata", "Broad.xml", package="GSEABase")
gss <- getBroadSets(fl)
gss
names(gss)


###################################################
### code chunk number 24: mapIds-GeneSetCollection
###################################################
gsc <- mapIdentifiers(gsc, EntrezIdentifier())
gsc
gsc[["GO:0005737"]]


###################################################
### code chunk number 25: ReportingTools
###################################################
## 'interesting' gene sets
idx <- sapply(gsc, function(x) length(geneIds(x))) > 2

library(ReportingTools)
gscReport <- HTMLReport(
        shortName="gsc_example",
        title="GSEABase Vignette GeneSetCollection",
        basePath=tempdir())
publish(gsc[idx], gscReport, annotation.db="org.Hs.eg")
url <- finish(gscReport)


###################################################
### code chunk number 26: ReportingTools-view (eval = FALSE)
###################################################
 browseURL(url)


