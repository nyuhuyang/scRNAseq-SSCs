#https://www.biostars.org/p/250927/
library(biomaRt)
# collect gene names from biomart
mart <- biomaRt::useMart("ensembl",dataset="mmusculus_gene_ensembl")
# Get ensembl gene ids and GO terms
GTOGO <- biomaRt::getBM(attributes = c("external_gene_name", "go_id"),
                        filters= "external_gene_name",
                        values = bg.genes,
                        mart = mart)
#examine result
head(GTOGO);dim(GTOGO)
#Remove blank entries
GTOGO <- GTOGO[GTOGO$go_id != '',]
# convert from table format to list format
geneID2GO <- by(GTOGO$go_id, GTOGO$external_gene_name, function(x) as.character(x))
# examine result
head(geneID2GO)


length(GTOGO$external_gene_name)
all.genes <- sort(unique(as.character(GTOGO$external_gene_name)))
length(all.genes)
fg.genes <- unique(Spermatogonia_markers) # Takes all the unique cell type specific genes
int.genes <- factor(as.integer(all.genes %in% fg.genes))
names(int.genes) = all.genes

go.obj <- new("topGOdata", ontology='BP',
              allGenes = int.genes,
              annot = annFUN.gene2GO,
              gene2GO = geneID2GO,
              nodeSize = 10
)
resultFisher <- runTest(go.obj, algorithm = "elim", statistic = "fisher")
resultKS <- runTest(go.obj, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(go.obj, algorithm = "elim", statistic = "ks")

results.tab <- GenTable(go.obj, classicFisher = resultFisher, 
                        classicKS = resultKS, elimKS = resultKS.elim,topNodes = 20,
                        orderBy = "elimKS", ranksOf = "classicFisher")
results.tab %>% kable() %>% kable_styling()
