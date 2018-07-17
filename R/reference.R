library(SingleR)
library(genefilter)

E_GEOD_43717 <- read.table(file = './data/GeneSets/E-GEOD-43717-query-results.tpms.tsv', sep = '\t', header = TRUE)
dup <- duplicated(E_GEOD_43717$Gene.Name)
E_GEOD_43717 <- E_GEOD_43717[!dup,]
rownames(E_GEOD_43717) <- E_GEOD_43717$Gene.Name
E_GEOD_43717 <- E_GEOD_43717[,-c(1:2)]
E_GEOD_43717[is.na(E_GEOD_43717)] <- 0
E_GEOD_43717 <- log1p(E_GEOD_43717)
name = 'Sertoli cells, spermatogonia, spermatocytes'
expr = as.matrix(E_GEOD_43717) # the expression matrix
types = gsub("\\."," ",colnames(E_GEOD_43717)) # a character list of the types. Samples from the same type should have the same name.
main_types = as.character(types) # a character list of the main types. 
ref_E_GEOD_43717 = list(name=name,data = expr, types=types, main_types=main_types)

# if using the de method, we can predefine the variable genes
ref_E_GEOD_43717$de.genes = CreateVariableGeneSet(expr,types,200)
ref_E_GEOD_43717$de.genes.main = CreateVariableGeneSet(expr,main_types,300)

# if using the sd method, we need to define an sd threshold
sd = rowSds(expr)
sd.thres = sort(sd, decreasing = T)[4000] # or any other threshold
ref_E_GEOD_43717$sd.thres = sd.thres

save(ref_E_GEOD_43717,file='./data/GeneSets/ref_E_GEOD_43717.RData') # it is best to name the object and the file with the same name.

source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)
listEnsembl()
Mouse = useEnsembl(biomart="ENSEMBL_MART_MOUSE")
