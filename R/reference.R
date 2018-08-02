library(SingleR)
library(genefilter)
library(dplyr)
source("../R/Seurat_functions.R")

#---------------------------
# E-GEOD-43717 (skip if GSE43717 is available)
#---------------------------
E_GEOD_43717 <- read.table(file = './data/GeneSets/E-GEOD-43717/E-GEOD-43717-query-results.tpms.tsv', sep = '\t', header = TRUE)
dup <- duplicated(E_GEOD_43717$Gene.Name)
E_GEOD_43717 <- E_GEOD_43717[!dup,]
rownames(E_GEOD_43717) <- E_GEOD_43717$Gene.Name
E_GEOD_43717 <- E_GEOD_43717[,-c(1:2)]
E_GEOD_43717[is.na(E_GEOD_43717)] <- 0
E_GEOD_43717 <- log1p(E_GEOD_43717)
boxplot(E_GEOD_43717)
name = 'Sertoli cells, spermatogonia, spermatocytes'
expr = as.matrix(E_GEOD_43717) # the expression matrix
types = gsub("\\."," ",colnames(E_GEOD_43717)) # a character list of the types. Samples from the same type should have the same name.
main_types = as.character(types) # a character list of the main types. 
ref_E_GEOD_43717 = list(name=name,data = expr, types=types, main_types=main_types)

# if using the de method, we can predefine the variable genes
ref_E_GEOD_43717$de.genes.main = CreateVariableGeneSet(expr,main_types,300)

# if using the sd method, we need to define an sd threshold
sd = rowSds(expr)
sd.thres = sort(sd, decreasing = T)[4000] # or any other threshold
ref_E_GEOD_43717$sd.thres = sd.thres

save(ref_E_GEOD_43717,file='./data/GeneSets/E-GEOD-43717/ref_E_GEOD_43717.RData') # it is best to name the object and the file with the same name.

############################
# GSE43717
############################
library(GEOquery)
library(dplyr)
eList <- getGEOSuppFiles("GSE43717",baseDir = paste0(getwd(),"/data/GeneSets")) #it takes some time
eList
list.files("./data/GeneSets/GSE43717")
untar("./data/GeneSets/GSE43717/GSE43717_RAW.tar", exdir = "./data/GeneSets/GSE43717/")
list.files("./data/GeneSets/GSE43717/")
FPKM_GSE43717_list <- grep("FPKM",list.files("./data/GeneSets/GSE43717/"),value = T)
FPKM_GSE43717_txt_list <-lapply(paste0("./data/GeneSets/GSE43717/",FPKM_GSE43717_list),read.table,sep = '\t', header = TRUE)
head(FPKM_GSE43717_txt_list[[1]])
FPKM_GSE43717_txt_list <-lapply(FPKM_GSE43717_txt_list,function(x) x[,c("gene_id","FPKM")])
FPKM_GSE43717_list <- gsub(".txt.gz","",FPKM_GSE43717_list)
for(i in 1:length(FPKM_GSE43717_txt_list)) colnames(FPKM_GSE43717_txt_list[[i]])[2] <- FPKM_GSE43717_list[i]
head(FPKM_GSE43717_txt_list[[1]])
FPKM_GSE43717_txt <- Reduce(merge,FPKM_GSE43717_txt_list)
head(FPKM_GSE43717_txt)
dim(FPKM_GSE43717_txt)
boxplot(FPKM_GSE43717_txt)

library(biomaRt)
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
filters = listFilters(ensembl)
filters[1:20,]
attributes = listAttributes(ensembl)
attributes[grepl("Gene name",attributes$description),]
attributes[grepl("Transcript stable ID",attributes$description),]
selected_attri <- c('ensembl_gene_id',
                  'ensembl_transcript_id',
                  'external_gene_name')
lapply(selected_attri,function(x) attributes[grepl(x,attributes$name),])
biomart_results <- getBM(attributes=selected_attri, 
              filters="ensembl_gene_id",
             values=FPKM_GSE43717_txt$gene_id,
      mart = ensembl)
colnames(biomart_results)[which(colnames(biomart_results)=="ensembl_gene_id")] <- "gene_id"
head(FPKM_GSE43717_txt)
head(biomart_results)
dim(FPKM_GSE43717_txt)
dim(biomart_results)
FPKM_GSE43717 <- inner_join(FPKM_GSE43717_txt,biomart_results,by="gene_id")
dim(FPKM_GSE43717)
dup <- duplicated(FPKM_GSE43717$external_gene_name)
FPKM_GSE43717 <- FPKM_GSE43717[!dup,]
dim(FPKM_GSE43717)
rownames(FPKM_GSE43717) <- FPKM_GSE43717$external_gene_name
FPKM_GSE43717 <- FPKM_GSE43717[,grep("FPKM",colnames(FPKM_GSE43717))]
boxplot(FPKM_GSE43717)
table(rowSums(FPKM_GSE43717)>1500000)
FPKM_GSE43717 <- FPKM_GSE43717[rowSums(FPKM_GSE43717)<1500000,] # remove 
dim(FPKM_GSE43717)
Log_FPKM_GSE43717 <- log1p(FPKM_GSE43717)
GSE43717_types <- c("Sertoli cells","Spermatogonia","Spermatocytes","Spermatids","Spermatozoa")
name = "GSE43717"

ref_GSE43717 = CreateSinglerReference(name = "GSE43717",
                                      expr = as.matrix(Log_FPKM_GSE43717),
                                      types = GSE43717_types, 
                                      main_types = GSE43717_types)
save(ref_GSE43717,file='./data/GeneSets/GSE43717/ref_GSE43717.RData') # it is best to name the object and the file with the same name.

############################
# GSE83264
############################
eList <- getGEOSuppFiles("GSE83264",baseDir = paste0(getwd(),"/data/GeneSets")) #it takes some time
eList
list.files("./data/GeneSets/GSE83264")
GSE83264 <- read.table(file = './data/GeneSets/GSE83264/GSE83264_mRNA_genes_merged_counts.txt.gz', sep = '\t', header = TRUE)
dup <- duplicated(GSE83264$gene)
table(dup)
GSE83264 <- GSE83264[!dup,]
rownames(GSE83264) <- GSE83264$gene
GSE83264 <- GSE83264[,-1]
GSE83264[is.na(GSE83264)] <- 0
boxplot(GSE83264)
dim(GSE83264)
head(GSE83264)
#convert counts to TPM
lname <- load("../SingleR/data/gene_lengths.RData")
lname
head(mouse_lengths)
length(mouse_lengths)
TPM_GSE83264 <- TPM(GSE83264, mouse_lengths)
Log_TPM_GSE83264 <- log1p(TPM_GSE83264)
par(mfrow=c(1,2))
boxplot(Log_TPM_GSE83264)
rownames(Log_TPM_GSE83264) = Hmisc::capitalize(tolower(rownames(Log_TPM_GSE83264)))

name = "whole testis at day 1-7 postpartum"
colnames(Log_TPM_GSE83264)
GSE83264_types = c("Day1 pp whole testis","Day1 pp whole testis","Day3 pp whole testis",
          "Day3 pp whole testis","Day7 pp whole testis","Day7 pp whole testis",
          "leptotene/zygotene spermatocytes","leptotene/zygotene spermatocytes",
          "pachytene spermatocytes","pachytene spermatocytes") # a character list of the types. Samples from the same type should have the same name.
GSE83264_main_types = as.character(GSE83264_types) # a character list of the main types. 

ref_GSE83264 = CreateSinglerReference(name = "GSE83264",
                                      expr = as.matrix(Log_TPM_GSE83264),
                                      types = GSE83264_types, 
                                      main_types = GSE83264_main_types)

save(ref_GSE83264,file='./data/GeneSets/GSE83264/ref_GSE83264.RData') # it is best to name the object and the file with the same name.

################################################
# merge hpca_blueprint_encode GSE43717 and GSE83264
################################################
# load hpca_blueprint_encode reference database
Iname = load(file='../SingleR/data/ref_Mouse.RData');Iname
Iname = load(file='./data/GeneSets/GSE43717/ref_GSE43717.RData');Iname
Iname = load(file='./data/GeneSets/GSE83264/ref_GSE83264.RData');Iname

#ref_GSE43717_GSE83264 <- reshape::merge_all(list(ref_immgen_mouse.rnaseq$data,
Ref_GSE43717 <- merge(ref_immgen_mouse.rnaseq$data, ref_GSE43717$data, by='row.names')
rownames(Ref_GSE43717) = Ref_GSE43717$Row.names
Ref_GSE43717 = Ref_GSE43717[,-1]
ref_GSE43717_GSE83264 <- merge(Ref_GSE43717, ref_GSE83264$data, by='row.names')
rownames(ref_GSE43717_GSE83264) = ref_GSE43717_GSE83264$Row.names
ref_GSE43717_GSE83264 = ref_GSE43717_GSE83264[,-1]

dim(ref_GSE43717_GSE83264)
testMMM(ref_GSE43717_GSE83264)
boxplot(ref_GSE43717_GSE83264) # too slow!!

types = c(ref_immgen_mouse.rnaseq$types,ref_GSE43717$types,ref_GSE83264$types)
main_types = c(ref_immgen_mouse.rnaseq$main_types,ref_GSE43717$main_types,
               ref_GSE83264$main_types)

Ref_GSE43717_GSE83264 = CreateSinglerReference(name = "immgen_mouse.rnaseq_GSE43717_GSE83264",
                                               expr = as.matrix(ref_GSE43717_GSE83264),
                                               types = types, 
                                               main_types = main_types)

save(Ref_GSE43717_GSE83264,file='./data/GeneSets/Ref_GSE43717_GSE83264.RData') # it is best to name the object and the file with the same name.

################################################
# merge GSE43717 and GSE83264
################################################
# load GSE43717 and GSE83264 reference database
Iname = load(file='./data/GeneSets/GSE43717/ref_GSE43717.RData');Iname
Iname = load(file='./data/GeneSets/GSE83264/ref_GSE83264.RData');Iname

ref_GSE43717_GSE83264 <- merge(ref_GSE43717$data, ref_GSE83264$data, by='row.names')
rownames(ref_GSE43717_GSE83264) = ref_GSE43717_GSE83264$Row.names
ref_GSE43717_GSE83264 = ref_GSE43717_GSE83264[,-1]

dim(ref_GSE43717_GSE83264)
testMMM(ref_GSE43717_GSE83264)
boxplot(ref_GSE43717_GSE83264)

ref_GSE43717_GSE83264 = CreateSinglerReference(name = "GSE43717_GSE83264",
                                           expr = as.matrix(ref_GSE43717_GSE83264),
                                           types = c(ref_GSE43717$types,ref_GSE83264$types), 
                                           main_types = c(ref_GSE43717$main_types,ref_GSE83264$main_types))

save(ref_GSE43717_GSE83264,file='./data/GeneSets/GSE43717_GSE83264.RData') # it is best to name the object and the file with the same name.
