########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
library(Seurat)
library(SingleR)
library(dplyr)
library(tidyr)
library(MAST)
library(kableExtra)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")

# 4.1.1 load data ==============
lname1 = load(file = "./data/SSCs_20180825.Rda");lname1
SSCs@meta.data$orig.ident = gsub("PND18pre","PND18",SSCs@meta.data$orig.ident)
#SSCs@meta.data$orig.ident = gsub("Ad-.*$","Adult",SSCs@meta.data$orig.ident)
#SSCs@meta.data$orig.ident = gsub("PND25|PND30","PND25-30",SSCs@meta.data$orig.ident)
table(SSCs@meta.data$orig.ident)
table(SSCs@ident)
major_cells <- c("Spermatogonia","Early Spermatocytes","Spermatocytes")
#-------
Spermato1 <- SubsetData(SSCs,ident.use = major_cells)
TSNEPlot.1(object = Spermato1,do.label = T, group.by = "orig.ident",
           do.return = TRUE, no.legend = T, #colors.use = singler.colors[4:6],
           pt.size = 1,label.size = 5,label.repel = T)+
        ggtitle("Spermatogonia and Spermatocytes at all time points")+
        theme(text = element_text(size=20),
              plot.title = element_text(hjust = 0.5, face = "bold"))
remove(Spermato1);GC()
# 4.1.1 heatmap ==============
"check the timepoint of sub population in Spermatocytes"
Spermato_list <- list()
for(i in 1:length(major_cells)){
        Spermato_list[[i]] <- SubsetData(SSCs,ident.use = major_cells[i])
        Spermato_list[[i]]@meta.data$orig.ident = paste(major_cells[i],
                                                         Spermato_list[[i]]@meta.data$orig.ident)
}
Spermato <- Reduce(MergeSeurat,Spermato_list)
remove(SSCs,Spermato_list);GC()
table(Spermato@meta.data$orig.ident)
# all time point ...1
Spermato@meta.data$orig.ident = gsub("Early Spermatocytes Ad-.*$",
                                     "Early Spermatocytes Adult",Spermato@meta.data$orig.ident)
Spermato@meta.data$orig.ident = gsub("Early Spermatocytes PND06|Early Spermatocytes PND14",
                                     "Early Spermatocytes PND06-14",Spermato@meta.data$orig.ident)
Spermato@meta.data$orig.ident = gsub("Early Spermatocytes PND25|Early Spermatocytes PND30",
                                     "Early Spermatocytes PND25-30",Spermato@meta.data$orig.ident)
Spermato <- SetAllIdent(object = Spermato, id = 'orig.ident')
Spermato = ScaleData(Spermato)
time_series <- c("PND06","PND14","PND18","PND25","PND30","Ad-depleteSp","Ad-Thy1")
# all time point ...1
time_points1 <- c(paste0("Spermatogonia ",time_series),
                  paste0("Early Spermatocytes ", c("PND06-14","PND18","PND25-30","Adult")),
                  paste0("Spermatocytes ",time_series[-1]))

#-------
All_markers = read.csv("./output/20180826/All_markers.csv", header = T,stringsAsFactors =F)
Spermato_markers <- All_markers[All_markers$major_cells %in% major_cells,] %>% 
        mutate(cluster = factor(major_cells))
Spermato_markers %>% head(20) %>% kable() %>% kable_styling()
DoHeatmap.1(Spermato,Spermato_markers,Top_n = 25, 
            group.order = time_points1, 
            ident.use = "Spermatogonia and Spermatocytes",
            group.label.rot = T,cex.row = 7,remove.key =T,title.size = 12)


# ===Create txt Expression Data for GSEA====
Spermato@ident = factor(Spermato@ident,levels = time_points1)
table(Spermato@ident)
#---
exprs_Spermato = AverageExpression(Spermato,use.scale = T)
#---
gde.all <- FindAllMarkers(object = Spermato,test.use = "MAST",
                          logfc.threshold = 0.001,
                          min.pct = 0, 
                          min.cells.gene = 0,
                          return.thresh = 1)
gde_Spermato <- gde.all[,c("avg_logFC","cluster","gene")]
exprs_Spermato <- spread(gde_Spermato,key = cluster, value = avg_logFC)
exprs_Spermato[is.na(exprs_Spermato)] = 0
#---
exprs_Spermato = Spermato@scale.data
idents <- as.data.frame(table(Spermato@ident))
time.rep <- as.integer(idents$Freq)
exprs_Spermato = exprs_Spermato[,1:(time.rep[1]+time.rep[2])]
exprs_Spermato = cbind.data.frame(rownames(exprs_Spermato),exprs_Spermato)
colnames(exprs_Spermato)[1] = "gene"
# Convert Mounse Gene names to Human gene for GSEA analysis
# https://www.r-bloggers.com/converting-mouse-to-human-gene-names-with-biomart-package/
Mouse2Human <- function(x, unique = FALSE){
        
        require("biomaRt")
        human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
        
        genesV2 = getLDS(attributes = c("mgi_symbol"), #filters = "mgi_symbol",
                         values = x , mart = mouse, attributesL = c("hgnc_symbol"), 
                         martL = human, uniqueRows=T)
        rm = duplicated(genesV2[,1]) | duplicated(genesV2[,2])
        genesV2 = genesV2[!rm,]
        # Print the first 6 genes found to the screen
        print(head(genesV2))
        return(genesV2)
}
MouseGene <- data.frame("NAME"= exprs_Spermato$gene,stringsAsFactors = FALSE)# didn't works for [,1]
GSEA <- Mouse2Human(MouseGene)
colnames(GSEA) = c("gene","NAME")
GSEA$DESCRIPTION = NA
GSEA_exp <- merge(GSEA,exprs_Spermato,by = "gene")
GSEA_exp = GSEA_exp[,-1]
path <- paste("./output",gsub("-","",Sys.Date()),sep = "/")
dir.create(path, recursive = T)
write.table(GSEA_exp, file= paste0(path,"/Spermatogonia_6_vs_14_zscore.txt"), sep="\t", row.names = F)
# Insert NAME at [1,1] using Edited
# remove all " using Edited
Spermatogonia_logFC <- GSEA_exp[,1:9]
write.table(Spermatogonia_logFC, file= paste0(path,"/Spermatogonia_logFC.txt"), sep="\t", row.names = F)
EarlySpermatocytes_logFC <- GSEA_exp[,c(1:2,10:13)]
write.table(EarlySpermatocytes_logFC, file= paste0(path,"/EarlySpermatocytes_logFC.txt"), sep="\t", row.names = F)
Spermatocytes_logFC <- GSEA_exp[,c(1:2,14:19)]
write.table(Spermatocytes_logFC, file= paste0(path,"/Spermatocytes_logFC.txt"), sep="\t", row.names = F)


# ===Create time_series phentype labels for GSEA ===
time.series <- c(6,14,18,25,30,35,40)
GSEA_cls <- list("#numeric",
                 paste0("#",major_cells[1]),
                 time.series,
                 paste0("#",major_cells[2]),
                 c("6-14",18,25-30,35),
                 paste0("#",major_cells[3]),
                 time.series[-1])
#GSEA.cls = lapply(GSEA.cls, function(x) sub(" ","_",x))
GSEA_cls <- t(list2df(GSEA_cls))
GSEA_cls[is.na(GSEA_cls)] <- ""
#GSEA_cls = noquote(GSEA_cls)
write.table(GSEA_cls,file= paste0(path,"/GSEA_exp.cls"), sep=" ")
# open with texteditor
# delete the first row;
# delete the first column;
# delete all NA and "
# delete the space before and after each row.



gsea_report <- read.delim(paste0(path,"/gsea_report_for_Spermatocytes_pos_1536630312315.txt"), row.names=1)
gsea_report %>% kable() %>% kable_styling()


# ===Create regular phentype labels for GSEA ===
sample_n <- ncol(exprs_Spermato)-1
classes_n <- 2
GSEA.cls <- list(c(sample_n,classes_n,1),
                 c("#","Spermatogonia_PND06", "Spermatogonia_PND14"),
                 c(rep(1,time.rep[1]),
                   rep(2,time.rep[2])))
GSEA.cls = lapply(GSEA.cls, function(x) sub(" ","",x))
GSEA_cls <- t(list2df(GSEA.cls))
GSEA_cls[is.na(GSEA_cls)] <- ""
write.table(GSEA_cls,file= paste0(path,"/Spermatogonia_6_vs_14.cls"), sep=" ")


# merge time point ..2
Spermato@meta.data$orig.ident = gsub("Spermatogonia PND18|Spermatogonia PND25|
                                     |Spermatogonia PND30|Spermatogonia Ad-depleteSp|Spermatogonia Ad-Thy1",
                                     "Spermatogonia PND18_and_later",Spermato@meta.data$orig.ident)
Spermato@meta.data$orig.ident = gsub("Early Spermatocytes PND18|Early Spermatocytes PND25-30|
                                     |Early Spermatocytes Adult",
                                     "Early Spermatocytes PND18_and_later",Spermato@meta.data$orig.ident)
Spermato@meta.data$orig.ident = gsub("Spermatocytes PND25|Spermatocytes PND30",
                                     "Spermatocytes PND25-30",Spermato@meta.data$orig.ident)
Spermato@meta.data$orig.ident = gsub("Spermatocytes Ad-depleteSp|Spermatocytes Ad-Thy1",
                                     "Spermatocytes Adult",Spermato@meta.data$orig.ident)


# merge time point ..2
time_points2 <- c("Spermatogonia PND06","Spermatogonia PND14","Spermatogonia PND18_and_later",
                  "Early Spermatocytes PND06-14","Early Spermatocytes PND18_and_later",
                  "Spermatocytes PND14","Spermatocytes PND18","Spermatocytes PND25-30",
                  "Spermatocytes Adult")


Spermatogonia_logFC <- GSEA_exp[,1:9]
write.table(Spermatogonia_logFC, file= paste0(path,"/Spermatogonia_logFC.txt"), sep="\t", row.names = F)
EarlySpermatocytes_logFC <- GSEA_exp[,c(1:2,10:13)]
write.table(EarlySpermatocytes_logFC, file= paste0(path,"/EarlySpermatocytes_logFC.txt"), sep="\t", row.names = F)
Spermatocytes_logFC <- GSEA_exp[,c(1:2,14:19)]
write.table(Spermatocytes_logFC, file= paste0(path,"/Spermatocytes_logFC.txt"), sep="\t", row.names = F)



# ===Create regular phentype labels for GSEA ===
sample_n <- ncol(exprs_Spermato)-1
classes_n <- 2
GSEA.cls <- list(c(sample_n,classes_n,1),
                 c("#","Spermatogonia_PND06", "Spermatogonia_PND14"),
                 c(rep(1,time.rep[1]),
                   rep(2,time.rep[2])))
GSEA.cls = lapply(GSEA.cls, function(x) sub(" ","",x))
GSEA_cls <- t(list2df(GSEA.cls))
GSEA_cls[is.na(GSEA_cls)] <- ""
write.table(GSEA_cls,file= paste0(path,"/Spermatogonia_6_vs_14.cls"), sep=" ")
### FUNCTON TO COMPUTE PATHWAY Z-SCORE --------------------------------------------------------------
# Usage:
# get.pathway.zcore(file.expr, file.pathway)
#
# Description:
# file.expr = path of file containing gene-expression matrix with Row as Genes and Columns as Samples.
# file.pathway = path of file containing pathway-gene relationship
#	- each line should contain information about a unique pathway
# 	- first-column must contain name of the pathway
#	- second-column must contain genes mapped to the corresponding pathway 
#		multiple genes must be separated by comma symbol (,)
### --------------------------------------------------------------------------------------------------
get.pathway.zcore <- function(file.expr, file.pathway){
        # LOAD LIBRARIES -----
        require("stringr")
        
        # LOAD Expression ----
        expr <- read.delim(file.expr, header=T, row.names=1, stringsAsFactors=F)
        
        # LOAD Pathway ------
        dat.path <- read.delim(file.pathway, header=T, stringsAsFactors=F)
        
        # Get Zscores -------
        pathz <- matrix(0, nrow=nrow(dat.path), ncol=4, dimnames=list(dat.path[,1], colnames(expr)[-1]))
        for(i in 1:nrow(dat.path)){
                pathway <- dat.path[i,1]
                genes <- str_split(dat.path[i,2], ",")[[1]]
                
                expr.path <- subset(expr, rownames(expr) %in% genes)[,-1]
                exprz <- getZscore(expr.path)
                pathz[i,] <- apply(exprz, 2, mean)
        }
        return(pathz)
}


### FUNCTION TO CALCULATE Z-SCORE ----------------------------------------
# Let, M be the gene-expression matrix 
# with Row as Genes and Columns as Samples.
#
# Zscore = (x_ij - m_i) / sd_i
# Here,
# i = gene
# j = sample
# x_ij = expression value of gene i in sample j
# m_i = mean of expression values in for gene i across all j
# sd_i = standard-deviation of expression values in for gene i across all j
#
# Usage:
# getZscore(dat)
#
# Description:
# dat = gene-expression matrix with Row as Genes and Columns as Samples.
### ------------------------------------------------------------------------
getZscore <- function(dat){
        z.dat <- matrix(nrow=nrow(dat), ncol=ncol(dat), dimnames=list(rownames(dat),colnames(dat)))
        
        dat.mean <- apply(dat, 1, mean)
        dat.stdev <- apply(dat, 1, sd)
        
        for(i in 1:nrow(dat)){
                x <- as.numeric(dat[i,])
                z.dat[i,] <- (x - dat.mean[i])/dat.stdev[i]
        }
        z.dat <- as.data.frame(z.dat)
        return(z.dat)
}
