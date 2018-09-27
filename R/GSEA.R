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
library(sva)
library(kableExtra)
source("./R/Seurat_functions.R")
source("./R/SingleR_functions.R")

# 5.1 load data ==============
lname1 = load(file = "./data/SSCs_20180825.Rda");lname1
SSCs@meta.data$orig.ident = gsub("PND18pre","PND18",SSCs@meta.data$orig.ident)
table(SSCs@meta.data$orig.ident)
table(SSCs@ident)
major_cells <- c("Spermatogonia","Early Spermatocytes","Spermatocytes")
#------- heatmap optional -----
Spermato1 <- SubsetData(SSCs,ident.use = major_cells)
TSNEPlot.1(object = Spermato1,do.label = T, group.by = "orig.ident",
           do.return = TRUE, no.legend = T, #colors.use = singler.colorss[4:6],
           pt.size = 1,label.size = 5,label.repel = T)+
        ggtitle("Spermatogonia and Spermatocytes at all time points")+
        theme(text = element_text(size=20),
              plot.title = element_text(hjust = 0.5, face = "bold"))
remove(Spermato1);GC()
# 5.2 filter, scale, merge time points ==============
"check the timepoint of sub population in Spermatocytes"
# 5.2.1 remove spermatid contaminated transcripts =======
#Spermatid_Genes_to_Filter <- readxl::read_excel("./doc/Spermatid Genes to Filter.xlsx")
#Spermatid_Genes = Spermatid_Genes_to_Filter$`Completely Spermatid Gene List To Remove (Spermatid/Other UMI Ratio > 20; and GSEA and tSNEs)`
#genes.use = rownames(SSCs@data); length(genes.use)
#table(Spermatid_Genes %in% genes.use)
#genes.use = genes.use[!(genes.use %in% Spermatid_Genes)]; length(genes.use)
#SSCs@raw.data = SSCs@raw.data[genes.use, ]

# 5.2.2 SubsetData, merge, batch correction again ===
Spermato_list <- list()
for(i in 1:length(major_cells)){
        Spermato_list[[i]] <- SubsetData(SSCs,ident.use = major_cells[i])
        Spermato_list[[i]]@meta.data$orig.ident = paste(major_cells[i],
                                                        Spermato_list[[i]]@meta.data$orig.ident)
}
Spermato <- Reduce(MergeSeurat,Spermato_list)
table(Spermato@meta.data$orig.ident)
dim(Spermato@data)
remove(SSCs,Spermato_list);GC()
# 5.2.3 batch correction again ===
Spermato_data <- Spermato@data
batch.effect=Spermato@meta.data$batch.effect
names(batch.effect) =  rownames(Spermato@meta.data)
m = as.matrix(Spermato@data)
m = m[rowSums(m)>0,]
Spermato@data = ComBat(m, batch.effect, prior.plots=FALSE, par.prior=TRUE)
Spermato <- ScaleData(object = Spermato,
                  model.use = "negbinom", do.par=T,
                  vars.to.regress = c("CC.Difference"),
                  display.progress = T)
Spermato@data = Spermato_data
remove(m,Spermato_data);GC()
Spermato@ident = factor(Spermato@ident,levels = time_points1)
table(Spermato@ident)
save(Spermato, file = "./output/Spermato_20180912.Rda")

# 5.2.4 merge # all time point ...1 =======
Spermato@meta.data$orig.ident = gsub("Early Spermatocytes Ad-.*$",
                                     "Early Spermatocytes Adult",Spermato@meta.data$orig.ident)
Spermato@meta.data$orig.ident = gsub("Early Spermatocytes PND06|Early Spermatocytes PND14",
                                     "Early Spermatocytes PND06-14",Spermato@meta.data$orig.ident)
Spermato@meta.data$orig.ident = gsub("Early Spermatocytes PND25|Early Spermatocytes PND30",
                                     "Early Spermatocytes PND25-30",Spermato@meta.data$orig.ident)
Spermato <- SetAllIdent(object = Spermato, id = 'orig.ident')
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
# 5.3 Create txt Expression Data for GSEA====

#---
exprs_Spermato = AverageExpression(Spermato,use.scale = T)
#---
gde.all <- FindAllMarkers(object = Spermato,test.use = "MAST",
                          logfc.threshold = 0.01,
                          min.pct = 0, 
                          min.cells.gene = 0,
                          return.thresh = 1)
path <- paste("./output",gsub("-","",Sys.Date()),sep = "/")
dir.create(path, recursive = T)
write.csv(gde.all,paste0(path,"/gde.all.csv"))

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
        MouseGene <- data.frame("NAME"= x$gene,stringsAsFactors = FALSE)# didn't works for [,1]
        
        require("biomaRt")
        human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
        
        genesV2 = getLDS(attributes = c("mgi_symbol"), #filters = "mgi_symbol",
                         values = MouseGene , mart = mouse, attributesL = c("hgnc_symbol"), 
                         martL = human, uniqueRows=T)
        rm = duplicated(genesV2[,1]) | duplicated(genesV2[,2])
        genesV2 = genesV2[!rm,]
        # Print the first 6 genes found to the screen
        colnames(genesV2) = c("gene","NAME")
        genesV2$DESCRIPTION = NA
        print(head(genesV2))
        GSEA_exp <- merge(genesV2,x,by = "gene")
        GSEA_exp = GSEA_exp[,-1]
        return(GSEA_exp)
}
GSEA_exp = Mouse2Human(exprs_Spermato)
GSEA_exp = GSEA_exp[,c("NAME","DESCRIPTION",time_points1)]
path <- paste("./output",gsub("-","",Sys.Date()),sep = "/")
dir.create(path, recursive = T)
write.table(GSEA_exp, file= paste0(path,"/GSEA_exp_SpermRemoved.txt"), sep="\t", row.names = F)
# Insert NAME at [1,1] using Edited
# remove all " using Edited

# 5.4 Create time_series phentype labels for GSEA ===
time.series <- c(6,14,18,25,30,35,40)
GSEA_cls <- list("#numeric",
                 paste0("#",major_cells[1]),
                 time.series,
                 paste0("#",major_cells[2]),
                 c(14,18,25,35),
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
gsea_report <- read.delim("~/Downloads/gsea_report_for_Spermatocytes_neg_1537120967457.txt", row.names=1)
gsea_report %>% kable() %>% kable_styling()
