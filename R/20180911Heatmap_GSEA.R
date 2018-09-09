########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
library(Seurat)
library(SingleR)
library(dplyr)
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
Spermato1 <- SubsetData(SSCs,ident.use = major_cells)
TSNEPlot.1(object = Spermato1,do.label = T, group.by = "orig.ident",
           do.return = TRUE, no.legend = T, #colors.use = singler.colors[4:6],
           pt.size = 1,label.size = 5,label.repel = T)+
        ggtitle("Spermatogonia and Spermatocytes at all time points")+
        theme(text = element_text(size=20),
              plot.title = element_text(hjust = 0.5, face = "bold"))

# 4.1.1 heatmap ==============
"check the timepoint of sub population in Spermatocytes"
Spermato_list <- list()
for(i in 1:length(major_cells)){
        Spermato_list[[i]] <- SubsetData(SSCs,ident.use = major_cells[i])
        Spermato_list[[i]]@meta.data$orig.ident = paste(major_cells[i],
                                                         Spermato_list[[i]]@meta.data$orig.ident)
}
Spermato <- Reduce(MergeSeurat,Spermato_list)
table(Spermato@meta.data$orig.ident)
# all time point ...1
Spermato@meta.data$orig.ident = gsub("Early Spermatocytes Ad-.*$",
                                     "Early Spermatocytes Adult",Spermato@meta.data$orig.ident)
Spermato@meta.data$orig.ident = gsub("Early Spermatocytes PND06|Early Spermatocytes PND14",
                                     "Early Spermatocytes PND06-14",Spermato@meta.data$orig.ident)
Spermato@meta.data$orig.ident = gsub("Early Spermatocytes PND25|Early Spermatocytes PND30",
                                     "Early Spermatocytes PND25-30",Spermato@meta.data$orig.ident)
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


Spermato <- SetAllIdent(object = Spermato, id = 'orig.ident')
table(Spermato@cell.names == Spermato1@cell.names)
Spermato@scale.data = Spermato1@scale.data
time_series <- c("PND06","PND14","PND18","PND25","PND30","Ad-depleteSp","Ad-Thy1")
# all time point ...0
time_points0 <- c(paste0("Spermatogonia ",time_series),
                  paste0("Early Spermatocytes ",time_series),
                  paste0("Spermatocytes ",time_series[-1]))

# all time point ...1
time_points1 <- c(paste0("Spermatogonia ",time_series),
                  paste0("Early Spermatocytes ", c("PND06-14","PND18","PND25-30","Adult")),
                  paste0("Spermatocytes ",time_series[-1]))
# merge time point ..2
time_points2 <- c("Spermatogonia PND06","Spermatogonia PND14","Spermatogonia PND18_and_later",
                 "Early Spermatocytes PND06-14","Early Spermatocytes PND18_and_later",
                 "Spermatocytes PND14","Spermatocytes PND18","Spermatocytes PND25-30",
                 "Spermatocytes Adult")

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
exprs_Spermato = AverageExpression(Spermato,use.scale = T)

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
MouseGene <- data.frame("NAME"=rownames(exprs_Spermato),stringsAsFactors = FALSE)# didn't works for [,1]
GSEA <- Mouse2Human(MouseGene)
colnames(GSEA) = c("MGI.symbol","NAME")
GSEA$DESCRIPTION = NA

exprs_Spermato$MGI.symbol = rownames(exprs_Spermato)
GSEA_exp <- merge(GSEA,exprs_Spermato,by = "MGI.symbol")
GSEA_exp = GSEA_exp[,-1]

path <- paste("./output",gsub("-","",Sys.Date()),sep = "/")
dir.create(path, recursive = T)
write.table(GSEA_exp, file= paste0(path,"/GSEA_exp.txt"), sep="\t", row.names = F)
# Insert NAME at [1,1] using Edited
# remove all " using Edited

# ===Create time_series phentype labels for GSEA ===
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



GSEA_exp <- read.delim(paste0(path,"/GSEA_exp.txt"), row.names=1)
gsea_Spermatogonia %>% kable() %>% kable_styling()



