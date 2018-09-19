########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(kableExtra)
library(SingleR)
library(reshape2)
library(MAST)
source("./R/Seurat_functions.R")

path <- paste("./output",gsub("-","",Sys.Date()),sep = "/")
dir.create(path, recursive = T)
# 4.1 load data & Rename ident, Compare DE across all major cell types=========================
# 4.1.1 load data ==============
lname1 = load(file = "./data/SSCs_20180825.Rda");lname1
SSCs@meta.data$orig.ident = gsub("PND18pre","PND18",SSCs@meta.data$orig.ident)
table(SSCs@meta.data$orig.ident)
table(SSCs@ident)

TSNEPlot.1(object = SSCs,do.label = T, group.by = "ident",
            do.return = TRUE, no.legend = T,colors.use = singler.colors,
            pt.size = 1,label.size = 5,label.repel = T)+
            ggtitle("Manually label cell types")+
            theme(text = element_text(size=20),
            plot.title = element_text(hjust = 0.5, face = "bold"))

# 4.1.2 Compare DE across all cell types ==============
All_markers <- FindAllMarkers.UMI(SSCs,test.use = "MAST")
kable(All_markers[1:20,]) %>% kable_styling()
rownames(All_markers) = NULL
All_markers <- All_markers %>% select("gene", everything()) # Moving the last column to the start

major_cells <- c("Spermatogonia","Early Spermatocytes","Spermatocytes",
                 "Round Spermatids","Spermatids","Sertoli cells",
                 "Endothelial & Hematopoietic cells", "Smooth muscle")
top <- group_top_mutate(df = All_markers, major_cells)
top %>% head(20) %>% kable() %>% kable_styling()
DoHeatmap.1(SSCs,top,Top_n = 15, 
            group.order = major_cells,ident.use = "all cell types",
            group.label.rot = T,cex.row = 5,remove.key =T)
top <- select(top,-cluster)
write.csv(top,"./output/20180826/All_markers.csv")

# 4.1.3 GO enrichment analysis ==============
All_markers = read.csv("./output/20180826/All_markers.csv", header = T,stringsAsFactors =F)
Spermatogonia_markers <- All_markers[All_markers$major_cells == "Spermatogonia",]
Spermatogonia_markers <- Spermatogonia_markers[abs(Spermatogonia_markers$avg_logFC) > 0.75,"gene"]
length(unique(Spermatogonia_markers))
library(biomaRt)
# collect gene names from biomart
mart <- biomaRt::useMart("ensembl",dataset="mmusculus_gene_ensembl")
# Get ensembl gene ids and GO terms
getBM <- biomaRt::getBM(attributes = c("external_gene_name", "go_id"),
                        filters= "external_gene_name",
                        values = rownames(SSCs@data),
                        mart = mart)

Spermatogonia.Go <-  topGOterms(int.genes = unique(Spermatogonia_markers),
                                bg.genes = rownames(SSCs@data),
                                organism =  "Mouse",
                                getBM = getBM)
Spermatogonia.Go %>% kable() %>% kable_styling()
#====== 4.2 SubsetData, Compare DE across all major cell types==================

# 4.2.1 SSCs-spermato ==============
SSCs_spermato <- SubsetData(SSCs,ident.use = major_cells[1:5])
# 4.2.1.1 remove spermatid contaminated transcripts =======
Spermatid_Genes_to_Filter <- readxl::read_excel("./doc/Spermatid Genes to Filter.xlsx")
Spermatid_Genes = Spermatid_Genes_to_Filter$`Completely Spermatid Gene List To Remove (Spermatid/Other UMI Ratio > 20; and GSEA and tSNEs)`
genes.use = rownames(SSCs_spermato@data); length(genes.use)
table(Spermatid_Genes %in% genes.use)
genes.use = genes.use[!(genes.use %in% Spermatid_Genes)]; length(genes.use)
SSCs_spermato@raw.data = SSCs_spermato@raw.data[genes.use, ]
SSCs_spermato@data = SSCs_spermato@data[genes.use, ]
#SSCs_spermato@scale.data = SSCs_spermato@scale.data[genes.use, ]
TSNEPlot.1(object = SSCs_spermato,do.label = T, group.by = "ident",
           do.return = TRUE, no.legend = T,colors.use = singler.colors[4:8],
           pt.size = 1,label.size = 5,label.repel = T)+
        ggtitle("Five germ cell types")+
        theme(text = element_text(size=20),
              plot.title = element_text(hjust = 0.5, face = "bold"))

SSCs_spermato_markers <- top[(top$cluster %in% major_cells[1:5]),]
DoHeatmap.1(SSCs_spermato,SSCs_spermato_markers,Top_n = 20, 
            group.order = major_cells[1:5],ident.use = "all germ cells",
            group.label.rot = T,cex.row = 6,remove.key =T)

SSCs_spermato_markers <- FindAllMarkers.UMI(SSCs_spermato,test.use = "MAST")
kable(SSCs_spermato_markers[1:20,]) %>% kable_styling()
rownames(SSCs_spermato_markers) = NULL
SSCs_spermato_markers <- SSCs_spermato_markers %>% select("gene", everything()) # Moving the last column to the start
write.csv(SSCs_spermato_markers,"./output/20180918/SSCs_spermato_markers.csv")

# 4.2.0. define pipeline ==============

#' pipeline renameIdent + FindAllMarkers + group_top + write.csv
#' @param object A Seurat object
#' @param ident.use Create a cell subset based on the provided identity classes, inherited from Seurat::SubsetData
#' @param ... Name-value pairs of expressions, inherited from ... in dplyr::mutate
#' @param merge TRUE/FALSE, merge orig.ident or not
#' @param replace.from character, orig.ident before merge
#' @param replace.to character, orig.ident after merge
#' @export top a re-arranged data frame, sorted by avg_logFC, then arrange by ...
#' @example 
#' Spermatogonia_list <- RenameIdent_FindAllMarkers(object = Spermatogonia,
#'                       ident.use = "Spermatogonia",
#'                       time_points,
#'                       merge = T,
#'                       replace.from = "PND18|PND25|PND30",
#'                       replace.to = "PND18-30",
#'                       Top_n=300)
RenameIdent_FindAllMarkers <- function(object,ident.use, ...,
                                       replace.from = NULL,
                                       replace.to,Top_n=300){
        # merge ident
        if(!is.null(replace.from)) {
                object@meta.data$orig.ident = sub(replace.from,
                                                  replace.to,
                                                  object@meta.data$orig.ident)
        }
        # Top_n = NULL, save all genes
        if(is.null(Top_n)) Top_n = nrow(object@data)
        # rename ident
        new.levels = deparse(substitute(...))
        new.order = assign(new.levels,...)
        ident.vector <- as.factor(object@meta.data$orig.ident)
        names(x = ident.vector) <- names(object@ident)
        object@ident = ident.vector
        object@ident <- factor(x = object@ident,
                               levels = new.order) 
        print(table(object@ident))
        
        g <- TSNEPlot(object = object,do.label = F, group.by = "ident", 
                      do.return = TRUE, no.legend = F, 
                      pt.size = 1,label.size = 8 )+
                ggtitle(ident.use)+
                theme(text = element_text(size=20),     							
                      plot.title = element_text(hjust = 0.5,size = 25, face = "bold")) 
        print(g)
        
        # FindAllMarkers
        object_develop <- FindAllMarkers.UMI(object,test.use = "MAST")
        object_develop <- group_top_mutate(df = object_develop, time_points,Top_n = Top_n)
        object_develop_remove_cluster <- object_develop %>% subset(select = -cluster)
        print(head(object_develop))
        path <- paste("./output",gsub("-","",Sys.Date()),sep = "/")
        dir.create(path, recursive = T)
        write.csv(object_develop_remove_cluster,file= paste0(path,"/",ident.use,"_develop.csv"))
        return(list(object, object_develop, g))
}

# 4.2.2 Spermatogonia ==============
Spermatogonia <- SubsetData(SSCs,ident.use = "Spermatogonia")
table(Spermatogonia@meta.data$orig.ident)
time_points <- c("PND06","PND14","PND18","PND25-30","Ad-depleteSp","Ad-Thy1")
Spermatogonia_list <- RenameIdent_FindAllMarkers(object = Spermatogonia,
                                                ident.use = "Spermatogonia",
                                                time_points,
                                                replace.from = "PND25|PND30",
                                                replace.to = "PND25-30",
                                                Top_n=300)
Spermatogonia_list[[3]]
table(Spermatogonia_list[[1]]@ident) %>% t() %>% kable() %>% kable_styling()
Spermatogonia_list[[2]] %>% group_by(cluster) %>% top_n(5,avg_logFC)%>% kable() %>% kable_styling()
DoHeatmap.1(Spermatogonia_list[[1]],Spermatogonia_list[[2]],Top_n = 15, 
            group.order = time_points,ident.use = "Spermatogonia",
            group.label.rot = T,cex.row = 8,remove.key =F)
VlnPlot(object = Spermatogonia_list[[1]], features.plot = c("percent.mito"), nCol = 1,
        group.by = "orig.ident", x.lab.rot = T, do.return = T)

# add genes for 8/24 email
marker = MouseGenes(Spermatogonia,c("Gfra1", "Zbtb16", "Nanos2", "Utf1", "Sall4",
                                    "Id4", "Sohlh1", "Kit", "Lin28a", "Pax7", "Dmrt1"))
top <-  Spermatogonia_list[[2]] %>% group_by(cluster) %>% top_n(15, avg_logFC)
DoHeatmap(object = Spermatogonia_list[[1]], 
          genes.use = c(marker,"Odf2", top$gene), slim.col.label = TRUE,
          group.order = time_points, group.label.rot = T,cex.row = 6,remove.key =T) +
        ggtitle(paste("Expression heatmap of top",15,
                      "differential expression genes in each time points for",
                      "Spermatogonia"))+
        theme(text = element_text(size=20),     							
              plot.title = element_text(hjust = 0.5,size=14))


# 4.2.3 Early_spermatocytes ==============
Early_Spermatocytes <- SubsetData(SSCs,ident.use = "Early Spermatocytes")
table(Early_Spermatocytes@meta.data$orig.ident)
cells.use = Early_Spermatocytes@cell.names[!(Early_Spermatocytes@meta.data$orig.ident %in% c("Ad-Thy1",
                                                                                             "PND06"))]
Early_Spermatocytes <- SubsetData(SSCs,cells.use = cells.use)
table(Early_Spermatocytes@meta.data$orig.ident)
time_points <- c("PND14","PND18","PND25-30","Ad-depleteSp")
Early_spermatocytes_list <- RenameIdent_FindAllMarkers(Early_Spermatocytes,
                                                       ident.use = "Early spermatocytes",
                                                       time_points,
                                                       replace.from = "PND25|PND30",
                                                       replace.to = "PND25-30",
                                                       Top_n=300)
table(Early_spermatocytes_list[[1]]@ident) %>% t() %>% kable() %>% kable_styling()
DoHeatmap.1(Early_spermatocytes_list[[1]],Early_spermatocytes_list[[2]],Top_n = 20, 
            group.order  = time_points,
            ident.use = "Early spermatocytes")


# 4.2.4 Spermatocytes ==============
Spermatocytes <- SubsetData(SSCs,ident.use = "Spermatocytes")
table(Spermatocytes@meta.data$orig.ident)
time_points <- c("PND14","PND18","PND25","PND30","Ad-depleteSp","Ad-Thy1")
Spermatocytes_list <- RenameIdent_FindAllMarkers(Spermatocytes,
                                                ident.use = "Spermatocytes",
                                                time_points,
                                                Top_n=400)
DoHeatmap.1(Spermatocytes_list[[1]],Spermatocytes_list[[2]],Top_n = 10, 
            group.order = time_points,ident.use = "Spermatocytes",title.size = 14)
table(Spermatocytes_list[[1]]@ident) %>% t() %>% kable() %>% kable_styling()


#  Spermatids-----
Spermatids <- SubsetData(SSCs,ident.use = "Spermatids")
table(Spermatids@meta.data$orig.ident)
time_points <- c("PND25_and_early","PND30","Ad-depleteSp","Ad-Thy1")
Spermatids_list <- RenameIdent_FindAllMarkers(Spermatids,
                                              ident.use <- "Spermatids",
                                              time_points,
                                              merge = T,
                                              replace.from = "PND18|PND25",
                                              replace.to = "PND25_and_early")
DoHeatmap.1(Spermatids_list[[1]],Spermatids_list[[2]],Top_n = 10, 
            group.order = time_points,
            ident.use = "Spermatids")
table(Spermatocytes_list[[1]]@ident) %>% t() %>% kable() %>% kable_styling()

# Step 3. find overlap gene===========================
Spermatogonia_develop = read.csv("./output/20180812/Spermatogonia_develop.csv", header = TRUE)
kable(Spermatogonia_develop[1:20,]) %>% kable_styling()
Spermatocytes_develop = read.csv("./output/20180812/Spermatocytes_develop.csv", header = TRUE)
kable(Spermatocytes_develop[1:20,]) %>% kable_styling()
All_markers = read.csv("./output/20180812/All_markers.csv", header = TRUE)
kable(All_markers[1:20,]) %>% kable_styling()

intersect_DoHeatmap <- function(object, markers, dev_markers, Top_n = 10, 
                                Time_points,ident.use,...){
        if(!is.null(ncol(markers))) { # df
                subset_markers <- markers[markers$cluster == ident.use,]
        } else if(is.null(ncol(markers))){ # vector
                subset_markers = data.frame(gene = markers)
        }
        print(table(dev_markers$gene %in% subset_markers$gene))
        subset_dev_mark_genes <- intersect(dev_markers$gene, 
                                           subset_markers$gene)
        subset_dev_mark <- dev_markers[(dev_markers$gene %in% 
                                                subset_dev_mark_genes),]
        write.csv(subset_dev_mark,file= paste0("./output/20180812/",
                                               ident.use,"_dev_mark.csv"))
        DoHeatmap.2(object = object, marker_df = subset_dev_mark,
                    Top_n = Top_n, Time_points = Time_points,
                    ident.use = ident.use, ...)
}
# Spermatogonia
intersect_DoHeatmap(Spermatogonia_list[[1]], All_markers, Spermatogonia_list[[2]],
                    Top_n = 20, c("PND06","PND14","PND18_and_later"),
                    ident.use = "Spermatogonia",group.label.rot =F)

# Early_spermatocytes
intersect_DoHeatmap(object = Early_spermatocytes_list[[1]], All_markers, 
                    dev_markers = Early_spermatocytes_list[[2]],
                    Top_n = 30, Time_points = c("PND14","PND18_and_later"),
                    ident.use = "Early spermatocytes",group.label.rot =F)

# Spermatocytes
intersect_DoHeatmap(object = Spermatocytes_list[[1]], All_markers, 
                    dev_markers = Spermatocytes_list[[2]],
                    Top_n = 15, Time_points = c("PND18","PND25","PND30","Ad-depleteSp","Ad-Thy1"),
                    ident.use = "Spermatocytes",group.label.rot =T)


# FeatureHeatmap==========
FeatureHeatmap_save <- function(object, gene.list, save.folder, Top_n= 100){
        if(!is.null(nrow(gene.list))){
                gene.list_top <- gene.list %>% 
                        group_by(cluster) %>% 
                        top_n(Top_n, wt = avg_logFC) 
                gene.list_top <- as.character(gene.list_top$gene)
        } else if(is.character(gene.list)) {
                gene.list_top <- na.omit(gene.list[1:Top_n])
        }
        path <- paste("./output",gsub("-","",Sys.Date()),"tsne_matrix",save.folder,sep = "/")
        dir.create(path, recursive = T)
        for (gene in gene.list_top) {
                temp_plot = FeatureHeatmap.1(object =object, gene)+
                        theme(text = element_text(size=15))
                ggsave(temp_plot, file=paste0(path,"/",gene,".jpeg"),
                       width = 35.28, height = 24.69, units = "cm"
                )
        }
}

Spermatogonia_dev_mark = read.csv("./output/20180812/Spermatogonia_dev_mark.csv", header = TRUE)
kable(Spermatogonia_dev_mark[1:20,]) %>% kable_styling()
FeatureHeatmap_save(object = Spermato, gene.list = Spermatogonia_dev_mark,
                    save.folder = "Spermatogonia")

Early_spermatocytes_dev_mark = read.csv("./output/20180812/Early_spermatocytes_dev_mark.csv", header = TRUE)
kable(Early_spermatocytes_dev_mark[1:20,]) %>% kable_styling()
FeatureHeatmap_save(object = Spermato, gene.list = Early_spermatocytes_dev_mark,
                    save.folder = "Early spermatocytes")

Spermatocytes_dev_mark = read.csv("./output/20180812/Spermatocytes_dev_mark.csv", header = TRUE)
kable(Spermatocytes_dev_mark[1:20,]) %>% kable_styling()
FeatureHeatmap_save(object = Spermato, gene.list = Spermatocytes_dev_mark,
                    save.folder = "Spermatocytes")

#Clullin Ring ubiquitin ligase
DCAFs <- read.csv("./output/Paula20180813/DCAFs.csv",stringsAsFactors=F)
F_Box <- read.csv("./output/Paula20180813/F_Box_Protein_IDs.csv",stringsAsFactors=F)
Cullin_proteins <- read.csv("./output/Paula20180813/Cullin_proteins_and_others.csv",
                            stringsAsFactors=F)

DCAFs_genes <- MouseGenes(SSCs_Spermato,c(DCAFs$Symbol,F_Box$Symbol,
                                          Cullin_proteins$Symbol),
                          unique = T)

dev_order <- c("Sertoli cells","Spermatogonia","Early spermatocytes","Spermatocytes",
               "Spermatids")
All_markers <- All_markers[All_markers$avg_logFC >0,]
intersect_DoHeatmap(object = SSCs_Spermato, markers = DCAFs_genes, 
                    dev_markers = All_markers,
                    Top_n = 15, Time_points = dev_order,
                    ident.use = "",group.label.rot =T,cex.row = 8)
#2018/08/16
Spermato_idents <- c("Spermatogonia","Early spermatocytes","Spermatocytes")
Spermato <- SubsetData(SSCs, ident.use = Spermato_idents)
Spermato@meta.data$orig.ident <- gsub("Ad-","zAd-",Spermato@meta.data$orig.ident)
FeatureHeatmap_save(object = Spermato, gene.list = c("Asrgl1", "Tbpl1", "Dmrtb1",
"Rbmxl2"),save.folder = "Spermatogonia_spermatocytes")
