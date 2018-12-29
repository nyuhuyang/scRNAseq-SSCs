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
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")

#4.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 4.1.1 load data
lname1 = load(file = "./data/SSCs_suplabel_GSE43717.Rda") ; lname1
table(SSCs@ident)
table(SSCs@ident,SSCs@meta.data$orig.ident) %>% kable() %>%  kable_styling()
g <- TSNEPlot.1(object = SSCs, no.legend = F, do.label = T, pt.size = 1,
         do.return = TRUE, label.repel = T,label.size = 4,colors.use = singler.colors)+
        ggtitle("Manually assign cell types based on marker genes")+
        theme(text = element_text(size=20),
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 
# 08/12 Manual marked ======================================
SingleFeaturePlot.1(object = SSCs, "Zbtb16", threshold = 0.1)
SingleFeaturePlot.1(object = SSCs, "Kit", threshold = 0.1)

# select cell manually.......
Spermatogonia_cells1 <- FeaturePlot(Spermatogonia, features.plot = "Zbtb16",pt.size = 2,
                             cols.use = c("orange","darkblue"), do.identify = T) # manually select
early_spermatocytes <- FeaturePlot(SSCs, features.plot = "Kit",pt.size = 2,
                                   cols.use = c("orange","darkblue"), do.identify = T) # manually select
Sertoli <- FeaturePlot(SSCs, features.plot = "Wt1",pt.size = 2,
                      cols.use = c("orange","darkblue"), do.identify = T) # manually select
Hbb_bt <- FeaturePlot(SSCs, features.plot = "Hbb-bt",pt.size = 2,
                                   cols.use = c("orange","darkblue"), do.identify = T) # manually select
Macrophages <- FeaturePlot(SSCs_Spermato, features.plot = "Laptm5",pt.size = 2,
                       cols.use = c("orange","darkblue"), do.identify = T) # manually select

table(Spermatogonia %in% early_spermatocytes)
table(Hbb_bt %in% early_spermatocytes)
early_spermatocytes = early_spermatocytes[-which(early_spermatocytes %in% Hbb_bt)]
Sertoli2 = read.csv("./output/20180812/remove_Sertoli_20180809.csv",header = 1, row.names = 1)
Sertoli = unique(c(Sertoli,as.vector(Sertoli$x)))
table(Hbb_bt %in% as.vector(Sertoli))

# renameidnet manually.......
Cell_names_20180812 <- list("Spermatogonia" = Spermatogonia_cells,
                            "Sertoli cells" = Sertoli,
                            "Erythrocytes" = Hbb_bt,
                            "Macrophages" = Macrophages,
                            "Early spermatocytes" = early_spermatocytes)
for(i in 1:length(Cell_names_20180812)){
        SSCs <- RenameIdent.1(SSCs, new.ident.name = names(Cell_names_20180812)[i],
                                  cells.use = Cell_names_20180812[[i]])
}
TSNEPlot.1(object = SSCs, no.legend = F, do.label = T, pt.size = 1,
           do.return = TRUE, label.repel = T,label.size = 4,
           colors.use = singler.colors)+
        ggtitle("TSNE plot of all samples")+
        theme(text = element_text(size=20),							
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold"))

Cell_names_20180812 <- list2df(Cell_names_20180812)
write.csv(Cell_names_20180812,"./output/20180812/Cell_names_20180812.csv",na = "")
save(SSCs, file = "./data/SSCs_20180812.Rda")

#.....
Cell_names_20180812 = read.csv("./output/20180812/Cell_names_20180812.csv",header = 1, row.names = 1)
Cell_names_20180812 = apply(Cell_names_20180812, 2, function(x) x[x!= ""]) # df2list 
lapply(Cell_names_20180812,length)
names(Cell_names_20180812) = c("Spermatogonia","Early spermatocytes","Sertoli cells",
                               "Erythrocytes","Macrophages")

remove_Sertoli <- FeaturePlot(SSCs_Spermato, features.plot = "Rpl4",do.identify = T) # manually select 
#write.csv(remove_Sertoli,"./output/20180810/remove_Sertoli_20180809.csv",na = "")
remove_Sertoli = unique(c(as.vector(Cell_names_20180812$'Sertoli cells'),remove_Sertoli))
Cell_names_20180812$'Sertoli cells' = remove_Sertoli
#-----
remove_Sertoli <- read.csv("./output/20180810/remove_Sertoli_20180809.csv",row.names = 1)
remove_Sertoli = as.character(remove_Sertoli$x)
SSCs <- RenameIdent.1(SSCs,new.ident.name = "Sertoli cells", cells.use = remove_Sertoli)
table(SSCs@ident)

TSNEPlot.1(object = SSCs, no.legend = F, do.label = F, pt.size = 1,
         do.return = TRUE, label.size = 5,colors.use = singler.colors)+
        ggtitle("TSNE plot of easy way out")+
        theme(text = element_text(size=20),     							
              plot.title = element_text(hjust = 0.5))

"IS A DAY 6 GFRA1-POSITIVE SPERMATOGONIUM THE SAME AS A DAY 18 GFRA1-POSITIVE SPERMATOGONIUM?"
"IS A DAY 14 SYCP3-POSITIVE SPERMATOCYTE THE SAME AS AN ADULTSYCP3-POSITIVE SPERMATOCYTE?"
# Step 1. find DE genes among different cell types========
# 4.1.1 load data
lname1 = load(file = "./data/SSCs_20180812.Rda") ; lname1
dev_order <- c("Sertoli cells","Spermatogonia","Early spermatocytes","Spermatocytes",
               "Spermatids")
SSCs_Spermato = SubsetData(SSCs,ident.use = dev_order)
remove(SSCs);GC()
table(SSCs_Spermato@ident)
All_markers <- FindAllMarkers.UMI(SSCs_Spermato)
kable(All_markers[1:20,]) %>% kable_styling()
rownames(All_markers) = NULL
All_markers <- All_markers %>% select("gene", everything()) # Moving the last column to the start
write.csv(All_markers,"./output/20180812/All_markers.csv")
#----
All_markers <- read.csv("./output/20180812/All_markers.csv",header = T,row.names = 1)
Top_n = 15
top <-  All_markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC) %>%
        mutate(cell_type = factor(cluster, levels = dev_order)) %>%
        arrange(cell_type)
top %>% kable() %>% kable_styling()
DoHeatmap(object = SSCs_Spermato, genes.use = top$gene, 
          group.order = dev_order,
            slim.col.label = TRUE, remove.key = T,cex.row = 6,
            group.cex = 13, rotate.key = T,group.label.rot = T)+
        ggtitle(paste("Expression heatmap of top",Top_n,
                      "differential expression genes in each group"))+
        theme(text = element_text(size=20),     							
              plot.title = element_text(hjust = 0.5))
##################
Spermato_idents <- c("Spermatogonia","Early spermatocytes","Spermatocytes")
Spermato_markers = All_markers[All_markers$cluster %in% Spermato_idents,]
kable(Spermato_markers [1:20,]) %>% kable_styling()
write.csv(Spermato_markers,"./output/20180812/Spermato_markers.csv",na = "")
# Step 2. find DE genes among different day==============
Subset_RenameIdent_FindAllMarkers <- function(object,ident.use,merge =FALSE){
        # SubsetData
        sub_object <- SubsetData(object,ident.use = ident.use)
        # RenameIdent
        if(merge) {
                sub_object@meta.data$orig.ident = sub("PND18|PND25|PND30|Ad-depleteSp|Ad-Thy1",
                            "PND18_and_later",sub_object@meta.data$orig.ident)
        }
        ident.vector <- as.factor(sub_object@meta.data$orig.ident)
        names(x = ident.vector) <- names(sub_object@ident)
        sub_object@ident = ident.vector
        print(table(sub_object@meta.data$orig.ident))
        
        g <- TSNEPlot(object = sub_object,do.label = F, group.by = "ident", 
                 do.return = TRUE, no.legend = F, 
                 pt.size = 1,label.size = 8 )+
                ggtitle(ident.use)+
                theme(text = element_text(size=20),     							
                      plot.title = element_text(hjust = 0.5,size = 25, face = "bold")) 
        print(g)
        # FindAllMarkers
        sub_object_develop <- FindAllMarkers.UMI(sub_object)
        rownames(sub_object_develop) =NULL
        sub_object_develop <- sub_object_develop %>% select("gene", everything()) # Moving the last column to the start
        head(sub_object_develop)
        write.csv(sub_object_develop,file= paste0("./output/20180812/",
                                                  ident.use,"_develop.csv"))
        return(list(sub_object, sub_object_develop, g))
}
DoHeatmap.2 <- function(object, marker_df,Top_n = 10, Time_points,ident.use,
                        group.label.rot =T,cex.row = 6,
                        group.cex = 13){
        if(!all(marker_df$cluster %in% Time_points) || 
           !all(object@ident %in% Time_points)){
                # remove time points not int Time_points
                object <- SubsetData(object, ident.use = Time_points)
                marker_df <- marker_df[(marker_df$cluster %in% Time_points),] 
        }
        top <-  marker_df %>% group_by(cluster) %>% top_n(Top_n, avg_logFC) %>%
                mutate(time_points = factor(cluster, levels = Time_points)) %>%
                arrange(time_points)
        DoHeatmap(object = object, genes.use = top$gene, 
                  group.order = Time_points,
                  slim.col.label = TRUE, remove.key = T,cex.row = cex.row,
                  group.cex = group.cex, rotate.key = T,group.label.rot = group.label.rot)+
                ggtitle(paste("Expression heatmap of top",Top_n,
                              "differential expression genes in each time points",
                              ident.use))+
                theme(text = element_text(size=20),     							
                      plot.title = element_text(hjust = 0.5))
}

# Step 2. find DE genes among different day==============
# Spermatogonia-----
Spermatogonia_list <- Subset_RenameIdent_FindAllMarkers(SSCs_Spermato,
                                                        ident.use <- "Spermatogonia",
                                                        merge = T)
DoHeatmap.2(Spermatogonia_list[[1]],Spermatogonia_list[[2]],Top_n = 20, 
            c("PND06","PND14","PND18_and_later"),ident.use = "Spermatogonia",
            group.label.rot = F)
table(Spermatogonia_list[[1]]@ident) %>% t() %>% kable() %>% kable_styling()

# Early_spermatocytes_list --------
Early_spermatocytes_list <- Subset_RenameIdent_FindAllMarkers(SSCs_Spermato,
                                                        ident.use <- "Early spermatocytes",
                                                        merge = T)
DoHeatmap.2(Early_spermatocytes_list[[1]],Early_spermatocytes_list[[2]],Top_n = 30, 
            Time_points <- c("PND14","PND18_and_later"),
            ident.use = "Early spermatocytes",group.label.rot =F)
table(Early_spermatocytes_list[[1]]@ident) %>% t() %>% kable() %>% kable_styling()

#  Spermatocytes-----
Spermatocytes_list <- Subset_RenameIdent_FindAllMarkers(SSCs_Spermato,
                                                        ident.use <- "Spermatocytes")
DoHeatmap.2(Spermatocytes_list[[1]],Spermatocytes_list[[2]],Top_n = 10, 
            Time_points <- c("PND18","PND25","PND30","Ad-depleteSp","Ad-Thy1"),
            ident.use <- "Spermatocytes")
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
                        top_n(Top_n, avg_logFC) 
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


#2018/08/08
# DoHeatmap==================
# Spermatogonia
Spermatogonia_develop = read.csv("./output/Spermatogonia_dev_mark_new.csv",header = T)
Spermatogonia_develop = Spermatogonia_develop[,-1]
top10 <- Spermatogonia_develop %>% group_by(cluster) %>%
        top_n(10, avg_logFC)
top10 %>% kable() %>% kable_styling()
DoHeatmap(object = Spermatogonia, genes.use = top10$gene, 
          slim.col.label = TRUE, remove.key = TRUE,group.label.rot = T,
          title = "Spermatogonia development markers")

# Spermatocytes
Spermatocytes_develop = read.csv("./output/Spermatocytes_dev_mark_new.csv",header = T)
Spermatocytes_develop = Spermatocytes_develop[,-1]
top10 <- Spermatocytes_develop %>% group_by(cluster) %>%
        top_n(10, avg_logFC)
top10 %>% kable() %>% kable_styling()
DoHeatmap(object = Spermatocytes, genes.use = top10$gene, 
          slim.col.label = TRUE, remove.key = TRUE,group.label.rot = T,
          title = "Spermatocytes development markers")


# extrace cell names from DoHeatmap
All_markers <- read.csv("./output/20180812/All_markers.csv",header = T,row.names = 1)
top <-  All_markers %>% group_by(cluster) %>% top_n(20, avg_logFC)

DoHeatmap(object = Spermatocytes, genes.use = top$gene, cex.row = 6,
          group.cex = 12,
          slim.col.label = TRUE, remove.key = TRUE,group.label.rot = T,
          title = "Marker genes in Spermatocytes")

"1-a. Do Gfra1 and Id4 markers indicate the same clusters as SSCs?" # no
p1 <- SingleFeaturePlot.1(object = SSCs_Spermato, "Rbmxl2", threshold = 0.1)
p1
p2 <- SingleFeaturePlot.1(object = SSCs, "Id4", threshold = 0.8)
plot_grid(p1, p2)
SingleFeaturePlot.1(object = SSCs, "Hspe1", threshold = 0.8)
SingleFeaturePlot.1(object = SSCs, "Hmgn2", threshold = 0.8)
par(mfrow=c(1,2))
hist(SSCs@data["Gfra1",],freq = T,breaks = 30)
hist(SSCs@data["Id4",],freq = T,breaks = 30)

"1-a. What is the total cell count per library in'SSC clusters', 
and percent of SSC/all [germline] cells per library? "
kable(table(SSCs@ident, SSCs@meta.data$orig.ident)) %>%
        kable_styling()

prop.table(x = table(SSCs@ident, SSCs@meta.data$orig.ident),margin = 2) %>%
        kable  %>% kable_styling()

"1-b. What genes distinguish 'spermatogonia clusters' from non-spermatogonia clusters?
(expecting GFRa1 and ID4 to be high in this list by definition)
"
table(SSCs@ident)
idents <- as.data.frame(table(SSCs@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("non-Spermatogonia",
                     "non-Spermatogonia",
                     "non-Spermatogonia",
                     "non-Spermatogonia",
                     "non-Spermatogonia",
                     "non-Spermatogonia",
                     "non-Spermatogonia",
                     "Spermatogonia")
SSCs@ident <- plyr::mapvalues(x = SSCs@ident,
                              from = old.ident.ids,
                              to = new.cluster.ids)
TSNEPlot(object = SSCs,do.label = F, group.by = "ident", 
         do.return = TRUE, no.legend = T,
         pt.size = 1,label.size = 8 )+
        ggtitle("All samples")+
        theme(text = element_text(size=20),     							
              plot.title = element_text(hjust = 0.5))
#test
p1 <- SingleFeaturePlot.1(object = SSCs, "Sycp3", threshold = 0.1)
p2 <- SingleFeaturePlot.1(object = SSCs, "Sycp3", threshold = 2)
plot_grid(p1, p2)

p1 <- SingleFeaturePlot.1(object = SSCs, "Gfra1", threshold = 0.1)
p2 <- SingleFeaturePlot.1(object = SSCs, "Gfra1", threshold = 1)
plot_grid(p1, p2)


par(mfrow=c(1,2))
hist(SSCs@data["Rsph1",],freq = T,breaks = 30)
hist(SSCs@data["Tex101",],freq = T,breaks = 30)
"1-b. For these genes, what is the average expression within cells i'SSC clusters',
and broken down by library (mouse age)? We are most interested in genes with
different average expression levels among libraries.
Of top interest are genes that are highly expressed in PND06 cells 
in 'SSC clusters' and not as highly expressed (or as often detected) 
in cells in'SSC clusters' from older mice.
"
SSCs@meta.data$orig.ident <- gsub("zAd-","Ad-",SSCs@meta.data$orig.ident)



ident.vector <- as.character(x = SSCs@ident)
names(x = ident.vector) <- names(SSCs@ident)
ident.vector <- SSCs@meta.data$orig.ident

test <- factor(x = ident.vector, levels = SSCs@meta.data$orig.ident)
sample_markers <- FindAllMarkers.UMI(SSCs, only.pos = T)


uSSCs_As_only <- MouseGenes(SSCs,c("ID4","PAX7","BMI1","EOMES","GFRA1","FGFR3"))# As undifferentiated spermatogonia only
g <- lapply(uSSCs_As_only,function(x) SingleFeaturePlot.1(object = SSCs, feature = x, threshold = 0.5))

# 2018-08-17
lname1 = load(file = "./data/SSCs_20180812.Rda") ; lname1
SSCs@meta.data$orig.ident <- gsub("Ad-","zAd-",SSCs@meta.data$orig.ident)
VlnPlot(SSCs, features.plot = "nUMI",group.by = "ident", x.lab.rot = T)
VlnPlot(SSCs, features.plot = "nUMI",group.by = "orig.ident", x.lab.rot = T)
# 2018-08-18
Spermato_order <- c("Spermatogonia","Early spermatocytes","Spermatocytes")
Spermato = SubsetData(SSCs,ident.use = Spermato_order)
TSNEPlot(Spermato)
Spermato@meta.data$orig.ident <- gsub("Ad-","zAd-",Spermato@meta.data$orig.ident)
FeatureHeatmap_save(object = Spermato, gene.list = c("Wt1","Hnrnpa1"),save.folder = "Spermato")
FeatureHeatmap_save(object = Spermato, gene.list = grep("Mir",rownames(Spermato@data),value = T),
                    save.folder = "Spermato")
grep("Mir",rownames(Spermato@raw.data),value = T)
