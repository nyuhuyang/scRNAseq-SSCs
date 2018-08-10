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
lname2 = load(file = "./output/singler_Ref_GSE43717.RData") ; lname2
table(SSCs@ident)
table(SSCs@ident,SSCs@meta.data$orig.ident) %>% kable() %>%  kable_styling()


TSNEPlot(object = Spermato, no.legend = T, do.label = TRUE, pt.size = 1,
         do.return = TRUE, label.size = 5)+
        ggtitle("TSNE plot of Spermatocytes and Spermatogonia")+
        theme(text = element_text(size=20),							
              plot.title = element_text(hjust = 0.5))
# 08/07 easy way out...
remove_Sertoli <- FeaturePlot(Spermato, features.plot = "Rpl4",do.identify = T) # manually select 
write.csv(remove_Sertoli,"./output/20180810/remove_Sertoli_20180809.csv",na = "")
#-----
remove_Sertoli <- read.csv("./output/20180810/remove_Sertoli_20180809.csv",row.names = 1)
remove_Sertoli = as.character(remove_Sertoli$x)
SSCs <- RenameIdent.1(SSCs,new.ident.name = "Sertoli cells", cells.use = remove_Sertoli)
table(SSCs@ident)

TSNEPlot(object = SSCs, no.legend = F, do.label = F, pt.size = 1,
         do.return = TRUE, label.size = 5)+
        ggtitle("TSNE plot of easy way out")+
        theme(text = element_text(size=20),     							
              plot.title = element_text(hjust = 0.5))

"IS A DAY 6 GFRA1-POSITIVE SPERMATOGONIUM THE SAME AS A DAY 18 GFRA1-POSITIVE SPERMATOGONIUM?"
"IS A DAY 14 SYCP3-POSITIVE SPERMATOCYTE THE SAME AS AN ADULTSYCP3-POSITIVE SPERMATOCYTE?"
# Step 1. find DE genes among different cell types========
SSCs_Spermato = SubsetData(SSCs,ident.use = c("Sertoli cells","Spermatocytes",
                                              "Spermatogonia","Spermatids"))
All_markers <- FindAllMarkers.UMI(SSCs_Spermato)
kable(All_markers[1:20,]) %>% kable_styling()
rownames(All_markers) = NULL
All_markers <- All_markers %>% select("gene", everything()) # Moving the last column to the start
write.csv(All_markers,"./output/20180810/All_markers.csv")
top <-  All_markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
DoHeatmap(object = SSCs_Spermato, genes.use = top$gene, 
            slim.col.label = TRUE, remove.key = T,cex.row = 8,
            group.cex = 13, rotate.key = T,group.label.rot = F)+
        ggtitle("Expression heatmap of top 20 differential expression genes in each group")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle

Spermatogonia_markers = All_markers[All_markers$cluster == "Spermatogonia",]
kable(Spermatogonia_markers[1:20,]) %>% kable_styling()
write.csv(Spermatogonia_markers,"./output/20180810/Spermatogonia_markers.csv")

Spermatocytes_markers = All_markers[All_markers$cluster == "Spermatocytes",]
kable(Spermatocytes_markers[1:20,]) %>% kable_styling()
write.csv(Spermatocytes_markers,"./output/20180810/Spermatocytes_markers.csv")

# Step 2. find DE genes among different day==============
# replace ident with orig.ident
# merge 
Spermatogonia = SubsetData(SSCs,ident.use = c("Spermatogonia"))
table(Spermatogonia@meta.data$orig.ident) %>% t() %>% kable() %>% kable_styling()
Spermatogonia@meta.data$orig.ident = gsub("PND18|PND25|PND30|Ad-depleteSp|Ad-Thy1",
                                 "PND18_and_later",Spermatogonia@meta.data$orig.ident)
table(Spermatogonia@meta.data$orig.ident) %>% t() %>% kable() %>% kable_styling()
ident.vector <- as.factor(Spermatogonia@meta.data$orig.ident)
names(x = ident.vector) <- names(Spermatogonia@ident)
Spermatogonia@ident = ident.vector

Spermatocytes = SubsetData(SSCs,ident.use = c("Spermatocytes"))
table(Spermatocytes@meta.data$orig.ident)
Spermatocytes@meta.data$orig.ident = gsub("Ad-","zAd-",Spermatocytes@meta.data$orig.ident)
table(Spermatocytes@meta.data$orig.ident) %>% t() %>% kable() %>% kable_styling()
ident.vector <- as.factor(Spermatocytes@meta.data$orig.ident)
names(x = ident.vector) <- names(Spermatocytes@ident)
Spermatocytes@ident = ident.vector

TSNEPlot(object = Spermatocytes,do.label = F, group.by = "ident", 
         do.return = TRUE, no.legend = F,
         pt.size = 1,label.size = 8 )+
        ggtitle("Spermatocytes")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5,size = 25, face = "bold")) 

p1 <- SingleFeaturePlot.1(object = Spermatogonia, "Zbtb16", threshold = 0.1)
p2 <- SingleFeaturePlot.1(object = Spermatocytes, "Zbtb16", threshold = 0.8)
plot_grid(p1, p2)

MouseGenes(SSCs,"Lhcgr")
Spermatogonia_develop <- FindAllMarkers.UMI(Spermatogonia)
rownames(Spermatogonia_develop) =NULL
Spermatogonia_develop <- Spermatogonia_develop %>% select("gene", everything()) # Moving the last column to the start
kable(Spermatogonia_develop[1:20,]) %>% kable_styling()
top <- Spermatogonia_develop %>% group_by(cluster) %>% top_n(20, avg_logFC) 
write.csv(Spermatogonia_develop,"./output/20180810/Spermatogonia_develop.csv")

Spermatocytes_develop <- FindAllMarkers.UMI(Spermatocytes)
rownames(Spermatocytes_develop) =NULL
Spermatocytes_develop <- Spermatocytes_develop %>% select("gene", everything()) # Moving the last column to the start
kable(Spermatocytes_develop[1:20,]) %>% kable_styling()
top <- Spermatocytes_develop %>% group_by(cluster) %>% top_n(10, avg_logFC) 
write.csv(Spermatocytes_develop,"./output/20180810/Spermatocytes_develop.csv")

DoHeatmap(object = Spermatogonia, genes.use = top$gene, 
          slim.col.label = TRUE, remove.key = T,cex.row = 8,group.cex = 12,
          rotate.key = T,group.label.rot = T)+
        ggtitle("Expression heatmap of top 20 Spermatogonia development genes in each time point")+
        theme(text = element_text(size=10),							
              plot.title = element_text(hjust = 0.5))

# Step 3. find overlap gene===========================
Spermatogonia_develop = read.csv("./output/20180810/Spermatogonia_develop.csv", header = TRUE)
kable(Spermatogonia_develop[1:20,]) %>% kable_styling()
Spermatocytes_develop = read.csv("./output/20180810/Spermatocytes_develop.csv", header = TRUE)
kable(Spermatocytes_develop[1:20,]) %>% kable_styling()
All_markers = read.csv("./output/20180810/All_markers.csv", header = TRUE)
kable(All_markers[1:20,]) %>% kable_styling()

# Spermatogonia
Spermatogonia_markers <- All_markers[All_markers$cluster == "Spermatogonia",]
table(Spermatogonia_develop$gene %in% Spermatogonia_markers$gene)
Spermatogonia_dev_mark_genes <- intersect(Spermatogonia_develop$gene, Spermatogonia_markers$gene)
Spermatogonia_dev_mark <- Spermatogonia_develop[(Spermatogonia_develop$gene %in% 
                                                        Spermatogonia_dev_mark_genes),]
top <- Spermatogonia_dev_mark %>% group_by(cluster) %>% top_n(20, avg_logFC) 
write.csv(Spermatogonia_dev_mark,"./output/20180810/Spermatogonia_dev_mark.csv")
DoHeatmap.1(object = Spermatogonia, genes.use = top$gene, 
            slim.col.label = TRUE, remove.key = T,cex.row = 8,group.cex = 12,
            rotate.key = T,group.label.rot = T)+
        ggtitle("Expression heatmap of top 20 Spermatogonia development & marker genes in each time point")+
        theme(text = element_text(size=10),							
              plot.title = element_text(hjust = 0.5))

# Spermatocytes
Spermatocytes_markers <- All_markers[All_markers$cluster == "Spermatocytes",]
table(Spermatocytes_develop$gene %in% Spermatocytes_markers$gene)
Spermatocytes_dev_mark_genes <- intersect(Spermatocytes_develop$gene, Spermatocytes_markers$gene)
Spermatocytes_dev_mark <- Spermatocytes_develop[(Spermatocytes_develop$gene %in% 
                                                         Spermatocytes_dev_mark_genes),]
top <- Spermatocytes_dev_mark %>% group_by(cluster) %>% top_n(10, avg_logFC) 
write.csv(Spermatocytes_dev_mark,"./output/20180810/Spermatocytes_dev_mark.csv")
DoHeatmap(object = Spermatocytes, genes.use = top$gene, 
          slim.col.label = TRUE, remove.key = T,cex.row = 8,group.cex = 12,
          rotate.key = T,group.label.rot = T)+
        ggtitle("Expression heatmap of top 10 Spermatocytes development & markers genes in each time points")+
        theme(text = element_text(size=10),							
              plot.title = element_text(hjust = 0.5))

# FeatureHeatmap==========
Spermato = SubsetData(SSCs,ident.use = c("Spermatocytes","Spermatogonia"))
Spermatogonia_top100 <- Spermatogonia_dev_mark %>% group_by(cluster) %>% top_n(100, avg_logFC) 
Spermatogonia_top100 <- as.character(Spermatogonia_top100$gene)
for (gene in Spermatogonia_top100) {
        temp_plot = Customize_Seurat_FeatureHeatmap(Spermato, gene)+
                theme(text = element_text(size=15))
        ggsave(temp_plot, file=paste0("./output/20180810/tsne_matrix/Spermatogonia/",
                                      gene,".jpeg"),
               width = 35.28, height = 24.69, units = "cm"
        )
}

Spermatocytes_top100 <- Spermatocytes_dev_mark %>% group_by(cluster) %>% top_n(100, avg_logFC) 
Spermatocytes_top100 <- as.character(Spermatocytes_top100$gene)
for (gene in Spermatocytes_top100) {
        temp_plot = Customize_Seurat_FeatureHeatmap(Spermato, gene)+
                theme(text = element_text(size=15))
        ggsave(temp_plot, file=paste0("./output/20180810/tsne_matrix/Spermatocytes/",
                                      gene,".jpeg"),
               width = 35.28, height = 24.69, units = "cm"
        )
}
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


"1-a. Do Gfra1 and Id4 markers indicate the same clusters as SSCs?" # no
p1 <- SingleFeaturePlot.1(object = SSCs, "Kit", threshold = 0.1)
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
        theme(text = element_text(size=20),     #larger text including legend title							
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
