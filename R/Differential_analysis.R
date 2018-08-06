########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(kableExtra)
source("../R/Seurat_functions.R")

#4.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 4.1.1 load data
lnames = load(file = "./data/SSCs_suplabel_GSE43717.Rda")
lnames
table(SSCs@ident)
table(SSCs@ident,SSCs@meta.data$orig.ident) %>% kable() %>%  kable_styling()
Spermatogonia = SubsetData(SSCs,ident.use = "Spermatogonia")
Spermatocytes = SubsetData(SSCs,ident.use = "Spermatocytes")
SSCs_pure = SubsetData(SSCs,ident.use = c("Spermatocytes","Spermatogonia"))

"IS A DAY 6 GFRA1-POSITIVE SPERMATOGONIUM THE SAME AS A DAY 18 GFRA1-POSITIVE SPERMATOGONIUM?"
"IS A DAY 14 SYCP3-POSITIVE SPERMATOCYTE THE SAME AS AN ADULTSYCP3-POSITIVE SPERMATOCYTE?"
# Step 1. find DE genes among different cell types========
All_markers <- FindAllMarkers.UMI(SSCs_pure,return.thresh = 0.1)
kable(All_markers[1:20,]) %>% kable_styling()
rownames(All_markers) = NULL
All_markers <- All_markers %>% select("gene", everything()) # Moving the last column to the start
write.csv(All_markers,"./output/All_markers_GSE43717.csv")

Spermatogonia_markers = All_markers[All_markers$cluster == "Spermatogonia",]
kable(Spermatogonia_markers[1:20,]) %>% kable_styling()
write.csv(Spermatogonia_markers,"./output/Spermatogonia_markers.csv")

Spermatocytes_markers = All_markers[All_markers$cluster == "Spermatocytes",]
kable(Spermatocytes_markers[1:20,]) %>% kable_styling()
write.csv(Spermatocytes_markers,"./output/Spermatocytes_markers.csv")

# Step 2. find DE genes among different day==============
# replace ident with orig.ident

Spermatogonia@meta.data$orig.ident = gsub("Ad-","zAd-",Spermatogonia@meta.data$orig.ident)
ident.vector <- as.factor(Spermatogonia@meta.data$orig.ident)
names(x = ident.vector) <- names(Spermatogonia@ident)
Spermatogonia@ident = ident.vector
Spermatogonia = SubsetData(Spermatogonia, ident.remove = "PND25")

Spermatocytes@meta.data$orig.ident = gsub("Ad-","zAd-",Spermatocytes@meta.data$orig.ident)
ident.vector <- as.factor(Spermatocytes@meta.data$orig.ident)
names(x = ident.vector) <- names(Spermatocytes@ident)
Spermatocytes@ident = ident.vector

TSNEPlot(object = SSCs_pure,do.label = F, group.by = "ident", 
         do.return = TRUE, no.legend = T,
         pt.size = 1,label.size = 8 )+
        ggtitle("Spermatogonia and Spermatocytes")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5,size = 25, face = "bold")) 

Spermatogonia_develop <- FindAllMarkers.UMI(Spermatogonia)
rownames(Spermatogonia_develop) =NULL
Spermatogonia_develop <- Spermatogonia_develop %>% select("gene", everything()) # Moving the last column to the start
kable(Spermatogonia_develop[1:20,]) %>% kable_styling()
write.csv(Spermatogonia_develop,"./output/Spermatogonia_develop_new.csv")

Spermatocytes_develop <- FindAllMarkers.UMI(Spermatocytes)
rownames(Spermatocytes_develop) =NULL
Spermatocytes_develop <- Spermatocytes_develop %>% select("gene", everything()) # Moving the last column to the start
kable(Spermatocytes_develop[1:20,]) %>% kable_styling()
write.csv(Spermatocytes_develop,"./output/Spermatocytes_develop.csv")

# Step 3. find overlap gene===========================
Spermatogonia_develop = read.csv("./output/Spermatogonia_develop.csv", header = TRUE)
kable(Spermatogonia_develop[1:20,]) %>% kable_styling()
Spermatocytes_develop = read.csv("./output/Spermatocytes_develop.csv", header = TRUE)
kable(Spermatocytes_develop[1:20,]) %>% kable_styling()
Spermatogonia_markers = read.csv("./output/Spermatogonia_markers.csv", header = TRUE)
kable(Spermatogonia_markers[1:20,]) %>% kable_styling()
Spermatocytes_markers = read.csv("./output/Spermatocytes_markers.csv", header = TRUE)
kable(Spermatocytes_markers[1:20,]) %>% kable_styling()

# Spermatogonia
table(Spermatogonia_develop$gene %in% Spermatogonia_markers$gene)
Spermatogonia_dev_mark <- Spermatogonia_develop[Spermatogonia_develop$gene %in%
                                                        Spermatogonia_markers$gene,]
kable(Spermatogonia_dev_mark[1:20,]) %>% kable_styling()
top3 <- Spermatogonia_dev_mark %>% group_by(cluster) %>% top_n(3, avg_logFC) 
top3 %>% kable() %>% kable_styling()
write.csv(Spermatogonia_dev_mark,"./output/Spermatogonia_dev_mark_new.csv")

# Spermatocytes
table(Spermatocytes_develop$gene %in% Spermatocytes_markers$gene)
Spermatocytes_dev_mark <- Spermatocytes_develop[Spermatocytes_develop$gene %in%
                                                        Spermatocytes_markers$gene,]
kable(Spermatocytes_dev_mark[1:20,]) %>% kable_styling()
top3 <- Spermatocytes_dev_mark %>% group_by(cluster) %>% top_n(3, avg_logFC) 
top3 %>% kable() %>% kable_styling()
write.csv(Spermatocytes_dev_mark,"./output/Spermatocytes_dev_mark_new.csv")

# FeatureHeatmap==========
SSCs_pure@meta.data$orig.ident = gsub("Ad-","zAd-",SSCs_pure@meta.data$orig.ident)
ident.vector <- as.factor(SSCs_pure@meta.data$orig.ident)
names(x = ident.vector) <- names(SSCs_pure@ident)
SSCs_pure@ident = ident.vector

Spermatocytes_GOI = read.delim("./doc/Spermatocytes_dev_mark_new_GOI.txt", header = T)
Spermatocytes_GOI$GOI.spermatocytes = as.character(Spermatocytes_GOI$GOI.spermatocytes)
for (i in 1:nrow(Spermatocytes_GOI)) {
        temp_plot = Customize_Seurat_FeatureHeatmap(SSCs_pure, Spermatocytes_GOI[i,1])+
                theme(text = element_text(size=15))
        ggsave(temp_plot, file=paste0("./output/tsne_matrix/Spermatocytes_GOI/",
                                      Spermatocytes_GOI[i,1],".jpeg"),
               width = 35.28, height = 24.69, units = "cm"
               )
}

Spermatogonia_GOI = read.delim("./doc/Spermatogonia_dev_mark_new_GOI.txt", header = T)
Spermatogonia_GOI$GOI.spermatogonia = as.character(Spermatogonia_GOI$GOI.spermatogonia)
for (i in 45:nrow(Spermatogonia_GOI)) {
        temp_plot = Customize_Seurat_FeatureHeatmap(SSCs_pure, Spermatogonia_GOI[i,1])+
                theme(text = element_text(size=15))
        ggsave(temp_plot, file=paste0("./output/tsne_matrix/Spermatogonia_GOI/",
                                      Spermatogonia_GOI[i,1],".jpeg"),
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
p1 <- SingleFeaturePlot.1(object = SSCs, "Gfra1", threshold = 0.5)
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
