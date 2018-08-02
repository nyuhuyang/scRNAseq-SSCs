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
# Keep SSC related cells only
lnames = load(file = "./data/SSCs_label.Rda")
lnames
table(SSCs@ident)
#SSCs = SubsetData(SSCs,ident.remove = c("B cells","Dendritic cells","Endothelial cells",
#                                        "Erythrocytes","Granulocytes","Macrophages",
#                                        "Microglia","Monocytes","NK cells","T cells"))
SSCs = RenameIdent(SSCs, old.ident.name = "Peritubular/\nepithelial cells/\nwhite blood cells",
                   new.ident.name = "Peritubular, epithelial cells and white blood cells")
TSNEPlot(object = SSCs,do.label = T, group.by = "ident", 
           do.return = TRUE, no.legend = T,
           pt.size = 1,label.size = 8 )+
        ggtitle("Testis cells only")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5,size = 25, face = "bold")) 


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
Spermatogonia_markers <- FindMarkers.1(SSCs,ident.1="Spermatogonia",only.pos = T)
kable(Spermatogonia_markers[1:20,]) %>%
        kable_styling()
write.csv(Spermatogonia_markers,"./output/Spermatogonia_markers.csv")
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
levels(SSCs@ident) = SSCs@meta.data$orig.ident
head(SSCs@ident)

new.levels <- old.levels <- levels(x = SSCs@ident)
ident.vector <- as.character(x = SSCs@ident)
names(x = ident.vector) <- names(SSCs@ident)
ident.vector <- SSCs@meta.data$orig.ident

test <- factor(x = ident.vector, levels = SSCs@meta.data$orig.ident)
sample_markers <- FindAllMarkers.UMI(SSCs, only.pos = T)


uSSCs_As_only <- MouseGenes(SSCs,c("ID4","PAX7","BMI1","EOMES","GFRA1","FGFR3"))# As undifferentiated spermatogonia only
g <- lapply(uSSCs_As_only,function(x) SingleFeaturePlot.1(object = SSCs, feature = x, threshold = 0.5))
