library(Seurat)
library(SingleR)
library(dplyr)
library(tidyr)
library(kableExtra)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
#====== 3.1 load data  ==========================================
lname1 = load(file = "./data/SSCs_20180822.Rda");lname1
lname2 = load(file = "./output/singler_20180822.RData");lname2
SSCs@meta.data$orig.ident = gsub("Ad-","zAd-",SSCs@meta.data$orig.ident)
SSCs@meta.data$orig.ident = gsub("PDN18","PDN18pre",SSCs@meta.data$orig.ident)
# visualize
v1 <- VlnPlot(SSCs, "nUMI", group.by = "ident",x.lab.rot = T, y.log=T, do.return = T)
v2 <- VlnPlot(SSCs, "nUMI", group.by = "orig.ident",x.lab.rot = T, y.log=T, do.return = T)
plot_grid(v1,v2)
dev_order <- c("Spermatogonia","Spermatocytes","Spermatids")
SSCs_Spermato = SubsetData(SSCs,ident.use = dev_order)
table(SSCs_Spermato@ident)
SSCs_Spermato@ident <- factor(x = SSCs_Spermato@ident, levels = dev_order) # Relevel object@ident
#save(SSCs_Spermato, file = "./output/SSCs_Spermato.Rda")
#lnames = load(file = "./output/SSCs_Spermato.Rda") ;lnames

SSCs_Spermato <- ScaleData(object = SSCs_Spermato, model.use = "linear",
                           use.umi = T, vars.to.regress = c("percent.mito"),do.par = T)
v1 <- VlnPlot(SSCs_Spermato, "nUMI", group.by = "ident",x.lab.rot = T, 
              y.log=T, do.return = T, use.scaled =T)
v2 <- VlnPlot(SSCs_Spermato, "nUMI", group.by = "orig.ident",x.lab.rot = T,
              y.log=T, do.return = T, use.scaled =T)
plot_grid(v1,v2)

TSNEPlot.1(object = SSCs,do.label = T, group.by = "ident", 
           do.return = TRUE, no.legend = T,colors.use = singler.colors,
           pt.size = 1,label.size = 5,label.repel = T,force=25)+
        ggtitle("Supervised cell type labeling by GSE43717")+
        theme(text = element_text(size=20),							
              plot.title = element_text(hjust = 0.5, face = "bold"))

color.ident = which(levels(SSCs@ident) %in% dev_order)

TSNEPlot.1(object = SSCs_Spermato,do.label = T, group.by = "ident", 
                do.return = TRUE, no.legend = T,colors.use = singler.colors[rev(color.ident)],
                pt.size = 1,label.size = 5,label.repel = T,force=25)+
        ggtitle("Supervised cell type labeling by GSE43717")+
        theme(text = element_text(size=20),							
              plot.title = element_text(hjust = 0.5, face = "bold"))

#====== 3.2 Store cell names by identities ==========================================
#' convert cell.names ~ ident into a list of cell.name by ident
#' @param cell.names cell names in seurat.object@cell.names
#' @param ident cluster name or cell type in seurat.object@ident
#' @export labels_id A list of cell.names with ident as names
factor2list <- function(cell.names, ident){
        
        head(labels <- data.frame(cell.names = cell.names,
                                  ident = ident,
                                  stringsAsFactors = F))
        labels_tab <- labels %>% table() %>% as.data.frame.matrix()
        cell.names = rownames(labels_tab)
        labels_tab <- apply(labels_tab,2, as.logical) %>% as.data.frame()
        labels_id <- apply(labels_tab,2, function(x) cell.names[x])
        
        return(labels_id)
}

singler.labels <- singler$singler[[1]]$SingleR.single.main$labels
head(singler_labels_id <- factor2list(cell.names = rownames(singler.labels),
                                 ident = singler.labels[,1]),3)
save(singler_labels_id, file="./output/singler_labels.RData")
head(SSC_labels_id <- factor2list(cell.names = names(SSCs@ident),
                                  ident = as.character(SSCs@ident)),3)

#====== 3.3 unsupervised clustering =========================
SSCs <- FindClusters(object = SSCs, reduction.type = "pca", dims.use = 1:30, resolution = 1, 
                     k.param = 30 ,force.recalc = T,
                     save.SNN = TRUE, n.start = 100, nn.eps = 0, print.output = FALSE)
#SSCs@meta.data$orig.ident <- gsub("PND18pre","PND18",SSCs@meta.data$orig.ident)
TSNEPlot.1(object = SSCs,do.label = T, group.by = "ident", 
           do.return = TRUE, no.legend = T,colors.use = singler.colors,
           pt.size = 1,label.size = 4,label.repel = T,force=1)+
        ggtitle("Unsupervised clustering")+
        theme(text = element_text(size=20),							
              plot.title = element_text(hjust = 0.5,size=18, face = "bold"))
save(SSCs, file = "./data/SSCs_20180822.Rda")
#====== 3.4 SingleFeaturePlot =========================
Somatic_markers <- MouseGenes(SSCs,c("Col1a2","Acta2","Vcam1","Insl3","Laptm5",
                                     "Hbb-bt","Ptgds","Wt1"))
SingleFeaturePlot.1(object = SSCs, "Acta2", threshold = 1)
SingleFeaturePlot.1(object = SSCs, "Col1a2")
SingleFeaturePlot.1(object = SSCs, "Vcam1")
SingleFeaturePlot.1(object = SSCs, "Insl3")
SingleFeaturePlot.1(object = SSCs, "Laptm5")
SingleFeaturePlot.1(object = SSCs, "Hbb-bt", threshold = 1.5)
SingleFeaturePlot.1(object = SSCs, "Ptgds")
SingleFeaturePlot.1(object = SSCs, "Wt1", threshold = 1.5)
SingleFeaturePlot.1(object = SSCs, "Vim", threshold = 1.5)

GermCell_markers <- MouseGenes(SSCs,c("Gfra1","Zbtb16","Sall4","Dmrt1","Dazl","Kit","Id4",
                      "Sycp3","Mtl5","Nxt1","Shcbp1l","Aurka","Lyzl1",
                      "Acrv1","Hemgn","Txndc8","Tssk6","Oaz3","Prm2"))
SingleFeaturePlot.1(object = SSCs, "Gfra1", threshold = 0.8)
SingleFeaturePlot.1(object = SSCs, "Zbtb16", threshold = 1.5)
SingleFeaturePlot.1(object = SSCs, "Sall4", threshold = 1)
SingleFeaturePlot.1(object = SSCs, "Dmrt1", threshold = 0.1)
SingleFeaturePlot.1(object = SSCs, "Dazl", threshold = 2.0)
SingleFeaturePlot.1(object = SSCs, "Kit", threshold = 0.1)
SingleFeaturePlot.1(object = SSCs, "Id4", threshold = 0.8)
SingleFeaturePlot.1(object = SSCs, "Sycp3", threshold = 1.5)
SingleFeaturePlot.1(object = SSCs, "Mtl5", threshold = 1.0)
SingleFeaturePlot.1(object = SSCs, "Nxt1", threshold = 1.8)
SingleFeaturePlot.1(object = SSCs, "Shcbp1l", threshold = 1.5)
SingleFeaturePlot.1(object = SSCs, "Aurka", threshold = 1.5)
SingleFeaturePlot.1(object = SSCs, "Lyzl1", threshold = 1.5)
SingleFeaturePlot.1(object = SSCs, "Acrv1", threshold = 1.5)
SingleFeaturePlot.1(object = SSCs, "Hemgn", threshold = 1.5)
SingleFeaturePlot.1(object = SSCs, "Txndc8", threshold = 1.5)
SingleFeaturePlot.1(object = SSCs, "Tssk6", threshold = 1.5)
SingleFeaturePlot.1(object = SSCs, "Oaz3", threshold = 1.5)
SingleFeaturePlot.1(object = SSCs, "Prm2", threshold = 3.5)

#====== 3.5 DotPlot =========================
marker_order <- c(8,5,3,12,22,7,4,14,1,9,19,18,10,2,11,6,24,15,21,23,16,17,13,0,20)
SSCs@ident <- factor(x = SSCs@ident, levels = rev(marker_order)) # Relevel object@ident
markers.to.plot <- c(GermCell_markers,Somatic_markers)
DotPlot(SSCs, genes.plot = rev(markers.to.plot),
        cols.use = c("blue","red"), x.lab.rot = T, plot.legend = T,
        dot.scale = 8, do.return = T)

#====== 3.6 Rename ident =========================
lname1 = load(file = "./data/SSCs_20180822.Rda");lname1
table(SSCs@ident)
idents <- as.data.frame(table(SSCs@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("Sertoli cells",
                     "Spermatocytes",
                     "Spermatids",
                     "Early Spermatocytes",
                     "Spermatocytes",
                     "Spermatogonia",
                     "Spermatids",
                     "Spermatocytes",
                     "Spermatogonia",
                     "Spermatocytes",
                     "Spermatids",
                     "Spermatids",
                     "Early Spermatocytes",
                     "Sertoli cells",
                     "Spermatocytes",
                     "Smooth muscle",
                     "Endothelial & Hematopoietic cells",
                     "Endothelial & Hematopoietic cells",
                     "Spermatocytes",
                     "Spermatocytes",
                     "Sertoli cells",
                     "Endothelial & Hematopoietic cells",
                     "Early Spermatocytes",
                     "Endothelial & Hematopoietic cells",
                     "Smooth muscle")

SSCs@ident <- plyr::mapvalues(x = SSCs@ident,
                              from = old.ident.ids,
                              to = new.cluster.ids)
major_cells <- c("Sertoli cells","Spermatids","Spermatocytes","Early Spermatocytes",
                 "Spermatogonia","Smooth muscle","Endothelial & Hematopoietic cells")
SSCs@ident <- factor(x = SSCs@ident, levels = major_cells) # Relevel object@ident

TSNEPlot.1(object = SSCs,do.label = T, group.by = "ident", 
           do.return = TRUE, no.legend = T,colors.use = singler.colors[10:16],
           pt.size = 1,label.size = 5,label.repel = T)+
        ggtitle("Manually label cell types")+
        theme(text = element_text(size=20),							
              plot.title = element_text(hjust = 0.5, face = "bold"))

#====== 3.7 SubsetData =========================
dev_order <- c("Spermatogonia","Early Spermatocytes","Spermatocytes","Spermatids")
SSCs_Spermato = SubsetData(SSCs,ident.use = dev_order)
SSCs_Spermato@ident <- factor(x = SSCs_Spermato@ident, levels = rev(dev_order))
TSNEPlot.1(object = SSCs_Spermato,do.label = T, group.by = "ident", 
           do.return = TRUE, no.legend = T,colors.use = singler.colors[11:14],
           pt.size = 1,label.size = 5,label.repel = T)+
        ggtitle("Germ cells only")+
        theme(text = element_text(size=20),							
              plot.title = element_text(hjust = 0.5, face = "bold"))
# total cell number
SSCs_Spermato@meta.data$orig.ident = gsub("PND18pre","PND18",SSCs_Spermato@meta.data$orig.ident)
counts <- table(as.vector(SSCs_Spermato@ident), SSCs_Spermato@meta.data$orig.ident)
kable(counts) %>% kable_styling()

par(mfrow=c(1, 1), mar=c(5, 4, 4, 2))
barplot(counts, main="Total numbers of germ cells vs. time-points",
        xlab="Time points", ylab="Cell numbers",
        col=singler.colors[c(13,11,12,14)],
        legend = rownames(counts),
        args.legend = list(x = "topleft", bty = "n"))

# total cell percentage
prop_table <- prop.table(counts, margin = 2)
kable(prop_table) %>% kable_styling()

par(mfrow=c(1, 1), mar=c(5, 4, 4, 8))
barplot(prop_table, main="Percentage of germ cells vs. time-points",
        xlab="Time points", ylab="Cell percentage",
        col=singler.colors[c(13,11,12,14)],
        legend = rownames(counts),
        args.legend = list(x = "topright", bty = "n", inset=c(-0.15, 0)))