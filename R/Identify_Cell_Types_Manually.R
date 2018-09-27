library(Seurat)
library(SingleR)
library(dplyr)
library(tidyr)
library(kableExtra)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
#====== 3.1 load data  ==========================================
lname1 = load(file = "./data/SSCs_20180825.Rda");lname1
#SSCs@meta.data$orig.ident = gsub("PDN18","PDN18pre",SSCs@meta.data$orig.ident)
TSNEPlot.1(object = SSCs,do.label = T, group.by = "ident", 
           do.return = TRUE, no.legend = T, colors.use = singler.colors,
           pt.size = 1,label.size = 5,label.repel = T,force=1)+
        ggtitle("TSNEPlot of all clusters")+
        theme(text = element_text(size=20),							
              plot.title = element_text(hjust = 0.5,size=18, face = "bold"))
#====== 3.2 Identify cell types by SingleFeaturePlot =========================
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
SingleFeaturePlot.1(object = SSCs, "Col1a2", threshold = 1.5)

GermCell_markers <- MouseGenes(SSCs,c("Gfra1","Zbtb16","Sall4","Dmrt1","Dazl","Kit","Cdca8",
                                      "Id4","Sycp3","Mtl5","Nxt1","Shcbp1l","Aurka","Lyzl1",
                                      "Acrv1","Hemgn","Txndc8","Tssk6","Oaz3","Prm2"))
SingleFeaturePlot.1(object = SSCs, "Gfra1", threshold = 0.1)
SingleFeaturePlot.1(object = SSCs, "Zbtb16", threshold = 1.5)
SingleFeaturePlot.1(object = SSCs, "Sall4", threshold = 1)
SingleFeaturePlot.1(object = SSCs, "Dmrt1", threshold = 0.1)
SingleFeaturePlot.1(object = SSCs, "Dazl", threshold = 2.0)
SingleFeaturePlot.1(object = SSCs, "Kit", threshold = 0.1)
SingleFeaturePlot.1(object = SSCs, "Id4", threshold = 0.8)
SingleFeaturePlot.1(object = SSCs, "Cdca8", threshold = 1.0)
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

#====== 3.3 DotPlot =========================
marker_order <- c(8,6,3,9,0,15,2,7,18,17,10,11,5,4,12,20,21,14,16,13,1,19)
SSCs@ident <- factor(x = SSCs@ident, levels = rev(marker_order))
DotPlot(SSCs, genes.plot = rev(markers.to.plot),
        cols.use = c("blue","red"), x.lab.rot = T, plot.legend = T,
        dot.scale = 8, do.return = T)

#====== 3.4 Rename ident =========================
lname1 = load(file = "./data/SSCs_20180822.Rda");lname1
table(SSCs@ident)
idents <- as.data.frame(table(SSCs@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("Spermatocytes",
                     "Sertoli cells",
                     "Spermatocytes",
                     "Spermatogonia",
                     "Spermatids",
                     "Spermatids",
                     "Spermatogonia",
                     "Spermatocytes",
                     "Spermatogonia",
                     "Early Spermatocytes",
                     "Round Spermatids",
                     "Spermatids",
                     "Smooth muscle",
                     "Sertoli cells",
                     "Endothelial & Hematopoietic cells",
                     "Spermatocytes",
                     "Endothelial & Hematopoietic cells",
                     "Spermatocytes",
                     "Spermatocytes",
                     "Sertoli cells",
                     "Endothelial & Hematopoietic cells",
                     "Endothelial & Hematopoietic cells")

SSCs@ident <- plyr::mapvalues(x = SSCs@ident,
                              from = old.ident.ids,
                              to = new.cluster.ids)
major_cells <- c("Endothelial & Hematopoietic cells","Sertoli cells","Smooth muscle",
                 "Spermatogonia","Early Spermatocytes","Spermatocytes",
                 "Round Spermatids","Spermatids")
SSCs@ident <- factor(x = SSCs@ident, levels = major_cells) # Relevel object@ident

TSNEPlot.1(object = SSCs,do.label = T, group.by = "ident", 
           do.return = TRUE, no.legend = T, colors.use = singler.colors,
           pt.size = 1,label.size = 5,label.repel = T)+
        ggtitle("Unsupervised clustering and discovery of cell types by marker genes")+
        theme(text = element_text(size=20),							
              plot.title = element_text(hjust = 0.5,size=18, face = "bold"))
SSCs <- StashIdent(object = SSCs, save.name = "Cell.Types")
save(SSCs, file = "./data/SSCs_20180926.Rda")
#====== 3.5 hand pick cells ident and Rename ident again =========================

# 3.5.1 hand pick cells ident===
#' convert cell.names ~ ident into a list of cell.name by each ident
#' @param cell.names cell names in seurat.object_at_cell.names
#' @param ident cluster name or cell type in seurat.object_at_ident
#' @export labels_id A list of cell.names with ident as names
#' @example 
# SSC_labels_id <- factor2list(cell.names = names(SSCs@ident),
#                              ident = as.vector(SSCs@ident))
#' SSC_labels_id %>% lapply(head)
# SSC_labels_id %>% lapply(length) %>% unlist # compare with table(SSCs@ident)
factor2list <- function(cell.names, ident){
        labels <- data.frame(cell.names = cell.names,
                                  ident = ident,
                                  stringsAsFactors = F)
        labels_tab <- labels %>% table() %>% as.data.frame.matrix()
        cell.names = rownames(labels_tab)
        labels_tab <- apply(labels_tab,2, as.logical) %>% as.data.frame()
        labels_id <- apply(labels_tab,2, function(x) cell.names[x])
        
        return(labels_id)
}


SSC_labels_id <- factor2list(cell.names = names(SSCs@ident),
                             ident = as.character(SSCs@ident))
SSC_labels_id %>% lapply(length) %>% unlist # compare with table(SSCs@ident)

# hand pick 
SSC_labels_id$`Smooth muscle` = FeaturePlot(SSCs,"Acta2", do.identify = T,
                                            cols.use = c("grey", "red"))
SSC_labels_id$`Endothelial & Hematopoietic cells` = FeaturePlot(SSCs,"Hbb-bt", do.identify = T,
                                            cols.use = c("grey", "red"))
SSC_labels_id$Spermatids = FeaturePlot(SSCs,"Txndc8", do.identify = T,
                                       cols.use = c("grey", "red"))

SSC_labels_id$Spermatogonia = c(SSC_labels_id$Spermatogonia,
                                FeaturePlot(SSCs,"Gfra1", do.identify = T,
                                       cols.use = c("grey", "red")))
SSC_labels_id$`Early Spermatocytes` = c(SSC_labels_id$`Early Spermatocytes`,
                                FeaturePlot(SSCs,"Kit", do.identify = T,
                                            cols.use = c("grey", "red")))

SSC_labels_id$Spermatocytes = c(SSC_labels_id$Spermatocytes,
                                FeaturePlot(SSCs,"Nxt1", do.identify = T,
                                            cols.use = c("grey", "red")))


# 3.5.2 Rename ident again ====

#' rename ident for certain cells only
#' One Perk, will double check the duplicated cell names
#' @param object Seurat object
#' @param new.ident.name one character
#' @param cells.use a vector of cell names
#' @export object Seurat object with new ident
#' @example 
# SSCs@ident = RenameIdent.2(object = SSCs, 
#                            new.ident.name = "Smooth muscle",
#                            cells.use = SSC_labels_id[["Smooth muscle"]])
RenameIdent.2 <- function(object, new.ident.name, cells.use){
        
        old.ident = object@ident
        old.in.new = old.ident[names(old.ident) %in% cells.use] %>%
                droplevels()
        old.not.in.new = old.ident[!(names(old.ident) %in% cells.use)] %>%
                droplevels()
        old.levels = levels(old.in.new)
        wrong.types = old.levels[!(old.levels %in% new.ident.name)]
        old.in.new[old.in.new %in% wrong.types] = new.ident.name
        new.ident.ids <- unlist(list(old.in.new, old.not.in.new)) # concatenate two factors
        object@ident = factor(x = new.ident.ids, levels = levels(object@ident))
        
        return(object)
}

new.types <- c("Spermatogonia")
for(new.type in new.types) {
        SSCs = RenameIdent.2(object = SSCs,
                               new.ident.name = new.type,
                               cells.use = SSC_labels_id[[new.type]])
        }

TSNEPlot.1(object = SSCs,do.label = T, group.by = "ident", 
           do.return = TRUE, no.legend = T, colors.use = singler.colors,
           pt.size = 1,label.size = 5,label.repel = T,force=1)+
        ggtitle("Unsupervised clustering and discovery of cell types by marker genes")+
        theme(text = element_text(size=20),							
              plot.title = element_text(hjust = 0.5,size=18, face = "bold"))

SSC_labels_id <- factor2list(cell.names = names(SSCs@ident),
                              ident = as.vector(SSCs@ident))
# back to hand pick if need more
saveRDS(SSC_labels_id, file = "./output/SSC_labels_20180825.Rda")
save(SSCs, file = "./data/SSCs_20180825.Rda")
SSC_labels_id = readRDS(file = "./output/SSC_labels_20180825.Rda")

# 20180925
ident.use <- data.frame("ident.use" = as.vector(SSCs@ident),
                        row.names = names(SSCs@ident))
SSCs <- AddMetaData(SSCs, ident.use,"ident.use")
gene.use <- rownames(SSCs@scale.data);length(gene.use)
SSCs@data = SSCs@data[gene.use,]
save(SSCs, file = "./data/SSCs_20180925.Rda")

#====== 3.6 SubsetData SSCs_Spermato =========================
# total cell number
dev_order <- c("Spermatogonia","Early Spermatocytes","Spermatocytes",
                 "Round Spermatids","Spermatids")
SSCs_Spermato <- SubsetData(SSCs, ident.use =  dev_order)
SSCs_Spermato@meta.data$orig.ident = sub("PND18pre","PND18",SSCs_Spermato@meta.data$orig.ident)
SSCs_Spermato@meta.data$orig.ident = sub("Ad-","zAd-",SSCs_Spermato@meta.data$orig.ident)
#SSCs_Spermato@meta.data$orig.ident = sub("zAd-","Ad-",SSCs_Spermato@meta.data$orig.ident)

TSNEPlot.1(object = SSCs_Spermato,do.label = T, group.by = "ident", 
           do.return = TRUE, no.legend = T, colors.use = singler.colors[4:8],
           pt.size = 1,label.size = 5,label.repel = T,force=1)+
        ggtitle("Germ cells only")+
        theme(text = element_text(size=20),							
              plot.title = element_text(hjust = 0.5,size=18, face = "bold"))

counts <- table(as.vector(SSCs_Spermato@ident), SSCs_Spermato@meta.data$orig.ident)
kable(counts) %>% kable_styling()

par(mfrow=c(1, 1), mar=c(5, 4, 4, 2))
barplot(counts, main="Total numbers of germ cells vs. time-points",
        xlab="Time points", ylab="Cell numbers",
        col=singler.colors[c(5,7,8,6,4)],
        legend = rownames(counts),
        args.legend = list(x = "topleft", bty = "n"))

# total cell percentage
prop_table <- prop.table(counts, margin = 2)
kable(prop_table) %>% kable_styling()

par(mfrow=c(1, 1), mar=c(5, 4, 4, 9))
barplot(prop_table, main="Percentage of germ cells vs. time-points",
        xlab="Time points", ylab="Cell percentage",
        col=singler.colors[c(5,7,8,6,4)],
        legend = rownames(counts),
        args.legend = list(x = "topright", bty = "n", inset=c(-0.2, 0)))
