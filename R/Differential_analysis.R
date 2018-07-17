########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
source("./R/Seurat_functions.R")

#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
# Rename ident
lnames = load(file = "./data/MCL_alignment.Rda")
lnames
table(MCL@ident)
idents <- as.data.frame(table(MCL@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("T cells",
                     "CD14 Monocytes",
                     "B cells",
                     "CD8 T cells",
                     "B cells",
                     "NK T cells",
                     "CD16 Monocytes",
                     "CD8 T cells",
                     "Dendritic Cells",
                     "Myeloid cells")

MCL@ident <- plyr::mapvalues(x = MCL@ident,
                             from = old.ident.ids,
                             to = new.cluster.ids)
MCL@meta.data$conditions <- sub("_",".",MCL@meta.data$conditions)
#===========answer to 6/1's email=====
Treg.markers.to.plot <- c("FOXP3","CTLA4","PDCD1","ENTPD1","CD38","ICOS",
                     "TNFSF9","TNFRSF9")
markers.to.plot <- HumanGenes(MCL,markers.to.plot,unique=T)
SplitDotPlotGG.1(object=MCL, genes.plot = rev(Treg.markers.to.plot),
                 cols.use = c("blue","red"), x.lab.rot = T, 
                 plot.legend = T, dot.scale = 8, do.return = T, 
                 grouping.var = "conditions")
#===========answer to 6/5's email=====
markers.to.plot <- c("CD79A","CD3G","CD8A","CD40","CXCR4")
markers.to.plot <- HumanGenes(MCL,markers.to.plot,unique=T)
SplitDotPlotGG.1(object=MCL, genes.plot = rev(markers.to.plot),
                 cols.use = c("blue","red"), x.lab.rot = T, 
                 plot.legend = T, dot.scale = 8, do.return = T, 
                 grouping.var = "conditions")
cellcycle.to.plot <- c("CCND1","CCND2","CCND3","CDK4","CDK6",
                       "PCNA","SOX11")
SplitDotPlotGG.1(object=MCL, genes.plot = rev(cellcycle.to.plot),
                 cols.use = c("blue","red"), x.lab.rot = T, 
                 plot.legend = T, dot.scale = 8, do.return = T, 
                 grouping.var = "conditions")

gde.all <- FindAllMarkersAcrossConditions(MCL)
write.csv(gde.all,"./doc/MCL_patient_vs_normal.csv")
#===========answer to 6/17's email=====
Featureplot <- function(x,object = object,...){
        p <- FeaturePlot(object = object, 
                         reduction.use = "tsne",
                         features.plot = x, min.cutoff = NA, 
                         cols.use = c("lightgrey","blue"), pt.size = 0.5,...)
        return(p)
}
cellcycle <- HumanGenes(MCL,c("CCND1","CCND2", "CCND3", "CDK4",
                              "CDK6","PCNA","SOX11"))
Treg <- HumanGenes(MCL,c("FOXP3","CTLA4","PDCD1","ENTPD1","CD38",
                         "ICOS","TNFSF9","TNFRSF9"))
markers.to.plot <- HumanGenes(MCL,c("CD79A","CD3G","CD8A","CD40","CXCR4"))
MCL.subsets <- SplitCells(MCL)
MCL.subsets[[3]]
MCL.patient <- MCL.subsets[[1]]
MCL.normal <- MCL.subsets[[2]]
# Featureplot
Featureplot(cellcycle,MCL.patient) # cellcycle
Featureplot(cellcycle,MCL.normal) # cellcycle
Featureplot(Treg,MCL.patient) # Treg
Featureplot(Treg,MCL.normal) # Treg
Featureplot(markers.to.plot,MCL.patient) 
Featureplot(markers.to.plot,MCL.normal) 

