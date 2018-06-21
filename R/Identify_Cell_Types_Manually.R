library(Seurat)
library(dplyr)
source("./R/Seurat_functions.R")

#====== 2.1 identify phenotype for each cluster  ==========================================
lnames = load(file = "./data/MCL_alignment.Rda")
lnames

Featureplot <- function(x,object = MCL,...){
    p <- FeaturePlot(object = object, 
                     reduction.use = "tsne",
                     features.plot = x, min.cutoff = NA, 
                     cols.use = c("lightgrey","blue"), pt.size = 0.5,...)
    return(p)
}
#------
Adipocytes <- HumanGenes(MCL,c("SLC36A2","P2RX5","MYF5","UCP1","TRIP4","ASCC1"))
Endothelium <- HumanGenes(MCL,c("Cdh5","Pecam1","Flt1","Plvap","Kdr","ptprb",
                                  "Vwf","EMCN","Car4"))
Epithelium <- HumanGenes(MCL,c("Epcam","KRT19","KRT5",
                                 "MUC1","SCGB3A2","SCGB1A1","SCGB3A1","SFTPB","FOXJ1","Rpe65",
                                 "Rlbp1","Msln","Upk3b","Lrrn4"))
RPE <- HumanGenes(MCL,c("Rpe65","Rlbp1"))
Fibroblast <- HumanGenes(MCL,c("FGF1","FGF9","SFRP1"))
#--Hematopoietic----
Hematopoietic <- HumanGenes(MCL,c("PTPRC","LAPTM5","SRGN"))
#------Myeloid----
megakaryocytes <-  HumanGenes(MCL,c("PPBP","GNG11"))
erythrocyte <-  HumanGenes(MCL,c("HBA2","HBB"))
MastCells <- HumanGenes(MCL,c("Cma1","Mcpt4","Tpsb2","Cpa3"))
Neutrophil <- HumanGenes(MCL,c("ADAM8","MSMO1","FUT4","FCGR3A","CEACAM8"))
CD14_Monocytes <-  HumanGenes(MCL,c("CD14","LYZ","S100A9","CCL2"))
CD16_Monocytes <- HumanGenes(MCL,c("FCGR3A","MS4A7","VMO1"))
Macrophages <- HumanGenes(MCL,c("LYZ","CD68","MARCO","Emr1"))
DendriticCells <- HumanGenes(MCL,c("Itgax","GPR183","CST3","HLA-DQA1","FCER1A","TSPAN13",
                                     "IL3RA","IGJ"))
Myeloid <-  HumanGenes(MCL,c(megakaryocytes,erythrocyte,MastCells,
                               Monocytes,Macrophages,DendriticCells))
#------Lymphoid----
Lymphoid <- HumanGenes(MCL,c("Cd19","CD79A","MS4A1",
                               "GNLY","Ncr1","CCL5","KLRD1","NKG7"))
# T cell
T_Cell <- HumanGenes(MCL,c("CD2","CD3G","CD3D","CD8A","IL2RA","FOXP3"))
Treg <- HumanGenes(MCL,c("FOXP3","CTLA4","PDCD1","ENTPD1","CD38",
                         "ICOS","TNFSF9","TNFRSF9"))
CD4_Naive_T <- HumanGenes(MCL,c("CD4","IL7R","GIMAP5","SELL","IL2RG"))
NK <- HumanGenes(MCL,c("NKG7","CCL5"))
Regulatory_T <-  HumanGenes(MCL,c("CD4","IL2RA","FOXP3"))
Natural_killer_T <-  HumanGenes(MCL,c("NCAM1","FCGR3A"))
# B cell
B_Cell <-HumanGenes(MCL,c("CD19","MS4A1","CD79A","CD40","CD22","FCER2",
                          "HLA-DRB1","CXCR4"))
B_StemCell <- HumanGenes(MCL,c("SPN","CD20"))
Pre_Pro_B <- HumanGenes(MCL,c("CD34","MME","CD38"))
Pro_B <- HumanGenes(MCL,c("MME","CD19","SPN","CD38","CD24","IL7","IL3RA"))
Pre_B <- HumanGenes(MCL,c("MME","CD19","MS4A1","CD24","CD38","IL7","IL3RA","IL4R"))
Immature_B <- HumanGenes(MCL,c("MME","CD19","MS4A1","CR2","CD40","CD24","CD38","IL4R"))
Transitional_B <- HumanGenes(MCL,c("CD19","MS4A1","CD5","CR2","CD24","CD38"))
Marginal_zone_B <- HumanGenes(MCL,c("CD1C","CD19","MS4A1","CR2","CD27"))
Regulatory_B <- HumanGenes(MCL,c("CD1D","CD5","CD19","CR2","CD24"))
Follicular_B <- HumanGenes(MCL,c("CD19","MS4A1","CR2","CD22","FCER2","CD24",
                                   "HLA-DRB1","HLA-DQB1","HLA-DRA","HLA-DQA1"))
Activated_B <- HumanGenes(MCL,c("CD27","CD19","MS4A1","IL2RA","TNFRSF8","CD69","CD80","CD86","FLT3"))
Germinal_center_B <- HumanGenes(MCL,c("MME","CD19","MS4A1","FCER2","CD27","CD38","TNFRSF17"))
Plasma_blast <- HumanGenes(MCL,c("CD19","CD38","CD27","TNFRSF17","HLA-DRB1"))
Plasma_cell_long_lived <- HumanGenes(MCL,c("CXCR4","CD27","CD38","CD138","CD269"))
Memory_B <- HumanGenes(MCL,c("CD19","MS4A1","CD40","CD27","CXCR4","CXCR5","ACKR3"))


Melanocytes <- HumanGenes(MCL,c("Pmel","Mlana"))
Mesenchymal <- HumanGenes(MCL,c("Pdgfrb","Vim","Has2","Dcn"))
Myelinating_Schwann_cells <- HumanGenes(MCL,c("MBP","MPZ"))
Pericytes <- HumanGenes(MCL,c("Pdgfrb","Cspg4","Anpep","Rgs5",
                                "Myh11","Mylk","Des","Vtn","Ifitm1"))
Smooth_muscle_cells <- HumanGenes(MCL,c("Acta2","Myh11"))
Stem_cell <- HumanGenes(MCL,c("POU5F1","FUT4","CD34","PROM1","ABCG2","Runx1","ATXN1",
                                "Nes","NCAM","NGFR"))
Stromal_fibroblasts <- HumanGenes(MCL,c("DCN","COL6A1","TIMP3","PDGFRA"))
Neurons <- HumanGenes(MCL,c("Ihh","Gli1", "Ptch1", "Hhip"))
cellcycle <- HumanGenes(MCL,c("CCND1","CCND2", "CCND3", "CDK4",
                              "CDK6","PCNA","SOX11"))
# Featureplot
Featureplot(Adipocytes) # Adipocytes
Featureplot(Endothelium) # Endothelial Cells
Featureplot(Epithelium) # Epithelium
Featureplot(c(RPE,Melanocytes,Myelinating_Schwann_cells)) # RPE, Melanocytes, Myelinating Schwann cells
Featureplot(Fibroblast) # Fibroblasts

#==================
Featureplot(Hematopoietic) # Hematopoietic cells
Featureplot(Myeloid[c(5:7,9:11,13:14,18:20)]) # Myeloid cells

Featureplot(erythrocyte)
Featureplot(MastCells)
Featureplot(Neutrophil)
Featureplot(c(CD14_Monocytes,CD16_Monocytes))
Featureplot(Macrophages)
Featureplot(DendriticCells)
#=====================
Featureplot(Lymphoid) # Lymphoid cells
Featureplot(NK)
# T cell
Featureplot(c(T_Cell[1:6],Natural_killer_T))
Featureplot(Treg)
Featureplot(CD4_Naive_T)
Featureplot(c(Regulatory_T,Natural_killer_T))
# B cell
Featureplot(B_Cell)
Featureplot(unique(c(B_StemCell,
                     Pre_Pro_B,
                     Pro_B,
                     Pre_B)))
Featureplot(unique(c(Immature_B,
                     Transitional_B)))
Featureplot(Marginal_zone_B)
Featureplot(unique(c(Regulatory_B,
                     Activated_B)))
Featureplot(Follicular_B)
Featureplot(Germinal_center_B)
Featureplot(unique(c(Plasma_blast,
                     Plasma_cell_long_lived,
                     Memory_B)))

Featureplot(Mesenchymal) # Mesenchymal cells
Featureplot(Pericytes) # Pericytes
Featureplot(Smooth_muscle_cells)
Featureplot(Stem_cell)
Featureplot(Stromal_fibroblasts)
Featureplot(Neurons)

markers.to.plot <- c(CD14_Monocytes[c(1:3)],DendriticCells[c(3,5)],
                     Macrophages[2],Natural_killer_T[2],
                     CD16_Monocytes[1:3],T_Cell[1:4],B_Cell[1:6])
markers.to.plot <- c("CD14","LYZ","S100A9","CST3","CD68","FCER1A","FCGR3A","MS4A7","VMO1",
                     "CD2","CD3G","CD3D","CD8A","CD19","MS4A1","CD79A","CD40","CD22",
                     "FCER2")
markers.to.plot <- HumanGenes(MCL,markers.to.plot,unique=T)
DotPlot(MCL, genes.plot = rev(markers.to.plot),
        cols.use = c("blue","red"), x.lab.rot = T, plot.legend = F,
        dot.scale = 8, do.return = T)
#  test FindAllMarkers
AllMarkers <- FindAllMarkers.UMI(MCL)
AllMarkers <- AllMarkers[AllMarkers$avg_logFC > 0,]
write.csv(AllMarkers,"./output/AllMarkers_clusters.csv")
# Rename ident
table(MCL@ident)
idents <- as.data.frame(table(MCL@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("T cells 0",
                     "CD14 Monocytes 1",
                     "B cells 2",
                     "CD8 T cells 3",
                     "B cells 4",
                     "NK T cells 5",
                     "CD16 Monocytes 6",
                     "CD8 T cells 7",
                     "Dendritic Cells 8",
                     "Myeloid cells 9")

MCL@ident <- plyr::mapvalues(x = MCL@ident,
                                    from = old.ident.ids,
                                    to = new.cluster.ids)
DotPlot(MCL, genes.plot = rev(markers.to.plot),
        cols.use = c("blue","red"), x.lab.rot = T, plot.legend = F,
        dot.scale = 8, do.return = T)
# MCL <- RenameIdentBack(MCL)

#====== 2.2 dot Plots ==========================================
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
DotPlot(MCL.normal, genes.plot = rev(markers.to.plot),
        cols.use = c("blue","red"), x.lab.rot = T, plot.legend = T,
        dot.scale = 8, do.return = T)
MCL@meta.data$conditions <- sub("_",".",MCL@meta.data$conditions)
SplitDotPlotGG(object=MCL, genes.plot = rev(markers.to.plot),
               cols.use = c("blue","red"), x.lab.rot = T, 
               plot.legend = T, dot.scale = 8, do.return = T, 
               grouping.var = "conditions")
FindMarkers()
TSNEPlot(object = MCL, no.legend = TRUE, do.label = TRUE,
         do.return = TRUE, label.size = 5)+
        ggtitle("TSNE plot of major cell types")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle

freq_table <- prop.table(x = table(MCL@ident, MCL@meta.data[, "conditions"]),
margin = 2)
barplot(height = freq_table)
freq_table
table(MCL@meta.data[, "conditions"])

#====== 2.4 Compare cell type changes across conditions  ==========================================
# the two patients profiled have very different composition
# Compare clusters for each dataset
SplitTSNEPlot(object = MCL)
#  test FindAllMarkers
AllMarkers_cells <- FindAllMarkers.UMI(MCL)
AllMarkers_cells <- AllMarkers_cells[AllMarkers_cells$avg_logFC > 0,]
write.csv(AllMarkers_cells,"./output/AllMarkers_cells.csv")

#===========answer to 6/1's email
MCL.subsets <- SplitCells(MCL)
MCL.subsets[[3]]
MCL.patient <- MCL.subsets[[1]]
MCL.normal <- MCL.subsets[[2]]
g1 <- Featureplot(Treg[1:2],object=MCL.patient,do.return=T)
g2 <- Featureplot(Treg[1:2],object=MCL.normal,do.return=T)
g <- list()
g[[1]] <- g1[[1]];g[[2]] <- g1[[2]];g[[3]] <- g2[[1]];g[[4]] <- g2[[2]];
do.call(plot_grid, g)
