library(Seurat)
library(dplyr)
source("./R/Seurat_functions.R")

#====== 2.1 identify phenotype for each cluster  ==========================================
lnames = load(file = "./data/SSCs_alignment_8.Rda")
lnames

Featureplot <- function(x,object = SSCs,...){
    p <- FeaturePlot(object = object, 
                     reduction.use = "tsne",
                     features.plot = x, min.cutoff = NA, 
                     cols.use = c("lightgrey","blue"), pt.size = 0.5,...)
    return(p)
}
# Cell Type Marker gene database
# Stem cell=======
Embryonic_SCs <- MouseGenes(SSCs,c("ALPL","ALPP","ALPI","ALPPL2","TNFRSF8","TDGF1",
                                   "PODXL","NR6A1","POU5F1","B3GALNT1",
                                   "FUT1","FUT2","KIT","TERT","TERC","TEP1",
                                   "NANOG","SOX2"))
Ectoderm <- MouseGenes(SSCs,c("NES","NCAM1","PAX6","VIM"))# VIM
Mesoderm <- MouseGenes(SSCs,c("BMP4","TBXT"))
Endoderm <- MouseGenes(SSCs,c("AFP","GATA4","HNF4A"))
Pluripotent_Stem_Cells <- unique(c(Embryonic_SCs,Ectoderm,Mesoderm,Endoderm))

mouse_lin <- MouseGenes(SSCs,c("CD2","CD3G","CD3D","CD4","CD5","CD8A","KLRB1","PTPRC","Ly76","Ly6g"))
human_lin <- MouseGenes(SSCs,c("CD3G","CD3D","CD14","FCGR3A","CD19","MS4A1","NCAM1"))
lineage_Cells <- unique(c(mouse_lin,human_lin))

MSC <- MouseGenes(SSCs,c("BMPR2","BMPR1A","BMPR1B","CD34","ATXN1","CD44","THY1",
                         mouse_lin,human_lin),unique = T)
HSC <- MouseGenes(SSCs,c("CD34","CD38","KIT","ATXN1","THY1",
                         mouse_lin,human_lin),unique = T)
Collagen <- MouseGenes(SSCs,c("COL2A1","COL1A1","COL1A2"),unique =T)


# Blood Vessel=====
Endothelium <- MouseGenes(SSCs,c("Cdh5","Pecam1","Flt1","Plvap","Kdr","ptprb",
                                 "Vwf","EMCN","Car4","VEGFA"))
Smooth_muscle_cells <- MouseGenes(SSCs,c("Acta2","MYH11"))
# Bone ==
Osteoblast <- MouseGenes(SSCs,c("ALPL","SPP1","IBSP","BGLAP"))

# Fat
Adipocytes <- MouseGenes(SSCs,c("TRIP4","ASCC1","SLC36A2","P2RX5","MYF5","UCP1"))
Fat_stem_cell <- MouseGenes(SSCs,c("FABP3","FABP1","FABP2","FABP7","FABP4",
                                   "FABP5","FABP12","FABP6","SLC27A1","SLC27A2",
                                   "SLC27A3","SLC27A4","SLC27A5","SLC27A6"))
# Gender
Y_chromosome <- MouseGenes(SSCs,c("UTY","DDX3Y"))
X_chromosome <- MouseGenes(SSCs,c("Xist"))

# Liver
Hepatocyte <- MouseGenes(SSCs,c("ALB","ITGB1"))

# Nervous System
Neuron

Epithelium <- MouseGenes(SSCs,c("Epcam","KRT19","KRT5",
                               "MUC1","SCGB3A2","SCGB1A1","SCGB3A1","SFTPB","FOXJ1","Rpe65",
                               "Rlbp1","Msln","Upk3b","Lrrn4"))
RPE <- MouseGenes(SSCs,c("Rpe65","Rlbp1"))
Fibroblast <- MouseGenes(SSCs,c("FGF1","FGF9","SFRP1"))

# Bone Marrow and Blood ===
Mesenchymal <- MouseGenes(SSCs,c("Pdgfrb","Vim","Has2","Dcn"))

#--Hematopoietic----
Hematopoietic <- MouseGenes(SSCs,c("PTPRC","LAPTM5","SRGN"))
#------Myeloid----
megakaryocytes <-  MouseGenes(SSCs,c("PPBP","GNG11"))
erythrocyte <-  MouseGenes(SSCs,c("HBA2","HBB"))
MastCells <- MouseGenes(SSCs,c("Cma1","Mcpt4","Tpsb2","Cpa3"))
Neutrophil <- MouseGenes(SSCs,c("ADAM8","MSMO1","FUT4","FCGR3A","CEACAM8"))
CD14_Monocytes <-  MouseGenes(SSCs,c("CD14","LYZ","S100A9","CCL2"))
CD16_Monocytes <- MouseGenes(SSCs,c("FCGR3A","MS4A7","VMO1"))
Monocytes <- unique(c(CD14_Monocytes,CD16_Monocytes))
Macrophages <- MouseGenes(SSCs,c("LYZ","CD68","MARCO","EMR1","ITGB2"))
DendriticCells <- MouseGenes(SSCs,c("Itgax","GPR183","CST3","HLA-DQA1","FCER1A","TSPAN13",
                                   "IL3RA","IGJ"))
Platelet <- MouseGenes(SSCs,c("VWF","CD36","ITGB3","GP9","GP6","GP1BA","SELP","PDGFRB")) # CD36 is not specific
Myeloid_all <-  MouseGenes(SSCs,c(megakaryocytes,erythrocyte,MastCells,
                             CD14_Monocytes,CD16_Monocytes,Macrophages,DendriticCells),unique=T)
#------Lymphoid----
Lymphoid <- MouseGenes(SSCs,c("Cd19","CD79A","MS4A1",
                             "GNLY","Ncr1","CCL5","KLRD1","NKG7"))
# T cell
T_Cell <- MouseGenes(SSCs,c("CD2","CD3G","CD3D","CD4","CD8A","IL2RA","FOXP3",
                           "IL7R","SELL","IL2RG","GIMAP5"))

Treg <- MouseGenes(SSCs,c("FOXP3","CD4","IL2RA","CTLA4","PDCD1","ENTPD1","CD38",
                         "ICOS","TNFSF9","TNFRSF9"))
CD4_Naive_T <- MouseGenes(SSCs,c("CD4","IL7R","GIMAP5","SELL","IL2RG"))
Natural_killer_T <- MouseGenes(SSCs,c("NKG7","CCL5","NCAM1","FCGR3A","Ncr1","KLRD1"))

T_Cell_all <- unique(c(T_Cell,Treg,CD4_Naive_T,Natural_killer_T))
# B cell
B_Cell <-MouseGenes(SSCs,c("CD19","MS4A1","CD79A","CD40","CD22","FCER2","HLA-DRB1",
                          "CXCR4","SOX11","CD5","PAX5","CD27","IL4R"))

B_StemCell <- MouseGenes(SSCs,c("SPN","CD20"))
Pre_Pro_B <- MouseGenes(SSCs,c("CD34","MME","CD38"))
Pro_B <- MouseGenes(SSCs,c("MME","CD19","SPN","CD38","CD24","IL7","IL3RA"))
Pre_B <- MouseGenes(SSCs,c("MME","CD19","MS4A1","CD24","CD38","IL7","IL3RA","IL4R"))
Immature_B <- MouseGenes(SSCs,c("MME","CD19","MS4A1","CR2","CD40","CD24","CD38","IL4R"))
Transitional_B <- MouseGenes(SSCs,c("CD19","MS4A1","CD5","CR2","CD24","CD38"))
Marginal_zone_B <- MouseGenes(SSCs,c("CD1C","CD19","MS4A1","CR2","CD27"))
Regulatory_B <- MouseGenes(SSCs,c("CD1D","CD5","CD19","CR2","CD24"))
Follicular_B <- MouseGenes(SSCs,c("CD19","MS4A1","CR2","CD22","FCER2","CD24",
                                 "HLA-DRB1","HLA-DQB1","HLA-DRA","HLA-DQA1"))
Activated_B <- MouseGenes(SSCs,c("CD27","CD19","MS4A1","IL2RA","TNFRSF8","CD69","CD80","CD86","FLT3"))
Germinal_center_B <- MouseGenes(SSCs,c("MME","CD19","MS4A1","FCER2","CD27","CD38","TNFRSF17"))
Plasma_blast <- MouseGenes(SSCs,c("CD19","CD38","CD27","TNFRSF17","HLA-DRB1"))
Plasma_cell_long_lived <- MouseGenes(SSCs,c("CXCR4","CD27","CD38","CD138","CD269"))
Memory_B <- MouseGenes(SSCs,c("CD19","MS4A1","CD40","CD27","CXCR4","CXCR5","ACKR3"))

B_Cell_all <- unique(c(B_Cell,B_StemCell,Pre_Pro_B,Pro_B,Pre_B,Immature_B,Transitional_B,
                       Marginal_zone_B,Regulatory_B,Follicular_B,Activated_B,Germinal_center_B,
                       Plasma_blast,Plasma_cell_long_lived,Memory_B))
Lymphoid_all <- unique(c(Lymphoid,T_Cell_all,B_Cell_all))
Hematopoietic_all <- unique(c(Hematopoietic,Myeloid_all,Lymphoid_all))

# other
Melanocytes <- MouseGenes(SSCs,c("Pmel","Mlana"))

Myelinating_Schwann_cells <- MouseGenes(SSCs,c("MBP","MPZ"))
Pericytes <- MouseGenes(SSCs,c("Pdgfrb","Cspg4","Anpep","Rgs5",
                              "Myh11","Mylk","Des","Vtn","Ifitm1"))
Stromal_fibroblasts <- MouseGenes(SSCs,c("DCN","COL6A1","TIMP3","PDGFRA"))
Neurons <- MouseGenes(SSCs,c("Ihh","Gli1", "Ptch1", "Hhip"))
CellCycle <- MouseGenes(SSCs,c("CCND1","CCND2","CCND3","CDK4","CDK6","PCNA","SOX11",
                               "RB1","E2F1","TK1","CCNA2","MKI67","CDK1"))

# undifferentiated spermatogonia
uSSCs_As_only <- MouseGenes(SSCs,c("ID4","PAX7","BMI1","EOMES","GFRA1","FGFR3"))# As undifferentiated spermatogonia only
uSSCs_As_pr_al4 <- MouseGenes(SSCs,c("NANOS2","UTF1","ZBTB16","SALL4","LIN28A",
                                     "FOXO1","DPPA4","UCHL1","UTF1"),unique = T)# expression  As, Apr and Aal4
uSSCs <- unique(c(uSSCs_As_only,uSSCs_As_pr_al4))
u_di_SSCs <- MouseGenes(SSCs,c("NEUROG3","NANOS3","SOHLH1","MAGEA4","KIT","CD9"),unique = T) # un/differentiated
other_SSCs <- MouseGenes(SSCs,c("EXOSC10","SLC22A2","DMRT1","MKI67", "SSX2B","SSX3","SSX4",
                                "ITGA6","NANOG","CD9", "EpCAM","ADGRA3","GDNF","ITGB1","Ret","HLA-A",
                                "DDX4","DAZL","STRA8","CD24A","Nanos3","EGR3","FHL1","SOX3",
                                "TAF4B","Bcl6b","Numb","Lrp4","SOHLH2","CDH1","GNL3","UTF1","CST3","Vim",
                                "CENPB","NGN3","PRDX5","LRRC6","Rps27a"),unique = T) 
# ZBTB16(PLZF), SLC22A2(OCT2),MKI67(Ki67), ITGA6(α6-integrin,CD49f)
# THY1(CD90), ADGRA3(GPR125), ITGB1(β1-integrin, CD29),Ret(c-ret), VASA(DDX4)
# SERPINC1(EE2 antigen),Pou5f1(Oct4),EGR3,FHL1(RBM),GNL3(Nucleostemin)
Spermatocyte <- MouseGenes(SSCs,c("KHDRBS1","SPAG16","CATSPER2","SLC25A31","SPAG6"))
spermatozoids <- MouseGenes(SSCs,c("THY1","STRC","DEL15Q15.3","DAZ1","DAZ2","DAZ3","DAZ3"))
defective_spermatozoa <- MouseGenes(SSCs,c("TXNDC8","TXNDC2","ALOX15","NME8")) #TXNDC8(SPTRX3), ALOX15(15-LOX)

#  test FindAllMarkers=============
AllMarkers <- FindAllMarkers.UMI(SSCs, logfc.threshold = 0.25,
                                 min.pct = 0.25,only.pos = T)
write.csv(AllMarkers,"./output/AllMarkers.csv")
AllMarkers <- readr::read_csv("output/AllMarkers.csv")

MarkerList <- list("Embryonic Stem Cells"=Embryonic_SCs,
                   "Ectoderm"=Ectoderm,
                   "Mesoderm"=Mesoderm,
                   "Endoderm"=Endoderm,
                   "lineage Cells"=lineage_Cells,
                   "Mesencyhmal stem cell"=MSC,
                   "Hematopoietic stem cell"=HSC,
                   "Collagen"=Collagen,
                   "Endothelium"=Endothelium,
                   "Smooth muscle cells"=Smooth_muscle_cells,
                   "Osteoblast"=Osteoblast,
                   "Adipocytes"=Adipocytes,
                   "Fat stem cell"=Fat_stem_cell,
                   #"Hepatocyte"=Hepatocyte,
                   "Epithelium"=Epithelium,
                   "Fibroblast"=Fibroblast,
                   "Mesenchymal cells"=Mesenchymal,
                   "Hematopoietic cells"=Hematopoietic,
                   #"Myeloid"=Myeloid_all,
                   #"Lymphoid"=Lymphoid,
                   #"T_Cell"=T_Cell_all,
                   #"B_Cell"=B_Cell_all,
                   "Pericytes"=Pericytes,
                   #"CellCycle"=CellCycle,
                   "undifferentiated spermatogonia"=uSSCs,
                   "un_differentiated spermatogonia"=u_di_SSCs,
                   "SSCs related"=other_SSCs,
                   "Spermatocyte"=Spermatocyte,
                   "defective spermatozoa"=defective_spermatozoa
)
df_Markers <- Marker2Types(MarkerList)
df_Markers <- inner_join(df_Markers,AllMarkers, by="gene")
df_Markers <- df_Markers[,!(colnames(df_Markers) %in% c("X1","pct.1","pct.2"))]
df_Markers <- df_Markers %>% select("cluster", everything())
df_Markers <- df_Markers[order(df_Markers$cluster),]
write.csv(df_Markers,"./output/Markers_CellTypes.csv")
# Featureplot=======================================
# create list of markers by "Cell_Type"
List <- Types2Markers(df_Markers)
names(List)
Featureplot("Thy1")
Featureplot(c(List$defective_spermatozoa,"Rps27a"))
FeaturePlot(SSCs, "Txndc8",do.hover=T)
Abnormal_sperm <- df_Markers$Cell_Type=="defective spermatozoa" & df_Markers$gene=="Txndc8"
Abnormal_sperm <- unique(df_Markers[Abnormal_sperm,"cluster"])
Abnormal_sperm
Abnormal_sperm <- SubsetData(SSCs,ident.use = Abnormal_sperm)
Abnormal_sperm_cells <- as.data.frame(table(Abnormal_sperm@meta.data$orig.ident))
All_cells <- as.data.frame(table(SSCs@meta.data$orig.ident))
Abnormal_sperm_cells$total <- All_cells$Freq
Abnormal_sperm_cells$percentage <- Abnormal_sperm_cells$Freq/Abnormal_sperm_cells$total
Abnormal_sperm_cells

Featureplot(c(List$Spermatocyte))
Featureplot(c("Khdrbs1","Txndc8","Spag16","Spag6"))
Normal_sperm <- df_Markers$Cell_Type=="Spermatocyte" & df_Markers$gene=="Spag16"
Normal_sperm <- unique(df_Markers[Normal_sperm,"cluster"])
Normal_sperm <- Normal_sperm[!(Normal_sperm %in% c(19,26,28))] # remove 19,26,28
Normal_sperm <- unique(c(Normal_sperm,9,12))
Normal_sperm
Normal_sperm_cells <- SubsetData(SSCs,ident.use = Normal_sperm)
Normal_sperm_cells <- as.data.frame(table(Normal_sperm_cells@meta.data$orig.ident))
All_cells <- as.data.frame(table(SSCs@meta.data$orig.ident))
Normal_sperm_cells$total <- All_cells$Freq
Normal_sperm_cells$percentage <- Normal_sperm_cells$Freq/Normal_sperm_cells$total
Normal_sperm_cells

WBS <- df_Markers$Cell_Type=="Hematopoietic cells"
WBS <- unique(df_Markers[WBS,"cluster"])
WBS <- WBS[!(WBS %in% Abnormal_sperm)]
WBS
Collagens <- df_Markers$Cell_Type=="Collagen" & df_Markers$gene=="Col1a1"
Collagens <- unique(df_Markers[Collagens,"cluster"])
Collagens <- Collagens[!(Collagens %in% WBS)]
Collagens

WBS <- df_Markers$Cell_Type=="Hematopoietic cells"
WBS <- unique(df_Markers[WBS,"cluster"])
WBS <- WBS[!(WBS %in% Abnormal_sperm)]
not_SCs <- sort(unique(c(Abnormal_sperm,Normal_sperm,Collagen,WBS)))
SCs <- SubsetData(SSCs,ident.remove = not_SCs)

p1 <- SingleFeaturePlot.1(object = SSCs, feature="Txndc8")
p2 <- SingleFeaturePlot.1(object = SCs, feature="Txndc8")
p3 <- SingleFeaturePlot.1(object = SSCs, feature="Spag16")
p4 <- SingleFeaturePlot.1(object = SCs, feature="Spag16")
p5 <- SingleFeaturePlot.1(object = SSCs, feature="Col1a2")
p6 <- SingleFeaturePlot.1(object = SCs, feature="Col1a2")
p7 <- SingleFeaturePlot.1(object = SSCs, feature="Laptm5")
p8 <- SingleFeaturePlot.1(object = SCs, feature="Laptm5")
plot_grid(p1, p2, p3, p4,p5,p6,p7,p8,ncol = 2)

p1 <- SingleFeaturePlot.1(object = SSCs, feature="Gfra1")
p2 <- SingleFeaturePlot.1(object = SCs, feature="Gfra1")
p3 <- SingleFeaturePlot.1(object = SSCs, feature="Dmrt1")
p4 <- SingleFeaturePlot.1(object = SCs, feature="Dmrt1")
plot_grid(p1, p2, p3, p4,ncol = 2)


Featureplot(uSSCs_As_only,object = not_SCs)
Featureplot(c(List$Adipocytes,List$Fat_stem_cell))

Featureplot(c(List$Embryonic_Stem_Cells,List$Ectoderm,List$Endoderm)[c(1,5:7)])
Featureplot(unique(c(List$Collagen, List$Smooth_muscle_cells,List$Pericytes)))
Featureplot(c(List$Epithelium,List$Fibroblast))
Featureplot(c(List$Spermatocyte))
Featureplot("Thy1")
Featureplot(List$Hematopoietic_cells)
Featureplot(List$Hematopoietic_stem_cell)
Featureplot(c(List$undifferentiated_spermatogonia,"Pax7","Eomes"))
Featureplot(List$un_differentiated_spermatogonia)
others <- List$SSCs_related[!(List$SSCs_related %in% c(List$undifferentiated_spermatogonia,
                                   List$un_differentiated_spermatogonia))]
Featureplot(others[1:9])
Featureplot(others[10:18])
Featureplot(spermatozoids)

#====== 2.2 dot Plots ==========================================
markers.to.plot <- c(Hematopoietic[c(2,1,3)],Spermatocyte[c(2,3)],
                     Collagen[-1],defective_spermatozoa[-3],
                     Spermatocyte[c(5,2,3)],"Dmrt1","Pou5f1",uSSCs[c(5,9:14)],
                     "Vim","Itgb1","Gata4","Adgra3","Fhl1", "Itga6","NANOG","SOX2")
markers.to.plot <- MouseGenes(SSCs,markers.to.plot,unique=T)
# Rename ident
table(SSCs@ident)
idents <- as.data.frame(table(SSCs@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("Other stem cells 0",
                     "Spermatogonial stem cells 1",
                     "Spermatogonial stem cells 2",
                     "Defective sperm 3",
                     "Collagen 4",
                     "Spermatocyte 5",
                     "Spermatocyte 6",
                     "Spermatocyte 7",
                     "Defective sperm 8",
                     "Spermatocyte 9",
                     "Defective sperm 10",
                     "Spermatogonial stem cells 11",
                     "Spermatocyte 12",
                     "Spermatocyte 13",
                     "Defective sperm 14",
                     "Spermatocyte 15",
                     "Spermatocyte 16",
                     "Other stem cells 17",
                     "Defective sperm 18",
                     "Defective sperm 19",
                     "unknown 20",
                     "Spermatocyte 21",
                     "Defective sperm 22",
                     "Spermatocyte 23",
                     "White blood cells 24",
                     "Other stem cells 25",
                     "Defective sperm 26",
                     "Other stem cells 27",
                     "Defective sperm 28",
                     "White blood cells 29",
                     "Spermatogonial stem cells 30",
                     "Defective sperm 31",
                     "unknown 32",
                     "White blood cells 33",
                     "unknown 34")

SSCs@ident <- plyr::mapvalues(x = SSCs@ident,
                                    from = old.ident.ids,
                                    to = new.cluster.ids)
DotPlot(SSCs, genes.plot = rev(markers.to.plot),
        cols.use = c("blue","red"), x.lab.rot = T, plot.legend = T,
        dot.scale = 8, do.return = T)
# SSCs <- RenameIdentBack(SSCs)


lnames = load(file = "./data/SSCs_alignment_8.Rda")
lnames
table(SSCs@ident)
idents <- as.data.frame(table(SSCs@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("Other stem cells",
                     "Spermatogonial stem cells",
                     "Spermatogonial stem cells",
                     "Defective sperm",
                     "Collagen",
                     "Spermatocyte",
                     "Spermatocyte",
                     "Spermatocyte",
                     "Defective sperm",
                     "Spermatocyte",
                     "Defective sperm",
                     "Spermatogonial stem cells",
                     "Spermatocyte",
                     "Spermatocyte",
                     "Defective sperm",
                     "Spermatocyte",
                     "Spermatocyte",
                     "Other stem cells",
                     "Defective sperm",
                     "Defective sperm",
                     "unknown",
                     "Spermatocyte",
                     "Defective sperm",
                     "Spermatocyte",
                     "White blood cells",
                     "Other stem cells",
                     "Defective sperm",
                     "Other stem cells",
                     "Defective sperm",
                     "White blood cells",
                     "Spermatogonial stem cells",
                     "Defective sperm",
                     "unknown",
                     "White blood cells",
                     "unknown")

SSCs@ident <- plyr::mapvalues(x = SSCs@ident,
                                    from = old.ident.ids,
                                    to = new.cluster.ids)
DotPlot(SSCs, genes.plot = rev(markers.to.plot),
        cols.use = c("blue","red"), x.lab.rot = T, plot.legend = T,
        dot.scale = 8, do.return = T)
p1 <- TSNEPlot(SSCs, do.return = T, pt.size = 1, group.by = "orig.ident")
p2 <- TSNEPlot(SSCs, do.return = T, pt.size = 1, group.by = "ident")
#png('./output/TSNESplot_alignment.png')
plot_grid(p1, p2)

TSNEPlot(object = SSCs, no.legend = F, do.label = F,
         do.return = TRUE, label.size = 5)+
        ggtitle("TSNE plot of major cell types")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle

SpermatogonialSC <- SubsetData(SSCs,ident.use = "Spermatogonial stem cells")
SpermatogonialSC_cells <- as.data.frame(table(SpermatogonialSC@meta.data$orig.ident))
All_cells <- as.data.frame(table(SSCs@meta.data$orig.ident))
SpermatogonialSC_cells$total <- All_cells$Freq
SpermatogonialSC_cells$percentage <- SpermatogonialSC_cells$Freq/SpermatogonialSC_cells$total
SpermatogonialSC_cells

x <- FeatureHeatmap(object = SSCs, features.plot = c("Txndc8","Spag16"),
                    group.by = "orig.ident", sep.scale = T, pt.size = 1, 
                    cols.use = c("lightgrey","blue"), pch.use = 20, do.return = T)

x$data <- x$data[order(x$data$expression),]
x
#c("Pou5f1",uSSCs[c(5,9:10)],
#====== 2.4 Compare cell type changes across conditions  ==========================================
# the two patients profiled have very different composition
# Compare clusters for each dataset
SplitTSNEPlot(object = SSCs, split.by = "conditions")
SplitTSNEPlot(object = SSCs, split.by = "orig.ident")
#  test FindAllMarkers
AllMarkers_cells <- FindAllMarkers.UMI(SSCs)
AllMarkers_cells <- AllMarkers_cells[AllMarkers_cells$avg_logFC > 0,]
write.csv(AllMarkers_cells,"./output/AllMarkers_cells.csv")