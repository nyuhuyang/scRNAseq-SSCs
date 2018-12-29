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
save(singler_labels_id, file="singler_labels.RData")
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

# Cell Type Marker gene database
# Stem cell=======
Embryonic_SCs <- MouseGenes(SSCs,c("ALPL","ALPP","ALPI","ALPPL2","TNFRSF8","TDGF1",
                                   "PODXL","NR6A1","POU5F1","B3GALNT1",
                                   "FUT1","FUT2","KIT","TERT","TERC","TEP1",
                                   "NANOG","SOX2","Dppa5","Klf4","zfp42","Spp1"))
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

Housekeeping <- MouseGenes(SSCs,c("Rnr2","Rpl4","Actb","Gnas","Tubb5","Kras","Calb1"))
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
Neuron2A <- MouseGenes(SSCs,c("Tmod2","Skil","Slc30a1","Erbb2ip","PCDHA@","Vgf","Gabrb3"))

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
Macrophages <- MouseGenes(SSCs,c("LYZ2","CD68","MARCO","EMR1","ITGB2"))
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
uSSCs <- unique(c(uSSCs_As_only,uSSCs_As_pr_al4,u_di_SSCs,other_SSCs))
u_di_SSCs <- MouseGenes(SSCs,c("NEUROG3","NANOS3","SOHLH1","MAGEA4","KIT","CD9"),unique = T) # un/differentiated
other_SSCs <- MouseGenes(SSCs,c("EXOSC10","SLC22A2","DMRT1","MKI67", "SSX2B","SSX3","SSX4",
                                "ITGA6","NANOG","CD9", "EpCAM","ADGRA3","GDNF","ITGB1","Ret","HLA-A",
                                "DDX4","DAZL","STRA8","CD24A","Nanos3","EGR3","FHL1","SOX3",
                                "TAF4B","Bcl6b","Numb","Lrp4","SOHLH2","CDH1","GNL3","UTF1","CST3","Vim",
                                "CENPB","NGN3","PRDX5","LRRC6","Rps27a"),unique = T) 
Sertoli <- MouseGenes(SSCs,c("TF","TFRC","TFR2"))
# ZBTB16(PLZF), SLC22A2(OCT2),MKI67(Ki67), ITGA6(α6-integrin,CD49f)
# THY1(CD90), ADGRA3(GPR125), ITGB1(β1-integrin, CD29),Ret(c-ret), VASA(DDX4)
# SERPINC1(EE2 antigen),Pou5f1(Oct4),EGR3,FHL1(RBM),GNL3(Nucleostemin)
Spermatocyte <- MouseGenes(SSCs,c("KHDRBS1","SPAG16","CATSPER2","SLC25A31","SPAG6"))
Spermatogonia <- MouseGenes(SSCs,c("CBL","ITGA6","ITGB1","CD9","PROM1",
                                   "CHEK2","DMRT1","DMRTB1","DSG2","ENO2",
                                   "ELAVL2","EPCAM","EXOSC10","FGFR3","FMR1",
                                   "GFRA1","ADGRA3","ID4","KIT","LIN28A","LIN28B",
                                   "MAGEA4","NANOG","NANOS2","NANOS3","PASD1",
                                   "PAX7","ZBTB16","POU2F2","POU5F1","SAGE1",
                                   "SALL4","SOX3","PHF13","SPOCD1","SSX1",
                                   "SSX2","SSX3","SSX4","SSX4B","THY1",
                                   "TRAPPC6A","TSPY1","TSPYL1","UCHL1","UTF1",
                                   "ZKSCAN2"))
defective_spermatozoa <- MouseGenes(SSCs,c("TXNDC8","TXNDC2","ALOX15","NME8")) #TXNDC8(SPTRX3), ALOX15(15-LOX)

# 
test <- MouseGenes(SSCs,c("ZBTB16")) #TXNDC8(SPTRX3), ALOX15(15-LOX)
Featureplot <- function(x,object = SSCs,...){
        p <- FeaturePlot(object = object, 
                         reduction.use = "tsne",
                         features.plot = x, min.cutoff = NA, 
                         cols.use = c("lightgrey","blue"), pt.size = 0.5,...)
        return(p)
}
Featureplot(test)
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
new.cluster.ids <- c("Sertoli cells 0",
                     "Spermatogonia 1",
                     "Spermatogonia/\nSpermatocytes 2",
                     "Spermatids 3",
                     "Unkonwn cells 4",
                     "Spermatocytes 5",
                     "Spermatocytes 6",
                     "Spermatocytes 7",
                     "Spermatids 8",
                     "Spermatocytes 9",
                     "Spermatids/\nRound spermatids 10",
                     "Spermatogonia 11",
                     "Spermatocytes 12",
                     "Spermatocytes 13",
                     "Spermatids 14",
                     "Spermatocytes 15",
                     "Spermatocytes 16",
                     "Spermatogonia 17",
                     "Spermatids 18",
                     "Round spermatids 19",
                     "Spermatocytes 20",
                     "Spermatocytes 21",
                     "Spermatids 22",
                     "Spermatocytes 23",
                     "Peritubular/\n epithelial cells/\n white blood cells 24",
                     "Spermatogonia 25",
                     "Spermatocytes 26",
                     "Spermatogonia 27",
                     "Spermatids 28",
                     "Peritubular/\n epithelial cells/\n white blood cells 29",
                     "Spermatogonia 30",
                     "Spermatocytes 31",
                     "Peritubular/\n epithelial cells/\n white blood cells 32",
                     "Peritubular/\n epithelial cells/\n white blood cells 33",
                     "Sertoli cells 34")

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
new.cluster.ids <- c("Sertoli cells",
                     "Spermatogonia",
                     "Spermatogonia/\nSpermatocytes",
                     "Spermatids",
                     "Unkonwn cells",
                     "Spermatocytes",
                     "Spermatocytes",
                     "Spermatocytes",
                     "Spermatids",
                     "Spermatocytes",
                     "Spermatids/\nRound spermatids",
                     "Spermatogonia",
                     "Spermatocytes",
                     "Spermatocytes",
                     "Spermatids",
                     "Spermatocytes",
                     "Spermatocytes",
                     "Spermatogonia",
                     "Spermatids",
                     "Round spermatids",
                     "Spermatocytes",
                     "Spermatocytes",
                     "Spermatids",
                     "Spermatocytes",
                     "Peritubular/\n epithelial cells/\n white blood cells",
                     "Spermatogonia",
                     "Spermatocytes",
                     "Spermatogonia",
                     "Spermatids",
                     "Peritubular/\n epithelial cells/\n white blood cells",
                     "Spermatogonia",
                     "Spermatocytes",
                     "Peritubular/\n epithelial cells/\n white blood cells",
                     "Peritubular/\n epithelial cells/\n white blood cells",
                     "Sertoli cells")


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
        ggtitle("Manually labeled cell types")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
        ggrepel::geom_label_repel(aes(label = ident))
SpermatogonialSC <- SubsetData(SSCs,ident.use = "Spermatogonial stem cells")
SpermatogonialSC_cells <- as.data.frame(table(SpermatogonialSC@meta.data$orig.ident))
All_cells <- as.data.frame(table(SSCs@meta.data$orig.ident))
SpermatogonialSC_cells$total <- All_cells$Freq
SpermatogonialSC_cells$percentage <- SpermatogonialSC_cells$Freq/SpermatogonialSC_cells$total
SpermatogonialSC_cells

# FeatureHeatmap
SSCs@meta.data$orig.ident <- gsub("Ad-","zAd-",SSCs@meta.data$orig.ident)
x <- FeatureHeatmap(object = SSCs, features.plot = c("Wt1","Gfra1","Zbtb16",
                                                    "Sycp3", "Txndc8"),
                    group.by = "orig.ident", sep.scale = T, pt.size = 0.5, 
                    cols.use = c("gray98", "red"), pch.use = 20, do.return = T)

x$data <- x$data[order(x$data$expression),]
customize_Seurat_FeatureHeatmap(x, alpha.use = 0.8,
                                scaled.expression.threshold = 0,
                                gradient.use = c("orangered", "red4"))
# histogram 
lnames = load(file = "./data/SSCs_alignment_8.Rda")
lnames
genes <- c("Txndc8","Spag16","Gfra1","Dmrt1")
CountsList <- list()
for(i in 1:length(genes)) CountsList[[i]] <- CountsbyIdent(object = SSCs,
                                                           subset.name = genes[i],
                                                           accept.low = 1)
CountsList[[length(genes)+1]] <- as.data.frame(table(SSCs@meta.data$orig.ident))
library(plyr)
CountsbyIdents <- join_all(CountsList, by='Var1', type='full')
CountsbyIdents[is.na(CountsbyIdents)] <- 0
rownames(CountsbyIdents) <- gsub("Ad-","zAd-",CountsbyIdents$Var1)
CountsbyIdents <- CountsbyIdents[,-1]

CountsbyIdents <- CountsbyIdents[sort(rownames(CountsbyIdents)),]
CountsbyIdents

percentbyIdents <- apply(CountsbyIdents,2,function(x){ x/CountsbyIdents$Freq } )
percentbyIdents <- data.frame(percentbyIdents)
percentbyIdents$samples <- rownames(percentbyIdents)
percentbyIdents <- percentbyIdents[,!(colnames(percentbyIdents) %in% "Freq")]
percentbyIdents
library(reshape2)
new_percentbyIdents <- melt(percentbyIdents,id=c("samples"))
colnames(new_percentbyIdents)[2] <- "Gene.name"
new_percentbyIdents
ggplot(new_percentbyIdents, aes(x = samples, y = value, color = Gene.name)) +
        geom_line(aes(group = Gene.name))+
        ggtitle("Cell percentage with gene expression in each sample")+#ggplot title
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5,size = 25, face = "bold"), #title in middle
              axis.text.x  = element_text(angle=30, vjust=0.5))+#rotate xlab
        guides(colour = guide_legend(override.aes = list(size=10)), #larger legend diagram 
               shape = guide_legend(override.aes = list(size=10))) #larger legend diagram 


######################
# Cell Type Marker gene database
lnames = load(file = "./data/SSCs_alignment_8.Rda")
lnames
SSCs_8 <- SSCs
remove(SSCs)
GC()
lnames = load(file = "./data/SSCs_alignment.Rda")
lnames
SingleFeaturePlot.1(object = SSCs, "Wt1", threshold = 1.5)
SingleFeaturePlot.1(object = SSCs, "Ar", threshold = 0.1)
SingleFeaturePlot.1(object = SSCs, "Acta2", threshold = 1)
SingleFeaturePlot.1(object = SSCs, "Vcam1")
SingleFeaturePlot.1(object = SSCs, "Laptm5")
SingleFeaturePlot.1(object = SSCs, "Vim", threshold = 0.1)
SingleFeaturePlot.1(object = SSCs, "Zbtb16", threshold = 0.5)
SingleFeaturePlot.1(object = SSCs, "Sycp3", threshold = 1.5)
SingleFeaturePlot.1(object = SSCs, "Acrv1", threshold = 1.5)
SingleFeaturePlot.1(object = SSCs, "Txndc8", threshold = 1.5)
SingleFeaturePlot.1(object = SSCs, "Hbb-bt", threshold = 1.5)
SingleFeaturePlot.1(object = SSCs, "Prm2", threshold = 1.5)
SingleFeaturePlot.1(object = SSCs, "Sall4", threshold = 1)
SingleFeaturePlot.1(object = SSCs, "Gfra1", threshold = 0.1)
SingleFeaturePlot.1(object = SSCs, "Dmrt1", threshold = 0.1)
SingleFeaturePlot.1(object = SSCs, "Id4", threshold = 1)
SingleFeaturePlot.1(object = SSCs, "Sall4", threshold = 0.1)
SingleFeaturePlot.1(object = SSCs, "Acta2", threshold = 1)
unknown <- FeaturePlot(object = SSCs, "Sall4",do.identify = T)
unknowns <- SubsetData(SSCs_8,cells.use = unknown)
table(unknowns@ident)
# rename
table(SSCs@ident)
idents <- as.data.frame(table(SSCs@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("Spermatogonia",
                     "Sertoli cells",
                     "Spermatogonia",
                     "Spermatids",
                     "Spermatocytes",
                     "Spermatocytes",
                     "Sertoli cells",
                     "Spermatocytes",
                     "Spermatocytes",
                     "Spermatocytes",
                     "Spermatids",
                     "Spermatocytes",
                     "Smooth Muscle",
                     "Spermatids",
                     "Peritubular/\n epithelial cells/\n white blood cells",
                     "Round spermatids",
                     "Spermatogonia",
                     "Spermatocytes",
                     "Spermatocytes",
                     "Spermatids",
                     "Spermatids",
                     "Peritubular/\n epithelial cells/\n white blood cells",
                     "Round spermatids",
                     "Spermatogonia",
                     "Spermatogonia",
                     "Peritubular/\n epithelial cells/\n white blood cells",
                     "Red blood cells",
                     "Sertoli cells",
                     "Spermatids",
                     "Spermatocytes",
                     "Peritubular/\n epithelial cells/\n white blood cells",
                     "Peritubular/\n epithelial cells/\n white blood cells",
                     "Spermatogonia",
                     "Peritubular/\n epithelial cells/\n white blood cells")

SSCs@ident <- plyr::mapvalues(x = SSCs@ident,
                              from = old.ident.ids,
                              to = new.cluster.ids)

TSNEPlot.1(object = SSCs,do.label = T, group.by = "ident", 
           do.return = TRUE, no.legend = T,
           pt.size = 1,label.size = 8 )+
        ggtitle("All samples")+
        theme(text = element_text(size=15),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5,size = 25, face = "bold")) 

# rename cells manually
Spermatogonia <- FeaturePlot(object = SSCs, "Zbtb16",do.identify = T,pt.size = 2)
Spermatocytes <- FeaturePlot(object = SSCs, "Sycp3",do.identify = T,pt.size = 2)
wbc <- FeaturePlot(object = SSCs, "Laptm5",do.identify = T,pt.size = 2)
Rbc <- FeaturePlot(object = SSCs, "Hbb-bt",do.identify = T,pt.size = 2)
Spermatids <- FeaturePlot(object = SSCs, "Acrv1",do.identify = T,pt.size = 2)

SSCs <- RenameIdent.1(SSCs, new.ident.name = "Spermatogonia",
                      cells.use = Spermatogonia) 
SSCs <- RenameIdent.1(SSCs, new.ident.name = "Spermatocytes",
                      cells.use = Spermatocytes)
SSCs <- RenameIdent.1(SSCs, new.ident.name = "Peritubular/\n epithelial cells/\n white blood cells",
                      cells.use = wbc) 
SSCs <- RenameIdent.1(SSCs, new.ident.name = "Red blood cells",
                      cells.use = Rbc) 
SSCs <- RenameIdent.1(SSCs, new.ident.name = "Spermatids",
                      cells.use = Spermatids) 
SSCs <- RenameIdent(SSCs, old.ident.name = "Peritubular/\n epithelial cells/\n white blood cells",
                    new.ident.name = "Peritubular/\nepithelial cells/\nwhite blood cells") 

TSNEPlot(object = SSCs,do.label = F, group.by = "ident", 
           do.return = TRUE, no.legend = T,
           pt.size = 1,label.size = 8 )+
        ggtitle("All samples")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5))

kable(table(SSCs@ident, SSCs@meta.data$orig.ident)) %>%
kable_styling()

save(SSCs, file = "./data/SSCs_label.Rda")

######################
# merge manual marker genes and FindAllMarkers

# for each cell type or cluster find which known markers are expressed in that cell type or cluster
# label each cell type or cluster using these markers
All_markers = read.csv("./output/All_markers_GSE43717.csv",header = T)
All_markers = All_markers[All_markers$avg_logFC>0,]
All_markers = All_markers[,-which(colnames(All_markers) == "X")]

All_markers[All_markers$gene %in% Spermatogonia,] %>% kable() %>% kable_styling()
SSCs_marker <- All_markers[All_markers$gene %in% Spermatogonia,"gene"]
SingleFeaturePlot.1(object = SSCs, "Itgb1", threshold = 1.5)
SingleFeaturePlot.1(object = SSCs, "Adgra3", threshold = 1)
SingleFeaturePlot.1(object = SSCs, "Prom1", threshold = 0.1)
SingleFeaturePlot.1(object = SSCs, "Dmrtb1", threshold = .5)
SingleFeaturePlot.1(object = SSCs, "Pou2f2", threshold = .5)
SingleFeaturePlot.1(object = SSCs, "Id4", threshold = 0.5)
SingleFeaturePlot.1(object = SSCs, "Trappc6a", threshold = 1.0)
SingleFeaturePlot.1(object = SSCs, "Cbl", threshold = 0.8)
SingleFeaturePlot.1(object = SSCs, "Uchl1", threshold = 1.5)
SingleFeaturePlot.1(object = SSCs, "Dmrt1", threshold = 1.)
SingleFeaturePlot.1(object = SSCs, "Phf13", threshold = 1.)
SingleFeaturePlot.1(object = SSCs, "Itga6", threshold = .5)
SingleFeaturePlot.1(object = SSCs, "Sall4", threshold = .5)
SingleFeaturePlot.1(object = SSCs, "Adgra3", threshold = 1.)
SingleFeaturePlot.1(object = SSCs, "Zbtb16", threshold = .1)
SingleFeaturePlot.1(object = SSCs, "Cd9", threshold = .5)
SingleFeaturePlot.1(object = SSCs, "Itgb1", threshold = 1.5)
SingleFeaturePlot.1(object = SSCs, "Fmr1", threshold = 1.5)
SingleFeaturePlot.1(object = SSCs, "Elavl2", threshold = .5)
SingleFeaturePlot.1(object = SSCs, "Nanos3", threshold = .3)
SingleFeaturePlot.1(object = SSCs, "Gfra1", threshold = .1)
SingleFeaturePlot.1(object = SSCs, "Epcam", threshold = 1.)
SingleFeaturePlot.1(object = SSCs, "Lyz2", threshold = 1.)
SingleFeaturePlot.1(object = SSCs, "Prm2", threshold = 3.0)
SingleFeaturePlot.1(object = SSCs, "Pou5f1", threshold = 0.1)
SingleFeaturePlot.1(object = SSCs, "Tspyl1", threshold = 0.1)
Spermatogonia


# spermatogonia development gene heatmap