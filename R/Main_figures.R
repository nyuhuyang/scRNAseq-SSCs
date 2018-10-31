####This code demonstrates how the main and supplementary figures were generated. 
library(Seurat)
library(SingleR)
library(dplyr)
library(kableExtra)
library(MAST)
library(plyr)
library(scales)
library(gridExtra)
library(grid)
library(cowplot)
source("../R/Seurat_functions.R")

path <- paste("./output",gsub("-","",Sys.Date()),sep = "/")
if(!dir.exists(path)) dir.create(path, recursive = T)

#################################################################
#
# Figure 1 B) Germ cell composition per library (bar chart)
#
#################################################################

lname1 = load(file = "./data/SSCs_20180926.Rda");lname1
SSCs@meta.data$orig.ident = gsub("PND18pre","PND18",SSCs@meta.data$orig.ident)
table(SSCs@meta.data$orig.ident)
table(SSCs@ident)
major_cells <- c("Spermatogonia","Early Spermatocytes","Spermatocytes",
                 "Round Spermatids","Spermatids","Sertoli cells",
                 "Endothelial & Hematopoietic cells", "Smooth muscle")
SSCs_spermato <- SubsetData(SSCs,ident.use = major_cells[1:5])
counts <- table(as.vector(SSCs_spermato@ident), SSCs_spermato@meta.data$orig.ident)
kable(counts) %>% kable_styling()
major_cells = major_cells[1:5]
time_series <- c("PND06","PND14","PND18",#"PND18pre",
                 "PND25","PND30","Ad-depleteSp","Ad-Thy1")
counts = counts[rev(major_cells),time_series]

jpeg(paste0(path,"/1_B_bar_chart.jpeg"), units="in", width=10, height=7,
     res=600)
par(mfrow=c(1, 1), mar=c(5, 4, 4, 2))
barplot(counts, main="Total numbers of germ cells vs. time-points",
        xlab="Time points", ylab="Cell numbers",
        col=singler.colors[8:4],
        legend = rownames(counts),
        args.legend = list(x = "topleft", bty = "n"))
dev.off()

#################################################################
#
# Figure 2 A) Composite tSNE showing all libraries
#
#################################################################
lname1 = load(file = "./data/SSCs_20180926.Rda");lname1
SSCs@meta.data$orig.ident = gsub("PND18pre","PND18",SSCs@meta.data$orig.ident)
jpeg(paste0(path,"/2A_tSNE.jpeg"), units="in", width=10, height=7,
     res=600)
TSNEPlot.1(object = SSCs,do.label = T, group.by = "orig.ident",
           do.return = TRUE, no.legend = T,#colors.use = singler.colors,
           pt.size = 1,label.size = 5,label.repel = T)+
        ggtitle("tSNE plot of all libraries")+
        theme(text = element_text(size=20),
              plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

#################################################################
#
# Figure 2 B) Composite tSNE showing all clusters
#
#################################################################

jpeg(paste0(path,"/2B_tSNE.jpeg"), units="in", width=10, height=7,
     res=600)
TSNEPlot.1(object = SSCs,do.label = T, group.by = "ident",
           do.return = TRUE, no.legend = T,colors.use = singler.colors,
           pt.size = 1,label.size = 5,label.repel = T, force = 5)+
        ggtitle("tSNE plot of all cell types")+
        theme(text = element_text(size=20),
              plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

#################################################################
#
# Figure 3 Marker gene heatmap for five germ cell types  =======
#
#################################################################
lname1 = load(file = "./data/SSCs_20180926.Rda");lname1
SSCs@meta.data$orig.ident = gsub("PND18pre","PND18",SSCs@meta.data$orig.ident)
table(SSCs@ident)
table(SSCs@meta.data$orig.ident)
major_cells <- c("Spermatogonia","Early Spermatocytes","Spermatocytes",
                 "Round Spermatids","Spermatids")
SSCs_spermato <- SubsetData(SSCs,ident.use = major_cells)

# 5.2.1 generate heatmap =======
# run 4.2.1 script
SSCs_spermato_markers = read.csv("./output/20181030/SSCs_spermato_markers.csv",
                                 header = T,stringsAsFactors =F)
SSCs_spermato_markers %>% head(20) %>% kable() %>% kable_styling()
top <-  SSCs_spermato_markers %>% group_by(cluster) %>% top_n(250, avg_logFC)
write.csv2(top[!duplicated(top$gene),],paste0(path,"/Fig_3_five_germ_cell_types.csv"))
g1 <- DoHeatmap.1(SSCs_spermato, SSCs_spermato_markers, 
                  Top_n = 250, 
                  col.low = "#0000FF",col.mid = "#FFFFFF",col.high = "#FF0000",
                  group.order = major_cells, 
                  ident.use = "five types of germ cells",
                  draw.line = TRUE, cex.row = 0,#cex.row = 6,
                  group.label.rot = T,remove.key =F, title.size = 12)

# 5.2.2. make color bar ========
# combine cell types and time points =======
# merge # all time point ...1 =======
SSCs_spermato@meta.data$Cell.Types = gsub("Early Spermatocytes",
                                     "Early_Spermatocytes",SSCs_spermato@meta.data$Cell.Types)
SSCs_spermato@meta.data$Cell.Types = gsub("Round Spermatids",
                                         "Round_Spermatids",SSCs_spermato@meta.data$Cell.Types)

SSCs_spermato@meta.data$Cell.Types = paste(SSCs_spermato@meta.data$Cell.Types,
                                          SSCs_spermato@meta.data$orig.ident)
table(SSCs_spermato@meta.data$Cell.Types)

SSCs_spermato <- SetAllIdent(object = SSCs_spermato, id = 'Cell.Types')
time_series <- c("PND06","PND14","PND18","PND25","PND30","Ad-depleteSp","Ad-Thy1")
# all time point ...1
time_points1 <- c(paste("Spermatogonia",time_series),
                  paste("Early_Spermatocytes", time_series),
                  paste("Spermatocytes",time_series[-1]),
                  paste("Round_Spermatids", time_series[-c(1,2)]),
                  paste("Spermatids", time_series[-c(1,2)]))
table(unique(as.vector(SSCs_spermato@ident)) %in% time_points1)
SSCs_spermato@ident = factor(SSCs_spermato@ident,levels = time_points1)
SSCs_spermato@ident = factor(SSCs_spermato@ident,levels = time_points1)
table(SSCs_spermato@ident)

# make bottom color bar==========
Cell.Types1 <- data.frame("Cell.Types" = as.vector(SSCs_spermato@ident),
                        row.names = names(SSCs_spermato@ident))
Cell.Types <- data.frame("Cell.Types" = as.vector(SSCs_spermato@meta.data$Cell.Types),
                          row.names = rownames(SSCs_spermato@meta.data))
table(rownames(Cell.Types1) == rownames(Cell.Types))
b_color_bar <- data.table::tstrsplit(Cell.Types1$Cell.Types,split = " ",keep = c(1,2))
b_color_bar = as.data.frame(b_color_bar,col.names = c("major_cells","time_points"))
b_color_bar$x <- 1:nrow(b_color_bar)
b_color_bar$time_points = factor(b_color_bar$time_points,levels = time_series)
head(b_color_bar)
table(b_color_bar[b_color_bar$major_cells == "Spermatocytes","time_points"])
table(b_color_bar$time_points)

color = c("#53B400","#00C094","#00B6EB","#A58AFF","#FB61D7","#F8766D","#C49A00")
MakeCorlorBar <- function(df, cell_type, color =NULL, remove.legend =F){
        df = df[(df$major_cells %in% cell_type),]
        g <- ggplot(data = df, aes(x, major_cells, fill = time_points)) +
                geom_tile()+
                theme_bw() + 
                theme(panel.border = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.line = element_blank(), 
                      axis.title.x=element_blank(),
                      axis.text.x=element_blank(),
                      axis.ticks.x=element_blank(),
                      axis.title.y=element_blank(),
                      axis.text.y=element_blank(),
                      axis.ticks.y=element_blank())
        if(remove.legend) g = g + theme(legend.position="none")
        if(!is.null(color)) {
                colors_fill = color_scheme(color = color,
                                           names = unique(df$time_points))
                g = g + scale_fill_manual(values = colors_fill)
        }
        g
}
# Make color scheme vector
# names must be ordered vector or factor
color_scheme <- function(color,names){
        df_names <- as.data.frame(table(names))
        df_names_Var <- df_names$Freq
        color_code <- as.character(unlist(mapply(rep, color, df_names_Var)))
        names(color_code) = names
        return(color_code)
        }

show_col(hue_pal()(7))
table(b_color_bar[b_color_bar$major_cells == "Spermatids","time_points"])
g_legend <- MakeCorlorBar(df = b_color_bar, cell_type = "Spermatogonia",
                          remove.legend = F, color = color)
g_Spermatogonia <- MakeCorlorBar(df = b_color_bar, cell_type = "Spermatogonia",
                                 color = color)
g_Early_Spermatocytes <- MakeCorlorBar(df = b_color_bar, cell_type = "Early_Spermatocytes",
                                 color = color)
g_Spermatocytes <- MakeCorlorBar(df = b_color_bar, cell_type = "Spermatocytes",
                                       color = color)
g_Round_Spermatids <- MakeCorlorBar(df = b_color_bar, cell_type = "Round_Spermatids",
                              color = color)
g_Spermatids <- MakeCorlorBar(df = b_color_bar, cell_type = "Spermatids",
                                 color = color)
library(ggpmisc)
library(gridExtra)
FindAllMarkers()
# save ggplot
path <- paste("./output",gsub("-","",Sys.Date()),sep = "/")
dir.create(path, recursive = T)

jpeg(paste0(path,"/3_heatmap_scale_blue.red.jpeg"), units="in", width=10, height=7,res=600)
g1 + theme(strip.text.x = element_text(margin=margin(t = 30, r = 0, b = 0, l = 0)))
dev.off()

# add legned
jpeg(paste0(path,"/3_heatmap_legend.jpeg"), units="in", width=10, height=7,res=600)
g_legend
dev.off()

# add major figure
jpeg(paste0(path,"/3_heatmap_scale.jpeg"), units="in", width=10, height=7,res=600)
g1 + theme(strip.text.x = element_text(margin=margin(t = 30, r = 0, b = 0, l = 0))) +
        annotation_custom(ggplotGrob(g_Spermatogonia),ymin = -3,ymax = 3)
dev.off()

# add Spermatogonia color bar
jpeg(paste0(path,"/3_heatmap_1.jpeg"), units="in", width=10, height=7,res=600)
g_Spermatogonia
dev.off()

# add Early_Spermatocytes color bar
jpeg(paste0(path,"/3_heatmap_2.jpeg"), units="in", width=10, height=7,res=600)
g_Early_Spermatocytes
dev.off()

# add g_Spermatocytes color bar
jpeg(paste0(path,"/3_heatmap_3.jpeg"), units="in", width=10, height=7,res=600)
g_Spermatocytes
dev.off()

# add g_Round_Spermatids color bar
jpeg(paste0(path,"/3_heatmap_4.jpeg"), units="in", width=10, height=7,res=600)
g_Round_Spermatids
dev.off()

# add g_Spermatids color bar
jpeg(paste0(path,"/3_heatmap_5.jpeg"), units="in", width=10, height=7,res=600)
g_Spermatids
dev.off()

#################################################################
#
# Figure 4 Marker gene, development gene, and discovered GOI heatmap for spermatogonia
#
#################################################################
# 6.4.1 load data =======
lname1 = load(file = "./data/SSCs_20180926.Rda"); lname1
lname1 = load(file = "./data/SSCs_20181001.Rda") #do.center = F, do.scale = T
lname1 = load(file = "./data/SSCs_20181030.Rda") #.scale with median
table(SSCs@ident)
SSCs@meta.data$orig.ident = gsub("PND18pre","PND18",SSCs@meta.data$orig.ident)
table(SSCs@meta.data$orig.ident)

Spermatogonia <- SubsetData(SSCs,ident.use = "Spermatogonia")

# 6.4.2 remove spermatid contaminated transcripts =======
Spermatid_Genes_to_Filter <- readxl::read_excel("./doc/Spermatid Genes to Filter.xlsx")
Spermatid_Genes = Spermatid_Genes_to_Filter$`Completely Spermatid Gene List To Remove (Spermatid/Other UMI Ratio > 20; and GSEA and tSNEs)`
genes.use = rownames(Spermatogonia@data); length(genes.use)
table(Spermatid_Genes %in% genes.use)
genes.use = genes.use[!(genes.use %in% Spermatid_Genes)]; length(genes.use)
Spermatogonia@data = Spermatogonia@data[genes.use, ]
Spermatogonia@scale.data = Spermatogonia@scale.data[genes.use, ]
# 6.4.3 merge ident and find development genes =======
table(Spermatogonia@meta.data$orig.ident)
Spermatogonia@meta.data$orig.ident = gsub("PND18|PND25|PND30",
                                     "PND18-30",Spermatogonia@meta.data$orig.ident)
Spermatogonia@meta.data$orig.ident = gsub("Ad-depleteSp|Ad-Thy1",
                                     "Adult",Spermatogonia@meta.data$orig.ident)
Spermatogonia <- SetAllIdent.1(object = Spermatogonia, id = 'orig.ident')
table(Spermatogonia@ident)
time_series <- c("PND06","PND14","PND18-30","Adult")
Spermatogonia@ident = factor(Spermatogonia@ident,levels = time_series)
TSNEPlot(Spermatogonia)
Spermatogonia_dev <- FindAllMarkers.UMI(Spermatogonia,only.pos = T,
                                        get.slot = "data")
table(Spermatogonia_dev$cluster)
write.csv(Spermatogonia_dev, file = paste0(path,"/Spermatogonia_dev_data_mean~.csv"))
Spermatogonia_dev = read.csv("./output/20181030/Spermatogonia_dev_data_mean~.csv")

# 6.4.4 generate heatmap =======
SSCs_spermato_markers = read.csv("./output/20181030/SSCs_spermato_markers.csv",
                                 header = T,stringsAsFactors =F)
table(SSCs_spermato_markers$cluster)
Spermatogonia_marker = SSCs_spermato_markers[(SSCs_spermato_markers$cluster 
                                              %in% "Spermatogonia"),"gene"]
Spermatogonia_dev2 = Spermatogonia_dev[(Spermatogonia_dev$pct.1>0.1),]
markers = MouseGenes(Spermatogonia,c(Spermatogonia_marker[1:20],"Gfra1", "Zbtb16",
                                     "Nanos2", "Utf1", "Sall4",
                                    "Id4", "Sohlh1", "Kit", "Lin28a","Dmrt1"))
markers = MouseGenes(Spermatogonia,c(Spermatogonia_marker[1:30]))
g2 <- DoHeatmap.1(Spermatogonia, marker_df = Spermatogonia_dev2, 
                  Top_n = 250, add.genes = markers,
                  col.low = "#FF00FF",col.mid = "#000000", col.high = "#FFFF00",
                  #col.low = "#FF00FF",col.mid = "#FFFFFF",col.high = "#FF0000",
                  group.order = time_series, 
                  ident.use = "Spermatogonia",
                  draw.line = TRUE, use.scaled = TRUE,
                  group.label.rot = T,cex.row = 0,remove.key =F,
                  title.size = 12)+
        theme(legend.text=element_text(size=7))
# make color bar
gene.use = rownames(Spermatogonia@scale.data);length(gene.use)
top <-  Spermatogonia_dev2 %>% group_by(cluster) %>% top_n(250, avg_logFC)
table(top$cluster)
df_markers = data.frame("gene" = markers,
                        "cluster" = "Spermatogonia\nMarkers")
Cell.Types = top[,c("gene","cluster")]
color_bar = rbind.data.frame(df_markers,Cell.Types)
colnames(color_bar) = c("gene","time_points")
color_bar = color_bar[!duplicated(color_bar$gene),]
color_bar = color_bar[color_bar$gene %in% gene.use,]
head(color_bar);dim(color_bar)
color_bar$major_cells = "Spermatogonia"
color_bar$time_points = factor(color_bar$time,levels = c("Spermatogonia\nMarkers",
                                                         time_series))
color_bar$x <- 1:nrow(color_bar)
head(color_bar)
write.csv2(color_bar,paste0(path,"/Fig_4_spermatogonia_GOI_mean~.csv"))
color = c("#386CB0","#53B400","#00C094","#A58AFF","#F8766D")#,"#FB61D7","#F8766D","#C49A00")

g_Spermatogonia <- MakeCorlorBar(df = color_bar, 
                                 cell_type = "Spermatogonia",color=color) +
                        coord_flip() + scale_x_reverse()

jpeg(paste0(path,"/4_Spermatogonia_dev_data_mean~.jpeg"), units="in", width=10, height=7,res=600)
g2
dev.off()

jpeg(paste0(path,"/4_Spermatogonia_mean~.jpeg"), units="in", width=3, height=7,res=600)
g_Spermatogonia
dev.off()
#################################################################
#
# Figure 5 Marker gene, development gene, and discovered GOI heatmap for spermatocytes
#
#################################################################
# 6.5.1 load data =======
lname1 = load(file = "./data/SSCs_20180926.Rda"); lname1
lname1 = load(file = "./data/SSCs_20181029.Rda") #.scale with median
table(SSCs@ident)
table(SSCs@meta.data$orig.ident)
#SSCs@meta.data$orig.ident = gsub("PND18pre","PND18",SSCs@meta.data$orig.ident)
Spermatocytes <- SubsetData(SSCs,ident.use = "Spermatocytes")

# 6.5.2 remove spermatid contaminated transcripts =======
Spermatid_Genes_to_Filter <- readxl::read_excel("./doc/Spermatid Genes to Filter.xlsx")
Spermatid_Genes = Spermatid_Genes_to_Filter$`Completely Spermatid Gene List To Remove (Spermatid/Other UMI Ratio > 20; and GSEA and tSNEs)`
genes.use = rownames(Spermatocytes@data); length(genes.use)
table(Spermatid_Genes %in% genes.use)
genes.use = genes.use[!(genes.use %in% Spermatid_Genes)]; length(genes.use)
Spermatocytes@data = Spermatocytes@data[genes.use, ]
Spermatocytes@scale.data = Spermatocytes@scale.data[genes.use, ]
# 6.5.3 merge ident and find development genes =======
table(Spermatocytes@meta.data$orig.ident)
table(Spermatocytes@ident)
Spermatocytes <- SetAllIdent.1(object = Spermatocytes, id = 'orig.ident')
time_series <- c("PND14","PND18","PND25","PND30","Ad-depleteSp","Ad-Thy1")
Spermatocytes@ident = factor(Spermatocytes@ident,levels = time_series)
TSNEPlot(Spermatocytes)

Spermatocytes_dev <- FindAllMarkers.UMI(Spermatocytes,only.pos = T,
                                        get.slot = "data")
table(Spermatocytes_dev$cluster)
write.csv(Spermatocytes_dev, file = paste0(path,"/Spermatocytes_dev_data_mean~.csv"))
#Spermatocytes_dev = read.csv2("./output/20181030/Spermatocytes_dev.data_mean~.csv")
# 6.5.4 generate heatmap =======
SSCs_spermato_markers = read.csv("./output/20180918/SSCs_spermato_markers.csv",
                                 header = T,stringsAsFactors =F)
table(SSCs_spermato_markers$cluster)
Spermatocytes_marker = SSCs_spermato_markers[(SSCs_spermato_markers$cluster 
                                              %in% "Spermatocytes"),"gene"]
Spermatocytes_dev2 = Spermatocytes_dev[(Spermatocytes_dev$pct.1>0.1),]

markers = MouseGenes(Spermatocytes,c(Spermatocytes_marker[1:30]))

g1 <- DoHeatmap.1(Spermatocytes, marker_df = Spermatocytes_dev2, 
                  Top_n = 250, add.genes = markers,
                  group.order = time_series, 
                  ident.use = "Spermatocytes",
                  draw.line = TRUE, use.scaled = TRUE,
                  group.label.rot = T,cex.row = 0,remove.key =F,
                  title.size = 12)+
        theme(legend.text=element_text(size=7))

# make color bar
gene.use = rownames(Spermatocytes@scale.data);length(gene.use)
top <-  Spermatocytes_dev2 %>% group_by(cluster) %>% top_n(250, avg_logFC)
table(top$cluster)
df_markers = data.frame("gene" = markers,
                        "cluster" = "Spermatocytes\nMarkers")
Cell.Types = top[,c("gene","cluster")]
color_bar = rbind.data.frame(df_markers,Cell.Types)
colnames(color_bar) = c("gene","time_points")
color_bar = color_bar[!duplicated(color_bar$gene),]
color_bar = color_bar[color_bar$gene %in% gene.use,]
head(color_bar);dim(color_bar)
color_bar$major_cells = "Spermatocytes"
color_bar$time_points = factor(color_bar$time,levels = c("Spermatocytes\nMarkers",
                                                         time_series))
color_bar$x <- 1:nrow(color_bar)
head(color_bar)
table(color_bar$time_points)
write.csv2(color_bar,paste0(path,"/Fig_5_Spermatocytes_GOI_mean~.csv"))

color = c("#BF5B17","#00C094","#00B6EB","#A58AFF","#FB61D7","#F8766D","#C49A00")

g_Spermatocytes <- MakeCorlorBar(df = color_bar, 
                                 cell_type = "Spermatocytes",color=color) +
        coord_flip() + scale_x_reverse()

jpeg(paste0(path,"/5_Spermatocytes_dev_data_mean~.jpeg"), units="in", width=10, height=7,res=600)
g1
dev.off()

jpeg(paste0(path,"/5_Spermatocytes_mean~.jpeg"), units="in", width=3, height=7,res=600)
g_Spermatocytes
dev.off()
#################################################################
#
# Figure S1 – Data Quality Graphs
#
#################################################################
SSCs_raw <- list()
SSCs_Seurat <- list()
samples <- c("PND06","PND14","PND18",#"PND18pre",
             "PND25","PND30","Ad-depleteSp","Ad-Thy1")
conditions <- c("first-wave","first-wave","first-wave",#"first-wave",
                "first-wave","first-wave","Adault-SSCs","Adault-SSCs")
for(i in 1:length(samples)){
        SSCs_raw[[i]] <- Read10X(data.dir = paste0("./data/",
                                                   samples[i],"/outs/filtered_gene_bc_matrices/mm10/"))
        colnames(SSCs_raw[[i]]) <- paste0(samples[i],"_",colnames(SSCs_raw[[i]]))
        SSCs_Seurat[[i]] <- CreateSeuratObject(SSCs_raw[[i]],
                                               min.cells = 3,
                                               min.genes = 200,
                                               names.delim = "_",
                                               project = "paula")
        SSCs_Seurat[[i]]@meta.data$conditions <- conditions[i]
}
SSCs <- Reduce(function(x, y) MergeSeurat(x, y, do.normalize = F), SSCs_Seurat)
remove(SSCs_raw);GC()
SSCs <- FilterCells(SSCs, subset.names = "nGene",
                    low.thresholds = 200,
                    high.thresholds = Inf) %>%
        NormalizeData() %>%
        ScaleData(display.progress = FALSE) %>%
        FindVariableGenes(do.plot = FALSE, display.progress = FALSE)
mito.genes <- grep(pattern = "^mt-", x = rownames(x = SSCs@data), value = TRUE)
percent.mito <- Matrix::colSums(SSCs@raw.data[mito.genes, ])/Matrix::colSums(SSCs@raw.data)
SSCs <- AddMetaData(object = SSCs, metadata = percent.mito, col.name = "percent.mito")

time_series <- c("PND06","PND14","PND18","PND25","PND30","Ad-depleteSp","Ad-Thy1")
SSCs@ident = factor(SSCs@ident,levels = time_series)
g1 <- VlnPlot.1(object = SSCs, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 1,
              point.size.use = 0, group.by = "ident",x.lab = "Raw data",
              x.lab.rot = T, do.return = T,return.plotlist = T)
SSCs1 <- FilterCells(object = SSCs, subset.names = c("nGene","nUMI","percent.mito"),
                    low.thresholds = c(500,2000, -Inf), 
                    high.thresholds = c(8000,125000, 0.15))
g2 <- VlnPlot.1(object = SSCs1, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 1,
              point.size.use = 0,x.lab = "After data filtration",
              x.lab.rot = T, do.return = T,return.plotlist = T)

jpeg(paste0(path,"/S1_QC.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(g1[[1]]+coord_flip()+
                     ylim(0,10000)+
                     theme(axis.title.y=element_blank(),
                           axis.text.y=element_blank(),
                           axis.ticks.y=element_blank()),
             g1[[2]]+coord_flip()+
                     scale_y_log10(limits =c(1000,300000))+
                     theme(axis.title.y=element_blank(),
                           axis.text.y=element_blank(),
                           axis.ticks.y=element_blank()),
             g1[[3]]+coord_flip()+
                     ylim(0,1)+
                     theme(axis.title.y=element_blank(),
                           axis.text.y=element_blank(),
                           axis.ticks.y=element_blank()),
             g2[[1]]+coord_flip()+
                     ylim(0,10000)+
                     theme(axis.title.y=element_blank(),
                           axis.text.y=element_blank(),
                           axis.ticks.y=element_blank()),
             g2[[2]]+coord_flip()+
                     scale_y_log10(limits =c(1000,300000))+
                     theme(axis.title.y=element_blank(),
                           axis.text.y=element_blank(),
                           axis.ticks.y=element_blank()),
             g2[[3]]+coord_flip()+
                     #ylim(0,1)+
                     theme(axis.title.y=element_blank(),
                           axis.text.y=element_blank(),
                           axis.ticks.y=element_blank()),
        nrow = 2,align="hv"
)
dev.off()

jpeg(paste0(path,"/S1_label.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(g1[[1]]+coord_flip(),
             g2[[1]]+coord_flip(),
             nrow = 2,align="hv")
dev.off()
#################################################################
#
# Figure S2A – tSNE plot of all clusters
#
#################################################################
lname1 = load(file = "./data/SSCs_20180926.Rda");lname1
SSCs <- SetAllIdent.1(object = SSCs, id = 'res.0.8')
jpeg(paste0(path,"/S2A_tSNE.jpeg"), units="in", width=10, height=7,
     res=600)
TSNEPlot.1(object = SSCs,do.label = T, group.by = "ident",
           do.return = TRUE, no.legend = T,#colors.use = singler.colors,
           pt.size = 1,label.size = 5,label.repel = T)+
        ggtitle("tSNE plot of all clusters")+
        theme(text = element_text(size=20),
              plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

#################################################################
#
# Figure S2B – dot plot of all clusters
#
#################################################################
GermCell_markers <- MouseGenes(SSCs,c("Gfra1","Zbtb16","Sall4","Dmrt1","Dazl","Kit","Cdca8",
                                      "Id4","Sycp3","Mtl5","Nxt1","Shcbp1l","Aurka","Lyzl1",
                                      "Acrv1","Hemgn","Txndc8","Tssk6","Oaz3","Prm2"))
Somatic_markers <- MouseGenes(SSCs,c("Col1a2","Acta2","Vcam1","Insl3","Laptm5",
                                     "Hbb-bt","Ptgds","Wt1"))
marker_order <- c(8,6,3,9,0,15,2,7,18,17,10,11,5,4,12,20,21,14,16,13,1,19)
SSCs@ident <- factor(x = SSCs@ident, levels = rev(marker_order)) # Relevel object@ident
markers.to.plot <- c(GermCell_markers,Somatic_markers)
jpeg(paste0(path,"/S2B_dotplot.jpeg"), units="in", width=10, height=7,
     res=600)
DotPlot(SSCs, genes.plot = rev(markers.to.plot),
        cols.use = c("blue","red"), x.lab.rot = T, plot.legend = T,
        dot.scale = 8, do.return = T)
dev.off()

#################################################################
#
# Figure S2B – dot plot of all cell types
#
#################################################################
SSCs <- SetAllIdent.1(object = SSCs, id = 'Cell.Types')
major_cells <- c("Spermatogonia","Early Spermatocytes","Spermatocytes",
                 "Round Spermatids","Spermatids","Smooth muscle",
                 "Endothelial & Hematopoietic cells", "Sertoli cells")
SSCs@ident <- factor(x = SSCs@ident, levels = rev(major_cells))
jpeg(paste0(path,"/S2C_dotplot.jpeg"), units="in", width=10, height=7,
     res=600)
DotPlot(SSCs, genes.plot = rev(markers.to.plot),
        cols.use = c("blue","red"), x.lab.rot = T, plot.legend = T,
        dot.scale = 8, do.return = T)
dev.off()

#################################################################
#
# Figure S3 Distribution profiles
#
#################################################################
lname1 = load(file = "./data/SSCs_20180926.Rda");lname1

g1 <- VlnPlot.1(object = SSCs, features.plot = "nGene", nCol = 1,x.lab = "",
              x.lab.rot = T, do.return = T, use.scaled = TRUE, group.by = "ident",
              point.size.use = 0,)+coord_flip() 
g2 <- VlnPlot.1(object = SSCs, features.plot = "nUMI", nCol = 1,,x.lab = "",
              x.lab.rot = T, do.return = T, use.scaled = TRUE, group.by = "ident",
              point.size.use = 0,)+coord_flip() 
g3 <- VlnPlot(object = SSCs, features.plot = "percent.mito", nCol = 1,
              x.lab.rot = T, do.return = T, use.scaled = TRUE, group.by = "ident",
              point.size.use = 0,)+coord_flip() 
jpeg(paste0(path,"/S3_Distribution_profiles.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(g1 +theme(axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank()),
          g2 +theme(axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank()),
          g3 + theme(axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank()),
          nrow = 1,align="hv")
dev.off()

jpeg(paste0(path,"/S3_label.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(g2,
          nrow = 1,align="hv")
dev.off()

#################################################################
#
# Figure S4 – tsne plot 
#
#################################################################

lname1 = load(file = "./data/SSCs_20180926.Rda");lname1
table(SSCs@ident)
major_cells <- c("Spermatogonia","Early Spermatocytes","Spermatocytes",
                 "Round Spermatids","Spermatids")
SSCs_spermato <- SubsetData(SSCs,ident.use = major_cells)
g <- list()
Somatic_markers <- MouseGenes(SSCs,c("Acta2","Col1a2","Vcam1","Insl3","Laptm5",
                                     "Hbb-bt","Ptgds","Wt1","Col1a1"))
GermCell_markers <- MouseGenes(SSCs,c("Gfra1","Zbtb16","Sall4","Dmrt1","Dazl","Kit","Cdca8",
                                      "Id4","Sycp3","Mtl5","Nxt1","Shcbp1l","Aurka","Lyzl1",
                                      "Acrv1","Hemgn","Txndc8","Tssk6","Oaz3","Prm2"))

g[[1]] = TSNEPlot.1(object = SSCs_spermato,do.label = T, group.by = "ident",
           do.return = TRUE, no.legend = T,colors.use = singler.colors,
           pt.size = 1,label.size = 5,label.repel = T, force = 5)+
        ggtitle("tSNE plot of five germ cells")+
        theme(text = element_text(size=20),
              plot.title = element_text(hjust = 0.5, face = "bold"))
# Somatic_markers 
g[[2]] = SingleFeaturePlot.1(object = SSCs, "Acta2", threshold = 1)
g[[3]] = SingleFeaturePlot.1(object = SSCs, "Col1a2")
g[[4]] = SingleFeaturePlot.1(object = SSCs, "Vcam1")
g[[5]] = SingleFeaturePlot.1(object = SSCs, "Insl3")
g[[6]] = SingleFeaturePlot.1(object = SSCs, "Laptm5")
g[[7]] = SingleFeaturePlot.1(object = SSCs, "Hbb-bt", threshold = 1.5)
g[[8]] = SingleFeaturePlot.1(object = SSCs, "Ptgds")
g[[9]] = SingleFeaturePlot.1(object = SSCs, "Wt1", threshold = 1.5)
g[[10]] = SingleFeaturePlot.1(object = SSCs, "Col1a1", threshold = 1.5)
# GermCell_markers
g[[11]] = SingleFeaturePlot.1(object = SSCs_spermato, "Gfra1", threshold = 0.1)
g[[12]] = SingleFeaturePlot.1(object = SSCs_spermato, "Zbtb16", threshold = 1.5)
g[[13]] = SingleFeaturePlot.1(object = SSCs_spermato, "Sall4", threshold = 1)
g[[14]] = SingleFeaturePlot.1(object = SSCs_spermato, "Dmrt1", threshold = 0.1)
g[[15]] = SingleFeaturePlot.1(object = SSCs_spermato, "Dazl", threshold = 2.0)
g[[16]] = SingleFeaturePlot.1(object = SSCs_spermato, "Kit", threshold = 0.1)
g[[17]] = SingleFeaturePlot.1(object = SSCs_spermato, "Id4", threshold = 0.8)
g[[18]] = SingleFeaturePlot.1(object = SSCs_spermato, "Cdca8", threshold = 1.0)
g[[19]] = SingleFeaturePlot.1(object = SSCs_spermato, "Sycp3", threshold = 1.5)
g[[20]] = SingleFeaturePlot.1(object = SSCs_spermato, "Mtl5", threshold = 1.0)
g[[21]] = SingleFeaturePlot.1(object = SSCs_spermato, "Nxt1", threshold = 1.8)
g[[22]] = SingleFeaturePlot.1(object = SSCs_spermato, "Shcbp1l", threshold = 1.5)
g[[23]] = SingleFeaturePlot.1(object = SSCs_spermato, "Aurka", threshold = 1.5)
g[[24]] = SingleFeaturePlot.1(object = SSCs_spermato, "Lyzl1", threshold = 1.5)
g[[25]] = SingleFeaturePlot.1(object = SSCs_spermato, "Acrv1", threshold = 1.5)
g[[26]] = SingleFeaturePlot.1(object = SSCs_spermato, "Hemgn", threshold = 1.5)
g[[27]] = SingleFeaturePlot.1(object = SSCs_spermato, "Txndc8", threshold = 1.5)
g[[28]] = SingleFeaturePlot.1(object = SSCs_spermato, "Tssk6", threshold = 1.5)
g[[29]] = SingleFeaturePlot.1(object = SSCs_spermato, "Oaz3", threshold = 1.5)
g[[30]] = SingleFeaturePlot.1(object = SSCs_spermato, "Prm2", threshold = 3.5)

file_names = c("five_germ_cells",Somatic_markers,GermCell_markers)
length(file_names)
for (i in c(1:30)) {
        ggplot2::ggsave(filename = paste0("S4_",file_names[i],".jpeg"),
                        plot = g[[i]], path = path,
               units="in", width=10, height=7, dpi=600
        )
}

#################################################################
#
# upload to GEO
#
#################################################################

path <- paste("./output",gsub("-","",Sys.Date()),sep = "/")
dir.create(path, recursive = T)

lname1 = load(file = "./data/SSCs_20180926.Rda");lname1
write.csv2(SSCs@meta.data,paste0(path,"/meta_data.csv"))



#################################################################
#
# the cell count numbers for PND18pre and PND18 including somatic cells
#
#################################################################

lname1 = load(file = "./data/SSCs_20180926.Rda");lname1
table(SSCs@meta.data$orig.ident)
table(SSCs@ident)
major_cells <- c("Spermatogonia","Early Spermatocytes","Spermatocytes",
                 "Round Spermatids","Spermatids","Sertoli cells",
                 "Endothelial & Hematopoietic cells", "Smooth muscle")
counts <- table(as.vector(SSCs@ident), SSCs@meta.data$orig.ident)
kable(counts) %>% kable_styling()
major_cells = major_cells[1:5]
time_series <- c("PND06","PND14","PND18",#"PND18pre",
                 "PND25","PND30","Ad-depleteSp","Ad-Thy1")
counts = counts[rev(major_cells),time_series]

jpeg(paste0(path,"/1_B_bar_chart.jpeg"), units="in", width=10, height=7,
     res=600)
par(mfrow=c(1, 1), mar=c(5, 4, 4, 2))
barplot(counts, main="Total numbers of germ cells vs. time-points",
        xlab="Time points", ylab="Cell numbers",
        col=singler.colors[8:4],
        legend = rownames(counts),
        args.legend = list(x = "topleft", bty = "n"))
dev.off()

#################################################################
#
# tSNE plot that is all grey dots EXCEPT for the PND18pre and PND18 samples in colored dots
#
#################################################################
lname1 = load(file = "./data/SSCs_20180926.Rda");lname1
SSCs <- SetAllIdent(SSCs, id = "orig.ident")
color_scheme <- gg_colors(SSCs)
color_scheme <- c("#d3d3d3","#d3d3d3","#d3d3d3","#d3d3d3",
                  "#00BFC4", "#C77CFF","#d3d3d3","#d3d3d3")
p <- TSNEPlot.1(object = SSCs,do.label = T, group.by = "orig.ident",
           do.return = TRUE, no.legend = T,colors.use = color_scheme,
           pt.size = 1,label.size = 5,label.repel = T)+
        ggtitle("tSNE plot of all libraries")+
        theme(text = element_text(size=20),
              plot.title = element_text(hjust = 0.5, face = "bold"))

jpeg(paste0(path,"/tSNE_PND18_and_pre.jpeg"), units="in", width=10, height=7,
     res=600)
p
dev.off()

PND18 <- SubsetData(SSCs, ident.use = c("PND18","PND18pre"))
p1 <- TSNEPlot.1(object = PND18,do.label = T, group.by = "orig.ident",
                do.return = TRUE, no.legend = T,colors.use = c("#00BFC4", "#C77CFF"),
                pt.size = 1,label.size = 5,label.repel = T)+
        ggtitle("tSNE plot of PND18 and PND18pre")+
        theme(text = element_text(size=20),
              plot.title = element_text(hjust = 0.5, face = "bold"))

jpeg(paste0(path,"/tSNE_PND18_and_pre~.jpeg"), units="in", width=10, height=7,
     res=600)
p1
dev.off()
