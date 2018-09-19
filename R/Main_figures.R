####This code demonstrates how the main and supplementary figures were generated. 
library(Seurat)
library(SingleR)
library(dplyr)
library(kableExtra)
library(MAST)
source("./R/Seurat_functions.R")

path <- paste("./output",gsub("-","",Sys.Date()),sep = "/")
dir.create(path, recursive = T)
#################################################################
#
# Figure 1 B) Germ cell composition per library (bar chart)
#
#################################################################

lname1 = load(file = "./data/SSCs_20180825.Rda");lname1
SSCs@meta.data$orig.ident = gsub("PND18pre","PND18",SSCs@meta.data$orig.ident)
table(SSCs@meta.data$orig.ident)
table(SSCs@ident)
major_cells <- c("Spermatogonia","Early Spermatocytes","Spermatocytes",
                 "Round Spermatids","Spermatids","Sertoli cells",
                 "Endothelial & Hematopoietic cells", "Smooth muscle")
SSCs_spermato <- SubsetData(SSCs,ident.use = major_cells[1:5])
SSCs_spermato@meta.data$orig.ident = sub("PND18pre","PND18",SSCs_spermato@meta.data$orig.ident)
counts <- table(as.vector(SSCs_spermato@ident), SSCs_spermato@meta.data$orig.ident)
kable(counts) %>% kable_styling()
major_cells = major_cells[1:5]
time_series <- c("PND06","PND14","PND18","PND25","PND30","Ad-depleteSp","Ad-Thy1")
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
lname1 = load(file = "./data/SSCs_20180825.Rda");lname1
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
# Figure 2 B) Composite tSNE showing all clusters
#
#################################################################
lname1 = load(file = "./data/SSCs_20180822.Rda");lname1
jpeg(paste0(path,"/2B_tSNE_1.jpeg"), units="in", width=10, height=7,
     res=600)
TSNEPlot.1(object = SSCs,do.label = T, group.by = "ident",
           do.return = TRUE, no.legend = T,colors.use = singler.colors,
           pt.size = 1,label.size = 5,label.repel = T, force = 5)+
        ggtitle("tSNE plot of all clusters")+
        theme(text = element_text(size=20),
              plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

#################################################################
#
# Figure 3 Marker gene heatmap for five germ cell types  =======
#
#################################################################
lname1 = load(file = "./data/SSCs_20180825.Rda");lname1
SSCs@meta.data$orig.ident = gsub("PND18pre","PND18",SSCs@meta.data$orig.ident)
major_cells <- c("Spermatogonia","Early Spermatocytes","Spermatocytes",
                 "Round Spermatids","Spermatids")
SSCs_spermato <- SubsetData(SSCs,ident.use = major_cells)

# 5.2.1 generate heatmap =======
# run 4.2.1 script
SSCs_spermato_markers = read.csv("./output/20180918/SSCs_spermato_markers.csv",
                                 header = T,stringsAsFactors =F)
SSCs_spermato_markers %>% head(20) %>% kable() %>% kable_styling()
marker = MouseGenes(SSCs_spermato,c("Gfra1", "Zbtb16", "Nanos2", "Utf1", "Sall4",
                                    "Id4", "Sohlh1", "Kit", "Lin28a", "Pax7", "Dmrt1"))
g1 <- DoHeatmap.1(SSCs_spermato,SSCs_spermato_markers, 
                  Top_n = 15, 
                  group.order = major_cells, 
                  ident.use = "five types of germ cells",
                  group.label.rot = T,cex.row = 6,remove.key =T,title.size = 12)

# 5.2.2. make color bar ========
# combine cell types and time points =======
ident.use <- data.frame("ident.use" = as.vector(SSCs_spermato@ident),
                        row.names = names(SSCs_spermato@ident))
SSCs_spermato <- AddMetaData(SSCs_spermato, ident.use,"ident.use")
SSCs_spermato@meta.data$ident.use = gsub("Early Spermatocytes",
                                     "Early_Spermatocytes",SSCs_spermato@meta.data$ident.use)
SSCs_spermato@meta.data$ident.use = gsub("Round Spermatids",
                                         "Round_Spermatids",SSCs_spermato@meta.data$ident.use)

SSCs_spermato@meta.data$ident.use = paste(SSCs_spermato@meta.data$ident.use,
                                          SSCs_spermato@meta.data$orig.ident)
table(SSCs_spermato@meta.data$ident.use)
# 5.2.4 merge # all time point ...1 =======
#SSCs_spermato@meta.data$ident.use = gsub("Early Spermatocytes Ad-.*$",
#                                     "Early Spermatocytes Adult",SSCs_spermato@meta.data$ident.use)
#SSCs_spermato@meta.data$ident.use = gsub("Early Spermatocytes PND06|Early Spermatocytes PND14",
#                                     "Early Spermatocytes PND06-14",SSCs_spermato@meta.data$ident.use)
#SSCs_spermato@meta.data$ident.use = gsub("Early Spermatocytes PND25|Early Spermatocytes PND30",
#                                     "Early Spermatocytes PND25-30",SSCs_spermato@meta.data$ident.use)
SSCs_spermato <- SetAllIdent(object = SSCs_spermato, id = 'ident.use')
time_series <- c("PND06","PND14","PND18","PND25","PND30","Ad-depleteSp","Ad-Thy1")
# all time point ...1
time_points1 <- c(paste("Spermatogonia",time_series),
                  paste("Early_Spermatocytes", time_series),
                  paste("Spermatocytes",time_series[-1]),
                  paste("Round_Spermatids", time_series[-c(1,2)]),
                  paste("Spermatids", time_series[-c(1,2)]))
table(unique(as.vector(SSCs_spermato@ident)) %in% time_points1)
SSCs_spermato@ident = factor(SSCs_spermato@ident,levels = time_points1)
table(SSCs_spermato@ident)
# make color bar
ident.use1 <- data.frame("ident.use" = as.vector(SSCs_spermato@ident),
                        row.names = names(SSCs_spermato@ident))
rownames(ident.use1) == rownames(ident.use)
color_bar <- data.table::tstrsplit(ident.use1$ident.use,split = " ",keep = c(1,2))
color_bar = as.data.frame(color_bar,col.names = c("major_cells","time_points"))
color_bar$x <- 1:nrow(color_bar)
color_bar$time_points = factor(color_bar$time,levels = time_series)
head(color_bar)
table(color_bar[color_bar$major_cells == "Spermatocytes","time_points"])
table(color_bar$time_points)
library(plyr)
color = c("#53B400","#00C094","#00B6EB","#A58AFF","#FB61D7","#F8766D","#C49A00")
MakeCorlorBar <- function(df, cell_type, color =NULL, remove.legend =T){
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
                                           names = df$time_points)
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

library(scales)
show_col(hue_pal()(7))
table(color_bar[color_bar$major_cells == "Spermatids","time_points"])
g_legend <- MakeCorlorBar(df = color_bar, cell_type = "Spermatogonia",
                          remove.legend = F, color = color)
g_Spermatogonia <- MakeCorlorBar(df = color_bar, cell_type = "Spermatogonia",
                                 color = color)
g_Early_Spermatocytes <- MakeCorlorBar(df = color_bar, cell_type = "Early_Spermatocytes",
                                 color = color)
g_Spermatocytes <- MakeCorlorBar(df = color_bar, cell_type = "Spermatocytes",
                                       color = color)
g_Spermatids <- MakeCorlorBar(df = color_bar, cell_type = "Spermatids",
                                 color = color)
library(ggpmisc)
library(gridExtra)
# save ggplot
path <- paste("./output",gsub("-","",Sys.Date()),sep = "/")
dir.create(path, recursive = T)

jpeg(paste0(path,"/3_heatmap_scale.jpeg"), units="in", width=10, height=7,res=600)
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
