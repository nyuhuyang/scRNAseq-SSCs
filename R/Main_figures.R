####This code demonstrates how the main and supplementary figures were generated. 
library(Seurat)
library(dplyr)
library(kableExtra)
library(MAST)
source("./R/Seurat_functions.R")


# Figure 3 Marker gene heatmap for five germ cell types  =======
lname1 = load(file = "./data/SSCs_20180825.Rda");lname1
SSCs@meta.data$orig.ident = gsub("PND18pre","PND18",SSCs@meta.data$orig.ident)
table(SSCs@meta.data$orig.ident)
table(SSCs@ident)

#====
Spermatid_Genes_to_Filter <- readxl::read_excel("./doc/Spermatid Genes to Filter.xlsx")
Spermatid_Genes = Spermatid_Genes_to_Filter$`Completely Spermatid Gene List To Remove (Spermatid/Other UMI Ratio > 20; and GSEA and tSNEs)`
genes.use = rownames(SSCs@data); length(genes.use)
table(genes.use %in% Spermatid_Genes)
genes.use = genes.use[!(genes.use %in% Spermatid_Genes)]; length(genes.use)
SSCs@raw.data = SSCs@raw.data[genes.use, ]
SSCs@data = SSCs@data[genes.use, ]

#====
All_markers <- FindAllMarkers.UMI(SSCs,test.use = "MAST")
kable(All_markers[1:20,]) %>% kable_styling()
rownames(All_markers) = NULL
All_markers <- All_markers %>% select("gene", everything()) # Moving the last column to the start
path <- paste("./output",gsub("-","",Sys.Date()),sep = "/")
dir.create(path, recursive = T)
write.csv(All_markers,paste0(path,"/All_markers.csv"))

#====
major_cells <- c("Spermatogonia","Early Spermatocytes","Spermatocytes")#,
                 #"Round Spermatids","Spermatids")
top <- group_top_mutate(df = All_markers[(All_markers$cluster %in% major_cells),], major_cells)
Spermato <- SubsetData(SSCs, ident.use = major_cells)
top %>% head(20) %>% kable() %>% kable_styling()
g1 <- DoHeatmap.1(Spermato,top,Top_n = 15, 
            group.order = major_cells,ident.use = "five germ cell types",
            group.label.rot = T,cex.row = 8,remove.key =T,do.plot = F)
# add color bar
g1 + theme(strip.text.x = element_text(margin=margin(t = 30, r = 0, b = 0, l = 0))) +
        annotation_custom(ggplotGrob(g2),ymin = -3,ymax = -0.1)
#All_markers = read.csv(paste0(path,"/All_markers.csv"),header = T)

axis.text.x = element_text(margin=margin(5,5,10,5,"pt"))

