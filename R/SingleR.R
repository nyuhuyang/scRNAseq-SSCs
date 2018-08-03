library(SingleR)
library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)
source("../R/Seurat_functions.R")

#====== 3.1 Create Singler Object  ==========================================
lnames = load(file = "./data/SSCs_label.Rda")
lnames
lname = load(file='./data/GeneSets/GSE43717_GSE83264.RData') 
lname
SSCs@meta.data$orig.ident = gsub("PND18pre","PND18",SSCs@meta.data$orig.ident)
ref_GSE43717_GSE83264$name
length(ref_GSE43717_GSE83264$types)
length(unique(ref_GSE43717_GSE83264$types))
length(unique(ref_GSE43717_GSE83264$main_types))
#pca
DimPlot(object = SSCs, reduction.use = "tsne", no.legend = TRUE,
        do.return = TRUE,vector.friendly = F, pt.size = 1,
        do.label = TRUE,label.size = 8, group.by = "ident") + 
        ggtitle("Cluster ID") + 
        theme(plot.title = element_text(hjust = 0.5))
singler = CreateSinglerObject(as.matrix(SSCs@data), annot = NULL, 
                              project.name=SSCs@project.name,
                              min.genes = 500,technology = "10X", species = "Mouse",
                              ref.list = list(GSE43717_GSE83264), normalize.gene.length = F, 
                              variable.genes = "de",
                              fine.tune = F, do.signatures = F, clusters = NULL)
GC()
singler$meta.data$orig.ident = SSCs@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler$meta.data$xy = SSCs@dr$tsne@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = SSCs@ident # the Seurat clusters (if 'clusters' not provided)
save(singler,file="./output/singler_SSCs.RData")
#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
lnames = load(file = "./output/singler_SSCs.RData")
lnames
singler$seurat = SSCs # (optional)
SingleR.DrawScatter(sc_data = singler$seurat@data,cell_id = 10, 
                    ref = immgen, sample_id = 232)

# Step 2: Multiple correlation coefficients per cell types are aggregated 
# to provide a single value per cell type per single-cell. 
# In the examples below we use the 80% percentile of correlation values.
# for visualization purposes we only present a subset of cell types (defined in labels.use)
out = SingleR.DrawBoxPlot(sc_data = singler$seurat@data,cell_id = 10, 
                          ref = immgen,main_types = T,
                          labels.use=c('B cells','T cells','DC','Macrophages','Monocytes','NK cells',
                                       'Mast cells','Neutrophils','Fibroblasts','Endothelial cells'))
print(out$plot)
SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single.main, top.n = Inf,
                    clusters = singler$meta.data$orig.ident)
#Or by all cell types (showing the top 50 cell types):
SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single, top.n = 50,
                    clusters = singler$meta.data$orig.ident)
SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single.main,top.n = 50,
                    normalize = F,clusters = singler$meta.data$orig.ident)
#Next, we can use the fine-tuned labels to color the t-SNE plot:

out = SingleR.PlotTsne.1(singler$singler[[1]]$SingleR.single,
                         singler$meta.data$xy,do.label=T,
                         do.letters = F,labels = singler$singler[[1]]$SingleR.single$labels,
                         label.size = 5, dot.size = 1,do.legend = F,alpha = 1,
                         label.repel = T,force=2)
out$p+  ggtitle("Supervised sub-cell type labeling by immgen, GSE43717 and GSE83264")+#ggplot title
        theme(text = element_text(size=20),     #larger text including legend title
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) #title in middle
# main types-------
out = SingleR.PlotTsne.1(singler$singler[[1]]$SingleR.single.main,
                         singler$meta.data$xy,do.label=T,
                         do.letters = F,labels = singler$singler[[1]]$SingleR.single.main$labels,
                         label.size = 5, dot.size = 1,do.legend = F,alpha = 1,
                         label.repel = T,force=2)
out$p+  ggtitle("Supervised cell type labeling by GSE43717 and GSE83264")+#ggplot title
        theme(text = element_text(size=20),     #larger text including legend title
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) #title in middle

output <- SplitSingleR.PlotTsne(singler = singler, split.by = "conditions",
                              return.plots=T,do.label=T,do.legend = F,alpha = 1,
                              label.repel = F, force=2)
output[[1]]
output[[2]]


#Finally, we can also view the labeling as a table compared to the original identities:

# cell number
kable(table(singler$singler[[1]]$SingleR.single.main$labels,
            singler$meta.data$orig.ident)) %>%
        kable_styling()
# cell percentage
prop.table(x = table(singler$singler[[1]]$SingleR.single.main$labels,
                     SSCs@meta.data$orig.ident),margin = 2) %>%
        kable  %>% kable_styling()
# total cell number
table(singler$meta.data$orig.ident) %>% t() %>% kable() %>% kable_styling()

# Rename ident
table(names(SSCs@ident) == rownames(singler$singler[[1]]$SingleR.single.main$labels))

ident.use <- as.factor(as.character(singler$singler[[1]]$SingleR.single.main$labels))
names(ident.use) = rownames(singler$singler[[1]]$SingleR.single.main$labels)
SSCs@ident <- ident.use

TSNEPlot(object = SSCs,do.label = F, group.by = "ident", 
         do.return = TRUE, no.legend = T,
         pt.size = 1,label.size = 8 )+
        ggtitle("Supervised cell type labeling by GSE43717 and GSE83264")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5))

save(SSCs, file = "./data/SSCs_suplabel.Rda")

cells_prop <- as.data.frame(table(singler$singler[[1]]$SingleR.single$labels,
                                              singler$meta.data$orig.ident))
cells_prop <- dcast(cells_prop,Var2~Var1)
rownames(cells_prop) = cells_prop$Var2
cells_prop <- cells_prop[,-1]
total  <- colSums(cells_prop)
Cells_prop <- cells_prop/total
Cells_prop %>% kable %>% kable_styling()

All_cells <- as.data.frame(table(SSCs@meta.data$orig.ident))
SpermatogonialSC_cells$total <- All_cells$Freq
SpermatogonialSC_cells$percentage <- SpermatogonialSC_cells$Freq/SpermatogonialSC_cells$total
SpermatogonialSC_cells