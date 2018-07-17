library(SingleR)
library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)
source("./R/Seurat_functions.R")

#====== 3.1 Create Singler Object  ==========================================
lnames = load(file = "./data/SSCs_alignment_8.Rda")
lnames
lname = load(file = "./data/GeneSets/ref_E_GEOD_43717.RData")
lname
singler = CreateSinglerObject(as.matrix(SSCs@data), annot = NULL, project.name=SSCs@project.name,
                              min.genes = 500,technology = "10X", species = "Mouse", citation = "",
                              ref.list = list(ref_E_GEOD_43717), normalize.gene.length = F, variable.genes = "de",
                              fine.tune = T, do.signatures = T, clusters = NULL)
GC()
singler$seurat = SSCs # (optional)
singler$meta.data$orig.ident = SSCs@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler$meta.data$xy = SSCs@dr$tsne@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = SSCs@ident # the Seurat clusters (if 'clusters' not provided)
save(singler,file="./data/singler_SSCs.RData")
#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
lnames = load(file = "./data/singler_SSCs.RData")
lnames
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
SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single,top.n = 50,
                    normalize = F,clusters = singler$meta.data$orig.ident)
#Next, we can use the fine-tuned labels to color the t-SNE plot:
       
out = SingleR.PlotTsne.1(singler$singler[[1]]$SingleR.single,
                       singler$meta.data$xy,do.label=T,
                       do.letters = F,labels = singler$singler[[1]]$SingleR.single$labels,
                       label.size = 4, dot.size = 3,do.legend = F,alpha = 1,
                       label.repel = T,force=2)
out$p

SplitSingleR.PlotTsne(singler = singler, split.by = "conditions",do.label=T,
                      do.letters = T,do.legend = FALSE,force=2)
output <- SplitSingleR.PlotTsne(singler = singler, split.by = "conditions",
                              return.plots=T,do.label=T,do.legend = F,alpha = 0.5,
                              label.repel = F, force=2)
output[[1]]
output[[2]]
#Finally, we can also view the labeling as a table compared to the original identities:

kable(table(singler$singler[[1]]$SingleR.single$labels,
            singler$meta.data$orig.ident)) %>%
        kable_styling()
kable(table(singler$meta.data$orig.ident,
            singler$singler[[1]]$SingleR.single.main$labels)) %>%
        kable_styling()
kable(table(singler$meta.data$orig.ident,singler$seurat@ident)) %>%
        kable_styling()
