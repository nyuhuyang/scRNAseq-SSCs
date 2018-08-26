library(SingleR)
library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
AverageExpression ()
#====== 2.1 load data  ==========================================
lname1 = load(file = "./data/SSCs_20180822.Rda");lname1
lname2 = load(file='./data/GeneSets/Ref_GSE43717.RData');lname2
Ref_GSE43717$name
length(Ref_GSE43717$types)
length(unique(Ref_GSE43717$types))
length(unique(Ref_GSE43717$main_types))
length(Ref_GSE43717$sd.thres)
DimPlot(object = SSCs, reduction.use = "tsne", no.legend = TRUE,
        do.return = TRUE,vector.friendly = F, pt.size = 1,
        do.label = TRUE,label.size = 8, group.by = "ident") + 
        ggtitle("Cluster ID") + 
        theme(plot.title = element_text(hjust = 0.5))
SSCs@project.name = "Paula"
#====== 2.2 Create Singler Object  ==========================================
singler = CreateSinglerObject(as.matrix(SSCs@data), annot = NULL, 
                              project.name=SSCs@project.name,
                              min.genes = 500,technology = "10X", species = "Mouse",
                              ref.list = list(Ref_GSE43717), normalize.gene.length = F, 
                              variable.genes = "de",
                              fine.tune = F, do.signatures = F, clusters = SSCs@ident)
GC()
singler$meta.data$orig.ident = SSCs@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler$meta.data$xy = SSCs@dr$tsne@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = SSCs@ident # the Seurat clusters (if 'clusters' not provided)
save(singler,file="./output/singler_20180822.RData")
#====== 2.3 SingleR specifications ==========================================
# Step 1: Spearman coefficient
lnames = load(file = "./output/singler_20180822.RData")
lnames
singler$seurat = NULL # (optional)
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
out+  ggtitle("Supervised sub-cell type labeling by immgen, GSE43717 and GSE83264")+
        theme(text = element_text(size=20),
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold"))
# main types-------
out = SingleR.PlotTsne.1(singler$singler[[1]]$SingleR.single.main,
                         singler$meta.data$xy,do.label=T,
                         do.letters = F,labels = singler$singler[[1]]$SingleR.single.main$labels,
                         label.size = 5, dot.size = 2,do.legend = F,alpha = 1,
                         label.repel = T,force=25)
out +  ggtitle("Supervised cell type labeling by immgen and GSE43717")+
        theme(text = element_text(size=20),
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold"))
#====== 2.4 compared to the original identities ==========================================
# cell number
singler$meta.data$orig.ident = gsub("Ad-","zAd-",singler$meta.data$orig.ident)
singler$meta.data$orig.ident = gsub("PND18pre","PND18",singler$meta.data$orig.ident)
counts <- table(singler$singler[[1]]$SingleR.single.main$labels,
                singler$meta.data$orig.ident)
kable(counts) %>% kable_styling()

barplot(counts, main="Total numbers of cell types vs. time-points",
        xlab="Number of Gears", col=singler.colors,
        legend = rownames(counts),
        args.legend = list(x = "top", bty = "n"))
# cell percentage
prop.table(x = table(singler$singler[[1]]$SingleR.single.main$labels,
                     SSCs@meta.data$orig.ident),margin = 2) %>%
        kable  %>% kable_styling()
# total cell number
table(singler$meta.data$orig.ident) %>% t() %>% kable() %>% kable_styling()

#====== 2.5 Rename ident ==========================================
table(names(SSCs@ident) == rownames(singler$singler[[1]]$SingleR.single.main$labels))

ident.use <- as.factor(as.character(singler$singler[[1]]$SingleR.single.main$labels))
names(ident.use) = rownames(singler$singler[[1]]$SingleR.single.main$labels)
SSCs@ident <- ident.use

TSNEPlot.1(object = SSCs,do.label = T, group.by = "ident", 
         do.return = TRUE, no.legend = T,colors.use = singler.colors,
         pt.size = 1,label.size = 5,label.repel = T,force=25)+
        ggtitle("Supervised cell type labeling by GSE43717")+
        theme(text = element_text(size=20),							
              plot.title = element_text(hjust = 0.5, face = "bold"))

save(SSCs, file = "./data/SSCs_20180822.Rda")

#------- 2.6 Store cell names by identities------------------------
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
save(singler_labels_id, file="./output/singler_labels.RData")
head(SSC_labels_id <- factor2list(cell.names = names(SSCs@ident),
                                  ident = as.character(SSCs@ident)),3)