library(Seurat)
library(dplyr)
library(plyr)
library(janitor)
library(pheatmap)
source("./R/Seurat_functions.R")

#====== 4.1 load data  ==========================================
lnames = load(file = "./data/MCL_alignment.Rda")
table(MCL@ident)
idents <- as.data.frame(table(MCL@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("0 T cells",
                     "1 CD14 Monocytes",
                     "2 B cells",
                     "3 CD8 T cells",
                     "4 B cells",
                     "5 NK T cells",
                     "6 CD16 Monocytes",
                     "7 CD8 T cells",
                     "8 Dendritic Cells",
                     "9 Myeloid cells")
MCL@ident <- plyr::mapvalues(x = MCL@ident,
                             from = old.ident.ids,
                             to = new.cluster.ids)

MCL.subsets <- SplitCells(MCL)
MCL.subsets[[3]]
MCL.patient <- MCL.subsets[[1]]
MCL.normal <- MCL.subsets[[2]]
#====== 4.2 Organize Immgen data ===========================================
ImmGenDataV1 <- read.csv(file = "../scRNAseq-Lymphoma/data/ImmGenDataV1.csv")
ImmGenDataV2 <- read.csv(file = "../scRNAseq-Lymphoma/data/ImmGenDataV2.csv")
ImmGenData <- full_join(ImmGenDataV1, ImmGenDataV2, by = "GeneSymbol")
colnames(ImmGenData) <- sub('X.', '', colnames(ImmGenData)) # remove all X.
ImmGenData <- ImmGenData %>% clean_names()
ImmGenData <- ImmGenData[,-c(1,3,213,214)] # remove probesetid and description

# organize cell type
CellTypes <- colnames(ImmGenData)
CellTypes <- sub('mlp_', 'sc_mlp_', CellTypes) # stem cells
CellTypes <- sub('pro_b_', 'b_pro_', CellTypes) # B cells
CellTypes <- sub('pre_b_', 'b_pre_', CellTypes) # B cells
CellTypes <- sub('b1', 'b_b1', CellTypes) # B cells
CellTypes <- sub('pre_t_', 't_pre_', CellTypes) # T cells
CellTypes <- sub('nkt_', 't_nkt_', CellTypes) # T cells
CellTypes[c(186:198,228:244)] <- paste0("stro_",CellTypes[c(186:198,228:244)]) # stromal cells
CellTypes[248:254] <- paste0("nk_",CellTypes[248:254]) # Innate Lymphocytes
CellTypes <- sub("_",".",CellTypes)
colnames(ImmGenData) <- CellTypes
# reorder
CellTypes <- c(CellTypes[1],sort(CellTypes[-1]))
ImmGenData <- ImmGenData[CellTypes]

# sum up gene expression
colnames(ImmGenData)[1] <- "genesymbol"
ImmGenData$genesymbol <- toupper(ImmGenData$genesymbol)
ImmGenData$genesymbol <- gsub(" ","",ImmGenData$genesymbol)
system.time(ImmGenData <- aggregate(. ~ genesymbol, data=ImmGenData, FUN=sum))
rownames(ImmGenData) <- ImmGenData$genesymbol

# calculate ImmGenData averageExp
ImmGenData_short <- ImmGenData
ImmGenData_short <- ImmGenData_short[,-1]
ImmGenData_short <- as.data.frame(t(ImmGenData_short))
Major_CellTypes <- sub('\\..*', '', rownames(ImmGenData_short))
Full_names <- data.frame("b" = "B_cells",
                          "ba" = "Basophils",
                          "dc" = "Dendritic_cells",
                          "eo" = "Eosinophils",
                          "gn" = "Neutrophils",
                          "mc" = "Mast_cells",
                          "mf" = "Macrophages",
                          "mo" = "Monocytes",
                          "nk" = "Innate_Lymphocytes",
                          "sc" = "Stem_cells",
                          "stro" = "Stromal_cells",
                          "t" = "T_cells",
                          "tgd" = "gd_T_cells")
rownames(Full_names) <- "Major_CellTypes"
ImmGenData_short$Major_CellTypes <- t(Full_names[1,match(Major_CellTypes,names(Full_names))])
ImmGenData_short[,21753:21756]
system.time(ImmGenData_short <- aggregate(. ~ Major_CellTypes, data = ImmGenData_short, FUN=mean))
rownames(ImmGenData_short) <- ImmGenData_short$Major_CellTypes
ImmGenData_short <- as.data.frame(t(ImmGenData_short[,-1]))
ImmGenData_short$genesymbol <- rownames(ImmGenData_short)
head(ImmGenData_short[,(ncol(ImmGenData_short)-3):ncol(ImmGenData_short)])
ImmGenData.summary <- ImmGenData_short
#====== 4.3 Identify Cell Types by Spearman correlation ==================================
Identify_Cell_Types_Spearman <- function(object, gendata, hv.gene.num = 2500, cluster_rows=F,
                                         cluster_cols = F,fontsize_row = 15,fontsize_col = 15,
                                         fontsize =20,title = ""){
        "
        Calculate Average Expression of each ident of seurat object,
        Calculate spearman correlation between above results and provided gendata dataset
        "
        if(class(object) != "seurat") {
                stop(paste("Error : ", object, " is not a seurat object"))
        }
        if(class(gendata) != "data.frame") {
                stop(paste("Error : ", gendata, " is not a data frame"))
        }
        object <- FindVariableGenes(object = object, mean.function = ExpMean, dispersion.function = LogVMR, 
                                 do.plot = FALSE)
        hv.genes <- head(rownames(object@hvg.info), hv.gene.num)
        object.AverageExp <- AverageExpression(object, genes.use = hv.genes)
        object.AverageExp$genesymbol <- toupper(rownames(object.AverageExp))
        object.Exp <- object.AverageExp[,c(ncol(object.AverageExp),
                                           1:(ncol(object.AverageExp)-1))]# move last column to the first
        print("Merge genes expression:")
        table(gendata$genesymbol %in% object.Exp$genesymbol)
        # merge ===============
        object.Exp_gendata <- inner_join(object.Exp, gendata, by = "genesymbol")
        object.Exp_gendata <- object.Exp_gendata[order(object.Exp_gendata$genesymbol),]
        rownames(object.Exp_gendata) <- object.Exp_gendata$genesymbol
        object.Exp_gendata <- object.Exp_gendata[,-1]
        # Spearman correlation normal ==================
        c <- cor(object.Exp_gendata, method="spearman") # or naive_matrix
        diag(c) <-NA
        ident_num <- length(levels(object@ident))
        object_c_gendata <- c[(ident_num+1):nrow(c),1:ident_num]
        pheatmap::pheatmap(object_c_gendata,cex=.9,
                 cluster_rows= cluster_rows,
                 cluster_cols = cluster_cols,
                 fontsize_row = fontsize_row,
                 fontsize_col = fontsize_col,
                 fontsize = fontsize,
                 main = title)
        actual_cell_types <- apply(object_c_gendata, 2, which.max)
        rename_ident <- data.frame("cluster.num" = 0:(ncol(object_c_gendata)-1),
                                   "old.ident.ids" = colnames(object_c_gendata),
                                   "new.cluster.ids" = rownames(object_c_gendata)[actual_cell_types])
        print(rename_ident)
        return(rename_ident)
}

# correlate with ImmGenData_short
MCL.normal.rename_ident <- Identify_Cell_Types_Spearman(object = MCL.normal, 
                                       gendata = ImmGenData.summary,
                                       cluster_rows=T,cluster_cols = T,
                                       title = "Spearman correlation: MCL.normal vs ImmGen Data")
MCL.patient.rename_ident <- Identify_Cell_Types_Spearman(object = MCL.patient, 
                                       gendata = ImmGenData.summary,
                                       cluster_rows=T,cluster_cols = T,
                                        title = "Spearman correlation: MCL.patient vs ImmGen Data")
# correlate with ImmGenData
MCL.normal.rename_ident <- Identify_Cell_Types_Spearman(object = MCL.normal, 
                                                           gendata = ImmGenData,
                                                           cluster_rows=T,cluster_cols = T,
                                                           title = "Spearman correlation: MCL.normal vs ImmGen Data")
MCL.patient.rename_ident <- Identify_Cell_Types_Spearman(object = MCL.patient, 
                                                       gendata = ImmGenData,
                                                       cluster_rows=T,cluster_cols = T,
                                                       title = "Spearman correlation: MCL.patient vs ImmGen Data")
#===== 4.4  A table with the number of cells of each cluster and subcluster ======
# We can also compare proportional shifts in the data. As can be seen in the barplot, 
table(MCL.normal@ident)/MCL.normal@data@Dim[2]
table(MCL.patient@ident)/MCL.patient@data@Dim[2]