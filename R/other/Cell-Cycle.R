library(Seurat)

# Read in the expression matrix The first row is a header row, the first
# column is rownames
lnames = load(file = "./data/SSCs_alignment.Rda")
# Also read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "~/Downloads/seurat_resources/regev_lab_cell_cycle_genes.txt")

# We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- MouseGenes(SSCs,cc.genes[1:43])
g2m.genes <- MouseGenes(SSCs,cc.genes[44:97])

SSCs <- FindVariableGenes(object = SSCs, mean.function = ExpMean, dispersion.function = LogVMR, 
                          do.plot = FALSE)
hv.genes <- head(rownames(SSCs@hvg.info), 1000)
SSCs <- RunPCA(object = SSCs, pc.genes = hv.genes, pcs.compute = 100, do.print = TRUE, 
               pcs.print = 1:5, genes.print = 5)
PCElbowPlot(object = SSCs, num.pc = 100)
PCHeatmap(SSCs, pc.use = c(1:3, 25:30), cells.use = 500, do.balanced = TRUE)
PCAPlot(object = SSCs)
# Assign Cell-Cycle Scores
SSCs <- CellCycleScoring(object = SSCs, s.genes = s.genes, g2m.genes = g2m.genes, 
                         set.ident = TRUE)
# view cell cycle scores and phase assignments
head(x = SSCs@meta.data)
# Visualize the distribution of cell cycle markers across
RidgePlot(object = SSCs, features.plot = MouseGenes(SSCs,c("PCNA", "TOP2A", "MCM6", "MKI67")), 
          nCol = 2)

# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells
# separate entirely by phase
SSCs <- RunPCA(object = SSCs, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)
PCAPlot(object = SSCs)
cd_genes <- list(c('CD79B','CD79A','CD19','CD180','CD200','CD3D','CD2',
                   'CD3E','CD7','CD8A','CD14','CD1C','CD68','CD9','CD247'
))
SSCs <- AddModuleScore(object = SSCs, genes.list = cd_genes,ctrl.size = 5,
        enrich.name = 'CD_Genes'
)
head(x = SSCs@meta.data)

# Regress out cell cycle scores during data scaling
SSCs <- ScaleData(object = SSCs, vars.to.regress = c("S.Score", "G2M.Score"), 
                    display.progress = T)
# Now, a PCA on the variable genes no longer returns components associated
# with cell cycle
SSCs <- RunPCA(object = SSCs, pc.genes = SSCs@var.genes, genes.print = 10)
PCAPlot(object = SSCs)

# Alternate Workflow
# regressing out the difference between the G2M and S phase scores
SSCs@meta.data$CC.Difference <- SSCs@meta.data$S.Score - SSCs@meta.data$G2M.Score
SSCs <- ScaleData(object = SSCs, vars.to.regress = "CC.Difference", 
                    display.progress = T)

# cell cycle effects strongly mitigated in PCA
SSCs <- RunPCA(object = SSCs, pc.genes = SSCs@var.genes, genes.print = 10)
PCAPlot(object = SSCs)

# when running a PCA on cell cycle genes, actively proliferating cells
# remain distinct from G1 cells however, within actively proliferating
# cells, G2M and S phase cells group together
SSCs <- RunPCA(object = SSCs, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)
PCAPlot(object = SSCs)
