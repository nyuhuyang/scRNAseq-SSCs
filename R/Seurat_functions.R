# TRUE Global Environment

#' Produce a contingency table with certain gene expression at different ident
#'
#' @param object seurat object
#' @param subset.name Parameter to subset on. Eg, the name of a gene, PC1, ... 
#' Any argument that can be retreived using FetchData
#' @export Cells_ident nX2 dataframe with orig.ident at column1, Freq at column2
#' @examples
#' CountsbyIdent(SSCs,"Gfra1")
CountsbyIdent <- function(object,subset.name,...){
        "Generate nX2 dataframe with orig.ident at column1, Freq at column2"
        cells.use <- WhichCells(object=object,subset.name=subset.name,...)
        Cells.use <- data.frame("cells"=cells.use,"ident"=sub('_.*', '', cells.use))
        Cells_ident <- as.data.frame(table(Cells.use$ident))
        colnames(Cells_ident)[2] <- subset.name
        return(Cells_ident)
}


#' a supprting function for SingleFeaturePlot.1 and FeatureHeatmap.1
#' Change ggplot color scale to increase contrast gradient
#' #https://github.com/satijalab/seurat/issues/235
#' @param p ggplot object
#' @param alpha.use Define transparency of points
#' @param gradient.use Change fill and colour gradient values
#' @param scaled.expression.threshold Define lower limit of scaled gene expression level
#' @export p ggplot object
ChangeColorScale <- function(p, alpha.use = 0.8,
                             gradient.use = c("yellow", "red"),
                             scaled.expression.threshold = 0) {
        # Order data by scaled gene expresion level
        # Compute maximum value in gene expression
        if (length(p$data$scaled.expression)>0){                # FeatureHeatmap.1
                p$data <- p$data[order(p$data$scaled.expression),]
                max.scaled.exp <- max(p$data$scaled.expression)
        } else if (length(p$data$gene)>0){                   # SingleFeaturePlot.1 
                p$data <- p$data[order(p$data$gene),] 
                max.scaled.exp <- max(p$data$gene) 
        }
        
        # Define lower limit of scaled gene expression level
        if (!is.null(scaled.expression.threshold)) {
                scaled.expression.threshold <- scaled.expression.threshold
        } else if (is.null(scaled.expression.threshold)) {
                if (length(p$data$scaled.expression)>0){
                        scaled.expression.threshold <- min(p$data$scaled.expression)
                } else if (length(p$data$gene)>0) {
                        scaled.expression.threshold <- min(p$data$gene)
                }
        }
        
        # Fill points using the scaled gene expression levels
        p$layers[[1]]$mapping$fill <- p$layers[[1]]$mapping$colour
        
        # Define transparency of points
        p$layers[[1]]$mapping$alpha <- alpha.use
        
        # Change fill and colour gradient values
        p = p + guides(colour = FALSE)
        p = p + scale_colour_gradientn(colours = gradient.use, guide = F,
                                       limits = c(scaled.expression.threshold,
                                                  max.scaled.exp),
                                       na.value = "grey") +
                scale_fill_gradientn(colours = gradient.use,
                                     name = expression(atop(Scaled, expression)),
                                     limits = c(scaled.expression.threshold,
                                                max.scaled.exp),
                                     na.value = "grey") +
                scale_alpha_continuous(range = alpha.use, guide = F)
        
        # Return plot
        return(p)
}


#' Convert data frame to list
#'
#' This function will convert a data frame to a list, even if they are unequal length
#'
#' @param df
#' @export
#' @examples
#' library(GSVAdata)
#' data(brainTxDbSets)
#' brainTxDbSets_df <- list2df(brainTxDbSets)
#' genelist <- df2list(brainTxDbSets_df)
df2list <- function(df){
        list <- lapply(df, as.vector) # as.vector! not as.character
        list <- lapply(list, function(x) x[!is.na(x)])
        names(list) <- names(df)
        return(list)
}


#' DoHeatmap.1, automatically group_top by cluster, order by Time_points
DoHeatmap.1 <- function(object, marker_df, add.genes = NULL, Top_n = 10, group.order = NULL,ident.use,
                        group.label.rot =T,cex.row = 8,remove.key =T,use.scaled = T,
                        group.cex = 13,title.size = 14,...){
        if (!missing(x = marker_df)) {
                top <-  marker_df %>% 
                        group_by(cluster) %>% 
                        top_n(Top_n, avg_logFC)
        } else { 
                top <- list()
                top$gene = NULL
        }
        if(!is.null(add.genes)) {
                top_gene = c(add.genes, top$gene)
        } else top_gene <- top$gene
        return(DoHeatmap(object = object, genes.use = top_gene,
                         group.order = group.order, use.scaled = use.scaled,
                         slim.col.label = TRUE, remove.key = remove.key,cex.row = cex.row,
                         group.cex = group.cex, rotate.key = T,group.label.rot = group.label.rot,
                         title = paste("Expression heatmap of top",Top_n,
                                       "differential expression genes in",
                                       ident.use),...))
}


#' Modified Seurat::FeatureHeatmap function to increase contrast gradient
#' 
#' #https://github.com/satijalab/seurat/issues/235
#' @param object Seurat object
#' @param features.plot gene name
#' @param gradient.use Change fill and colour gradient values
#' @param scaled.expression.threshold Define lower limit of scaled gene expression level
#' @export p ggplot object
#' @example FeatureHeatmap.1(SSCs,"Gfra1")
FeatureHeatmap.1 <- function(object, features.plot, mouse.genes = TRUE,
                             group.by ="orig.ident",  alpha.use = 1,
                             pt.size = 0.5, scaled.expression.threshold = 0.1,
                             gradient.use = c("orangered", "red4"),
                             pch.use = 20) {
        
        if (mouse.genes) features.plot = MouseGenes(object =  object, unique = T,
                                                    marker.genes = features.plot)
        else features.plot = HumanGenes(object =  object, unique = T,
                                        marker.genes = features.plot)
        
        x <- FeatureHeatmap(object = object, features.plot = features.plot,
                            group.by = group.by, sep.scale = T, pt.size = pt.size, 
                            cols.use = c("gray98", "red"), pch.use = pch.use, do.return = T)
        x$data <- x$data[order(x$data$expression),]
        p <- ChangeColorScale(x, alpha.use = alpha.use,
                              scaled.expression.threshold = scaled.expression.threshold,
                              gradient.use = gradient.use)
        return(p)
        
}


#' Modified Seurat::DimPlot function to add ggrepel::geom_text_repel and ggrepel::geom_label_repel
#' A supporting function for TSNEPlot.1
#' @param object Seurat object
#' @param text.repel Adds text directly to the plot.
#' @param label.repel Draws a rectangle underneath the text, making it easier to read.
#' @param force Force of repulsion between overlapping text labels. Defaults to 1.
#' @param ... all other paramethers are the same as Seurat::DimPlot
#' @export p/p3 ggplot object
#' @example DimPlot.1(SSCs, reduction.use = "tsne")
DimPlot.1 <- function (object, reduction.use = "pca", dim.1 = 1, dim.2 = 2,
                       cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
                       cols.use = NULL, group.by = "ident", pt.shape = NULL, do.hover = FALSE,
                       data.hover = "ident", do.identify = FALSE, do.label = FALSE,
                       label.size = 4, no.legend = FALSE, coord.fixed = FALSE,
                       no.axes = FALSE, dark.theme = FALSE, plot.order = NULL,
                       cells.highlight = NULL, cols.highlight = "red", sizes.highlight = 1,
                       plot.title = NULL, vector.friendly = FALSE, png.file = NULL,
                       png.arguments = c(10, 10, 100), na.value = "grey50",
                       text.repel = TRUE, label.repel = FALSE,force=1, ...)
{
        if (vector.friendly) {
                previous_call <- blank_call <- png_call <- match.call()
                blank_call$pt.size <- -1
                blank_call$do.return <- TRUE
                blank_call$vector.friendly <- FALSE
                png_call$no.axes <- TRUE
                png_call$no.legend <- TRUE
                png_call$do.return <- TRUE
                png_call$vector.friendly <- FALSE
                png_call$plot.title <- NULL
                blank_plot <- eval(blank_call, sys.frame(sys.parent()))
                png_plot <- eval(png_call, sys.frame(sys.parent()))
                png.file <- Seurat:::SetIfNull(x = png.file, default = paste0(tempfile(),
                                                                              ".png"))
                ggsave(filename = png.file, plot = png_plot, width = png.arguments[1],
                       height = png.arguments[2], dpi = png.arguments[3])
                to_return <- AugmentPlot(plot1 = blank_plot, imgFile = png.file)
                file.remove(png.file)
                if (do.return) {
                        return(to_return)
                }
                else {
                        print(to_return)
                }
        }
        embeddings.use <- GetDimReduction(object = object, reduction.type = reduction.use,
                                          slot = "cell.embeddings")
        if (length(x = embeddings.use) == 0) {
                stop(paste(reduction.use, "has not been run for this object yet."))
        }
        cells.use <- Seurat:::SetIfNull(x = cells.use, default = colnames(x = object@data))
        dim.code <- GetDimReduction(object = object, reduction.type = reduction.use,
                                    slot = "key")
        dim.codes <- paste0(dim.code, c(dim.1, dim.2))
        data.plot <- as.data.frame(x = embeddings.use)
        cells.use <- intersect(x = cells.use, y = rownames(x = data.plot))
        data.plot <- data.plot[cells.use, dim.codes]
        ident.use <- as.factor(x = object@ident[cells.use])
        if (group.by != "ident") {
                ident.use <- as.factor(x = FetchData(object = object,
                                                     vars.all = group.by)[cells.use, 1])
        }
        data.plot$ident <- ident.use
        data.plot$x <- data.plot[, dim.codes[1]]
        data.plot$y <- data.plot[, dim.codes[2]]
        data.plot$pt.size <- pt.size
        if (!is.null(x = cells.highlight)) {
                if (is.character(x = cells.highlight)) {
                        cells.highlight <- list(cells.highlight)
                }
                else if (is.data.frame(x = cells.highlight) || !is.list(x = cells.highlight)) {
                        cells.highlight <- as.list(x = cells.highlight)
                }
                cells.highlight <- lapply(X = cells.highlight, FUN = function(cells) {
                        cells.return <- if (is.character(x = cells)) {
                                cells[cells %in% rownames(x = data.plot)]
                        }
                        else {
                                cells <- as.numeric(x = cells)
                                cells <- cells[cells <= nrow(x = data.plot)]
                                rownames(x = data.plot)[cells]
                        }
                        return(cells.return)
                })
                cells.highlight <- Filter(f = length, x = cells.highlight)
                if (length(x = cells.highlight) > 0) {
                        if (!no.legend) {
                                no.legend <- is.null(x = names(x = cells.highlight))
                        }
                        names.highlight <- if (is.null(x = names(x = cells.highlight))) {
                                paste0("Group_", 1L:length(x = cells.highlight))
                        }
                        else {
                                names(x = cells.highlight)
                        }
                        sizes.highlight <- rep_len(x = sizes.highlight,
                                                   length.out = length(x = cells.highlight))
                        cols.highlight <- rep_len(x = cols.highlight, length.out = length(x = cells.highlight))
                        highlight <- rep_len(x = NA_character_, length.out = nrow(x = data.plot))
                        if (is.null(x = cols.use)) {
                                cols.use <- "black"
                        }
                        cols.use <- c(cols.use[1], cols.highlight)
                        size <- rep_len(x = pt.size, length.out = nrow(x = data.plot))
                        for (i in 1:length(x = cells.highlight)) {
                                cells.check <- cells.highlight[[i]]
                                index.check <- match(x = cells.check, rownames(x = data.plot))
                                highlight[index.check] <- names.highlight[i]
                                size[index.check] <- sizes.highlight[i]
                        }
                        plot.order <- sort(x = unique(x = highlight), na.last = TRUE)
                        plot.order[is.na(x = plot.order)] <- "Unselected"
                        highlight[is.na(x = highlight)] <- "Unselected"
                        highlight <- as.factor(x = highlight)
                        data.plot$ident <- highlight
                        data.plot$pt.size <- size
                        if (dark.theme) {
                                cols.use[1] <- "white"
                        }
                }
        }
        if (!is.null(x = plot.order)) {
                if (any(!plot.order %in% data.plot$ident)) {
                        stop("invalid ident in plot.order")
                }
                plot.order <- rev(x = c(plot.order, setdiff(x = unique(x = data.plot$ident),
                                                            y = plot.order)))
                data.plot$ident <- factor(x = data.plot$ident, levels = plot.order)
                data.plot <- data.plot[order(data.plot$ident), ]
        }
        p <- ggplot(data = data.plot, mapping = aes(x = x, y = y, colour =ident)) +
                geom_point(mapping = aes(colour = factor(x = ident),
                                         size = pt.size))
        if (!is.null(x = pt.shape)) {
                shape.val <- FetchData(object = object, vars.all = pt.shape)[cells.use,
                                                                             1]
                if (is.numeric(shape.val)) {
                        shape.val <- cut(x = shape.val, breaks = 5)
                }
                data.plot[, "pt.shape"] <- shape.val
                p <- ggplot(data = data.plot, mapping = aes(x = x, y = y)) +
                        geom_point(mapping = aes(colour = factor(x = ident),
                                                 shape = factor(x = pt.shape), size = pt.size))
        }
        if (!is.null(x = cols.use)) {
                p <- p + scale_colour_manual(values = cols.use, na.value = na.value)
        }
        if (coord.fixed) {
                p <- p + coord_fixed()
        }
        p <- p + guides(size = FALSE)
        p2 <- p + xlab(label = dim.codes[[1]]) + ylab(label = dim.codes[[2]]) +
                scale_size(range = c(min(data.plot$pt.size), max(data.plot$pt.size)))
        p3 <- p2 + Seurat:::SetXAxisGG() + Seurat:::SetYAxisGG() +
                Seurat:::SetLegendPointsGG(x = 6) +
                Seurat:::SetLegendTextGG(x = 12) + Seurat:::no.legend.title + theme_bw() +
                Seurat:::NoGrid()
        if (dark.theme) {
                p <- p + DarkTheme()
                p3 <- p3 + DarkTheme()
        }
        p3 <- p3 + theme(legend.title = element_blank())
        if (!is.null(plot.title)) {
                p3 <- p3 + ggtitle(plot.title) + theme(plot.title = element_text(hjust = 0.5))
        }
        if (do.label) {
                centers <- data.plot %>% dplyr::group_by(ident) %>%
                        dplyr::summarize(x = median(x), y = median(y))
                p3 = p3 + geom_point(data = centers, aes(x = x, y = y),
                                     size = 0, alpha = 0)
                if (label.repel == TRUE) {
                        p3 = p3 + ggrepel::geom_label_repel(data = centers,
                                                            aes(label = ident),
                                                            size = label.size,
                                                            force = force)
                }
                else if (text.repel == TRUE){
                        p3 = p3 + ggrepel::geom_text_repel(data = centers,
                                                           aes(label = ident),
                                                           size = label.size,
                                                           force = force,
                                                           color = "black")
                }
                if (label.repel == FALSE & text.repel == FALSE) {
                        p3 = p3 + geom_text(data = centers,
                                            aes(label = ident),
                                            size = label.size,
                                            color = "black")
                }
                p3 = p3 + guides(colour = FALSE)
                x.range = layer_scales(p)$x$range$range
                add_to_x = sum(abs(x.range)) * 0.03
                p3 = p3 + xlim(x.range[1] - add_to_x, x.range[2] + add_to_x)
        }
        if (no.legend) {
                p3 <- p3 + theme(legend.position = "none")
        }
        if (no.axes) {
                p3 <- p3 + theme(axis.line = element_blank(), axis.text.x = element_blank(),
                                 axis.text.y = element_blank(), axis.ticks = element_blank(),
                                 axis.title.x = element_blank(), axis.title.y = element_blank(),
                                 panel.background = element_blank(), panel.border = element_blank(),
                                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 plot.background = element_blank())
        }
        if (do.identify || do.hover) {
                if (do.bare) {
                        plot.use <- p
                }
                else {
                        plot.use <- p3
                }
                if (do.hover) {
                        if (is.null(x = data.hover)) {
                                features.info <- NULL
                        }
                        else {
                                features.info <- FetchData(object = object,
                                                           vars.all = data.hover)
                        }
                        return(HoverLocator(plot = plot.use, data.plot = data.plot,
                                            features.info = features.info, dark.theme = dark.theme))
                }
                else if (do.identify) {
                        return(FeatureLocator(plot = plot.use, data.plot = data.plot,
                                              dark.theme = dark.theme, ...))
                }
        }
        if (do.return) {
                if (do.bare) {
                        return(p)
                }
                else {
                        return(p3)
                }
        }
        if (do.bare) {
                print(p)
        }
        else {
                print(p3)
        }
}


#' Combine FindAllMarkers and calculate average UMI
#' Modified Seurat::FindAllMarkers function to add average UMI for group1 (UMI.1) and group 2 (UMI.2)
#' @param ... all paramethers are the same as Seurat::FindAllMarkers
#' @param get.slot option to get different data slot like get.slot = "scale.data"
#' #' @export gde.all data frame
#' @example FindAllMarkers.UMI(SSCs)
FindAllMarkers.UMI <- function (object, genes.use = NULL, logfc.threshold = 0.25, 
                                test.use = "MAST", min.pct = 0.1, min.diff.pct = -Inf, 
                                print.bar = TRUE, only.pos = FALSE, max.cells.per.ident = Inf, 
                                return.thresh = 0.01, do.print = FALSE, random.seed = 1, 
                                min.cells = 3, latent.vars = "nUMI", assay.type = "RNA", 
                                get.slot = "data",...)
{
        data.1 <- GetAssayData(object = object, assay.type = assay.type, 
                               slot = get.slot)
        genes.use <- Seurat:::SetIfNull(x = genes.use, default = rownames(x = data.1))
        if ((test.use == "roc") && (return.thresh == 0.01)) {
                return.thresh = 0.7
        }
        idents.all <- sort(x = unique(x = object@ident))
        genes.de <- list()
        avg_UMI <- list()
        for (i in 1:length(x = idents.all)) {
                genes.de[[i]] <- tryCatch({
                        FindMarkers.UMI(object = object, assay.type = assay.type, 
                                        ident.1 = idents.all[i], ident.2 = NULL, genes.use = genes.use, 
                                        logfc.threshold = logfc.threshold, test.use = test.use, 
                                        min.pct = min.pct, min.diff.pct = min.diff.pct, 
                                        print.bar = print.bar, min.cells = min.cells, 
                                        latent.vars = latent.vars, max.cells.per.ident = max.cells.per.ident, 
                                        get.slot = get.slot, ...)
                })
                if (do.print) {
                        print(paste("Calculating cluster", idents.all[i]))
                }
        }
        gde.all <- data.frame()
        for (i in 1:length(x = idents.all)) {
                if (is.null(x = unlist(x = genes.de[i]))) {
                        next
                }
                gde <- genes.de[[i]]
                if (nrow(x = gde) > 0) {
                        if (test.use == "roc") {
                                gde <- subset(x = gde, subset = (myAUC > return.thresh | 
                                                                         myAUC < (1 - return.thresh)))
                        }
                        else {
                                gde <- gde[order(gde$p_val, -gde$avg_logFC), 
                                           ]
                                gde <- subset(x = gde, subset = p_val < return.thresh)
                        }
                        if (nrow(x = gde) > 0) {
                                gde$cluster <- idents.all[i]
                                gde$gene <- rownames(x = gde)
                        }
                        if (nrow(x = gde) > 0) {
                                gde.all <- rbind(gde.all, gde)
                        }
                }
        }
        if (only.pos) {
                return(subset(x = gde.all, subset = avg_logFC > 0))
        }
        rownames(x = gde.all) <- make.unique(names = as.character(x = gde.all$gene))
        return(gde.all)
}


#' Calculate average UMI and attach to FindMarkers results 
#' Modified Seurat::FindMarkers function to add average UMI for group1 (UMI.1) and group 2 (UMI.2)
#' @param ... all paramethers are the same as Seurat::FindMarkers
#' @export gde.all data frame
#' @example FindMarkers.UMI(SSCs,ident.1 = "1")
FindMarkers.UMI <- function (object, ident.1, ident.2 = NULL, genes.use = NULL,
                             logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1,
                             min.diff.pct = -Inf, print.bar = TRUE, only.pos = FALSE,
                             max.cells.per.ident = Inf, random.seed = 1, latent.vars = "nUMI",
                             min.cells.gene = 3, pseudocount.use = 1, assay.type = "RNA",
                             get.slot = "data",...)
{
        genes.de <- FindMarkers.1(object = object, assay.type = assay.type,
                                  ident.1 = ident.1, ident.2 = ident.2, genes.use = genes.use,
                                  logfc.threshold = logfc.threshold, test.use = test.use,
                                  min.pct = min.pct, min.diff.pct = min.diff.pct,
                                  print.bar = print.bar, min.cells.gene = min.cells.gene,
                                  latent.vars = latent.vars, max.cells.per.ident = max.cells.per.ident,
                                  get.slot = get.slot, ...)
        UMI.1 <-rowMeans(as.matrix(x = object@raw.data[, WhichCells(object = object,
                                                                ident = ident.1)]))
        UMI.2 <-rowMeans(as.matrix(x = object@raw.data[, WhichCells(object = object,
                                                                ident.remove = ident.1)]))
        avg_UMI <-data.frame(UMI.1, UMI.2)
        genes.de <- cbind(genes.de,avg_UMI[match(rownames(genes.de),rownames(avg_UMI)),])
        
        return(genes.de)
}

#' Enable FindMarkers to get different data slot
FindMarkers.1 <- function(object, ident.1, ident.2 = NULL, genes.use = NULL, 
                          logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                          min.diff.pct = -Inf, print.bar = TRUE, only.pos = FALSE, 
                          max.cells.per.ident = Inf, random.seed = 1, latent.vars = NULL, 
                          min.cells.gene = 3, min.cells.group = 3, pseudocount.use = 1, 
                          assay.type = "RNA", get.slot = "data", ...) 
{
        data.use <- GetAssayData(object = object, assay.type = assay.type, 
                                 slot = get.slot)
        genes.use <- Seurat:::SetIfNull(x = genes.use, default = rownames(x = data.use))
        methods.noprefiliter <- c("DESeq2", "zingeR")
        if (test.use %in% methods.noprefiliter) {
                genes.use <- rownames(x = data.use)
                min.diff.pct <- -Inf
                logfc.threshold <- 0
        }
        if (length(x = as.vector(x = ident.1) > 1) && any(as.character(x = ident.1) %in% 
                                                          object@cell.names)) {
                cells.1 <- intersect(x = ident.1, y = object@cell.names)
        }
        else {
                cells.1 <- WhichCells(object = object, ident = ident.1)
        }
        if (length(x = as.vector(x = ident.2) > 1) && any(as.character(x = ident.2) %in% 
                                                          object@cell.names)) {
                cells.2 <- intersect(x = ident.2, y = object@cell.names)
        }
        else {
                if (is.null(x = ident.2)) {
                        cells.2 <- WhichCells(object = object, cells.use = setdiff(object@cell.names, 
                                                                                   cells.1))
                }
                else {
                        cells.2 <- WhichCells(object = object, ident = ident.2)
                }
        }
        cells.2 <- setdiff(x = cells.2, y = cells.1)
        if (length(x = cells.1) == 0) {
                message(paste("Cell group 1 is empty - no cells with identity class", 
                              ident.1))
                return(NULL)
        }
        if (length(x = cells.2) == 0) {
                message(paste("Cell group 2 is empty - no cells with identity class", 
                              ident.2))
                return(NULL)
        }
        if (length(cells.1) < min.cells.group) {
                stop(paste("Cell group 1 has fewer than", as.character(min.cells.group), 
                           "cells in identity class", ident.1))
        }
        if (length(cells.2) < min.cells.group) {
                stop(paste("Cell group 2 has fewer than", as.character(min.cells.group), 
                           " cells in identity class", ident.2))
        }
        thresh.min <- 0
        data.temp1 <- round(x = apply(X = data.use[genes.use, cells.1, 
                                                   drop = F], MARGIN = 1, FUN = function(x) {
                                                           return(sum(x > thresh.min)/length(x = x))
                                                   }), digits = 3)
        data.temp2 <- round(x = apply(X = data.use[genes.use, cells.2, 
                                                   drop = F], MARGIN = 1, FUN = function(x) {
                                                           return(sum(x > thresh.min)/length(x = x))
                                                   }), digits = 3)
        data.alpha <- cbind(data.temp1, data.temp2)
        colnames(x = data.alpha) <- c("pct.1", "pct.2")
        alpha.min <- apply(X = data.alpha, MARGIN = 1, FUN = max)
        names(x = alpha.min) <- rownames(x = data.alpha)
        genes.use <- names(x = which(x = alpha.min > min.pct))
        if (length(x = genes.use) == 0) {
                stop("No genes pass min.pct threshold")
        }
        alpha.diff <- alpha.min - apply(X = data.alpha, MARGIN = 1, 
                                        FUN = min)
        genes.use <- names(x = which(x = alpha.min > min.pct & alpha.diff > 
                                             min.diff.pct))
        if (length(x = genes.use) == 0) {
                stop("No genes pass min.diff.pct threshold")
        }
        data.1 <- apply(X = data.use[genes.use, cells.1, drop = F], 
                        MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) + 
                                                                  pseudocount.use))
        data.2 <- apply(X = data.use[genes.use, cells.2, drop = F], 
                        MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) + 
                                                                  pseudocount.use))
        total.diff <- (data.1 - data.2)
        if (!only.pos) 
                genes.diff <- names(x = which(x = abs(x = total.diff) > 
                                                      logfc.threshold))
        if (only.pos) 
                genes.diff <- names(x = which(x = total.diff > logfc.threshold))
        genes.use <- intersect(x = genes.use, y = genes.diff)
        if (length(x = genes.use) == 0) {
                stop("No genes pass logfc.threshold threshold")
        }
        if (max.cells.per.ident < Inf) {
                set.seed(seed = random.seed)
                if (length(cells.1) > max.cells.per.ident) 
                        cells.1 = sample(x = cells.1, size = max.cells.per.ident)
                if (length(cells.2) > max.cells.per.ident) 
                        cells.2 = sample(x = cells.2, size = max.cells.per.ident)
        }
        if (!(test.use %in% c("negbinom", "poisson", "MAST")) && 
            !is.null(x = latent.vars)) {
                warning("'latent.vars' is only used for 'negbinom', 'poisson', and 'MAST' tests")
        }
        if (test.use == "bimod") {
                to.return <- DiffExpTest(object = object, assay.type = assay.type, 
                                         cells.1 = cells.1, cells.2 = cells.2, genes.use = genes.use, 
                                         print.bar = print.bar)
        }
        if (test.use == "roc") {
                to.return <- MarkerTest(object = object, assay.type = assay.type, 
                                        cells.1 = cells.1, cells.2 = cells.2, genes.use = genes.use, 
                                        print.bar = print.bar)
        }
        if (test.use == "t") {
                to.return <- DiffTTest(object = object, assay.type = assay.type, 
                                       cells.1 = cells.1, cells.2 = cells.2, genes.use = genes.use, 
                                       print.bar = print.bar)
        }
        if (test.use == "tobit") {
                to.return <- TobitTest(object = object, assay.type = assay.type, 
                                       cells.1 = cells.1, cells.2 = cells.2, genes.use = genes.use, 
                                       print.bar = print.bar)
        }
        if (test.use == "negbinom") {
                to.return <- NegBinomDETest(object = object, assay.type = assay.type, 
                                            cells.1 = cells.1, cells.2 = cells.2, genes.use = genes.use, 
                                            latent.vars = latent.vars, print.bar = print.bar, 
                                            min.cells = min.cells.gene)
        }
        if (test.use == "poisson") {
                to.return <- PoissonDETest(object = object, assay.type = assay.type, 
                                           cells.1 = cells.1, cells.2 = cells.2, genes.use = genes.use, 
                                           latent.vars = latent.vars, print.bar = print.bar, 
                                           min.cells = min.cells.gene)
        }
        if (test.use == "MAST") {
                to.return <- MASTDETest(object = object, assay.type = assay.type, 
                                        cells.1 = cells.1, cells.2 = cells.2, genes.use = genes.use, 
                                        latent.vars = latent.vars, ...)
        }
        if (test.use == "wilcox") {
                to.return <- WilcoxDETest(object = object, assay.type = assay.type, 
                                          cells.1 = cells.1, cells.2 = cells.2, genes.use = genes.use, 
                                          print.bar = print.bar, ...)
        }
        if (test.use == "LR") {
                to.return <- LRDETest(object = object, assay.type = assay.type, 
                                      cells.1 = cells.1, cells.2 = cells.2, genes.use = genes.use, 
                                      print.bar = print.bar, ...)
        }
        if (test.use == "DESeq2") {
                to.return <- DESeq2DETest(object = object, assay.type = assay.type, 
                                          cells.1 = cells.1, cells.2 = cells.2, genes.use = genes.use, 
                                          ...)
        }
        to.return[, "avg_logFC"] <- total.diff[rownames(x = to.return)]
        to.return <- cbind(to.return, data.alpha[rownames(x = to.return), 
                                                 , drop = FALSE])
        to.return$p_val_adj = p.adjust(p = to.return$p_val, method = "bonferroni", 
                                       n = nrow(x = GetAssayData(object = object, assay.type = assay.type, 
                                                                 slot = "data")))
        if (test.use == "roc") {
                to.return <- to.return[order(-to.return$power, -to.return$avg_logFC), 
                                       ]
        }
        else {
                to.return <- to.return[order(to.return$p_val, -to.return$avg_logFC), 
                                       ]
        }
        if (only.pos) {
                to.return <- subset(x = to.return, subset = avg_logFC > 
                                            0)
        }
        return(to.return)
}


FPKM <- function(counts, lengths) {
        rownames(counts) = tolower(rownames(counts))
        names(lengths) = tolower(names(lengths))
        A = intersect(rownames(counts), names(lengths))
        counts = counts[A, ]
        lengths = lengths[A]
        rate = counts/lengths
        sweep(rate,2,FUN="/",STATS=colSums(counts)) * 1e+06
}

#=====Clean memory======================
GC <- function()
{
        while (gc()[2, 4] != gc()[2, 4] | gc()[1, 4] != gc()[1,
                                                             4]) {
        }
}

#' Extract ColorHexa from Seurat TSNE plot
#' @param object aligned seurat object with ident
#' @param return.vector TRUE/return full color vector, FALSE/return color levels
#' @param cells.use only return ColorHexa of selected cells
#' @param ... other TSNEPlot inputs 
#' @export colors: color vector named by cell ID
gg_colors <- function(object = object, return.vector=FALSE, cells.use = NULL,
                      no.legend = TRUE, do.label = TRUE,
                      do.return = TRUE, label.size = 6, gg_title="", ...){
        
        g1 <- Seurat::TSNEPlot(object = object, no.legend = no.legend,
                               do.label = do.label,do.return = do.return,
                               label.size = label.size,...)+
                ggtitle(gg_title)+
                theme(text = element_text(size=15),     #larger text including legend title							
                      plot.title = element_text(hjust = 0.5)) #title in middle
        print(g1)
        g <- ggplot2::ggplot_build(g1)
        #        print(unique(g$data[[1]]["colour"]))
        colors <- unlist(g$data[[1]]["colour"])
        
        #select color by cells ID
        cells <- Seurat::WhichCells(object)
        names(colors) <- cells
        if(!is.null(cells.use)) colors <- colors[cells.use]
        
        if(return.vector) {
                print(head(colors));print(length(colors))
                return(colors)
        } else {
                colors <- as.data.frame(table(colors))
                colors <- colors$colors
                print(colors)
                return(as.character(colors))
        }
}



# select ggplot color
gg_color_hue <- function(n) {
        hues = seq(15, 375, length = n + 1)
        grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}
#gg_color_hue(4)


#' group_by + top_n + mutate + re-arrange data frame
#' @param df  a data frame from FindAllMarkers
#' @param ... Name-value pairs of expressions, same as ... in dplyr::mutate
#' @param Top_n number of rows to return, same as n in dplyr::top_n. Default is NULL, return all rows.
#' @export top a re-arranged data frame, sorted by avg_logFC, then arrange by ...
#' @example 
#' major_cells <- c("Spermatogonia","Early Spermatocytes","Spermatocytes","Spermatids")
#' top <- group_top(df = All_markers, major_cells)
group_top_mutate <- function(df, ..., Top_n = 500){
        rownames(df) = NULL
        df <- df %>% dplyr::select("gene", everything()) # Moving the last column to the start
        new.col = deparse(substitute(...))
        new.order = assign(new.col,...)
        if(class(df) != "data.frame") df = as.data.frame(df)
        top <-  df %>% 
                dplyr::select("gene", everything()) %>%
                group_by(cluster) %>% 
                top_n(Top_n, wt = avg_logFC) %>%
                mutate(new.col = factor(cluster, levels = new.order)) %>%
                arrange(new.col)
        colnames(top)[which(colnames(top) == "new.col")] = new.col
        return(as.data.frame(top))
}


# modified GenePlot
# GenePlot.1(do.hover = TRUE) will return ggplot 
GenePlot.1 <- function (object, gene1, gene2, cell.ids = NULL, col.use = NULL, 
                        pch.use = 16, pt.size = 2, use.imputed = FALSE, use.scaled = FALSE, 
                        use.raw = FALSE, do.hover = TRUE, data.hover = "ident",  no.legend = TRUE,
                        do.identify = FALSE, dark.theme = FALSE, do.spline = FALSE, 
                        spline.span = 0.75, ...) 
{
        cell.ids <- Seurat:::SetIfNull(x = cell.ids, default = object@cell.names)
        data.use <- as.data.frame(x = FetchData(object = object, 
                                                vars.all = c(gene1, gene2), cells.use = cell.ids, use.imputed = use.imputed, 
                                                use.scaled = use.scaled, use.raw = use.raw))
        data.plot <- data.use[cell.ids, c(gene1, gene2)]
        names(x = data.plot) <- c("x", "y")
        ident.use <- as.factor(x = object@ident[cell.ids])
        if (length(x = col.use) > 1) {
                col.use <- col.use[as.numeric(x = ident.use)]
        }
        else {
                col.use <- Seurat:::SetIfNull(x = col.use, default = as.numeric(x = ident.use))
        }
        gene.cor <- round(x = cor(x = data.plot$x, y = data.plot$y), 
                          digits = 2)
        if (dark.theme) {
                par(bg = "black")
                col.use <- sapply(X = col.use, FUN = function(color) ifelse(test = all(col2rgb(color) == 
                                                                                               0), yes = "white", no = color))
                axes = FALSE
                col.lab = "white"
        }
        else {
                axes = TRUE
                col.lab = "black"
        }
        
        if (dark.theme) {
                axis(side = 1, at = NULL, labels = TRUE, col.axis = col.lab, 
                     col = col.lab)
                axis(side = 2, at = NULL, labels = TRUE, col.axis = col.lab, 
                     col = col.lab)
        }
        if (do.spline) {
                spline.fit <- smooth.spline(x = data.plot$x, y = data.plot$y, 
                                            df = 4)
                loess.fit <- loess(formula = y ~ x, data = data.plot, 
                                   span = spline.span)
                points(x = data.plot$x, y = loess.fit$fitted, col = "darkblue")
        }
        if (do.identify | do.hover) {
                p <- ggplot2::ggplot(data = data.plot, mapping = aes(x = x, 
                                                                     y = y))
                p <- p + geom_point(mapping = aes(x = x,y=y, color = col.use), size = pt.size, 
                                    shape = pch.use )
                p <- p + labs(title = gene.cor, x = gene1, y = gene2)
                if (no.legend) {
                        p <- p + theme(legend.position = "none")
                }
                return(p)
        }
}

# HumanGenes
# turn list of gene character to uniformed Human gene list format
HumanGenes <- function(object, marker.genes, unique=FALSE){
        # marker.genes =c("Cdh5,Pecam1,Flt1,Vwf,Plvap,Kdr") for example
        if(missing(object)) 
                stop("A seurat object must be provided first")
        if(class(object) != "seurat") 
                stop("A seurat object must be provided first")
        if(missing(marker.genes)) 
                stop("A list of marker genes must be provided")
        if(object@var.genes[1]==Hmisc::capitalize(tolower(object@var.genes[1])))
                stop("This is mouse genome, use MouseGenes() instead!")
        
        marker.genes <- as.character(marker.genes)
        marker.genes <- unlist(strsplit(marker.genes,","))
        marker.genes <- toupper(marker.genes)
        print(paste("Before filtration:",length(marker.genes)))
        marker.genes <- CaseMatch(search = marker.genes, match = rownames(x = object@raw.data))
        if(unique) marker.genes <- unique(marker.genes)
        print(paste("After filtration:",length(marker.genes)))
        return(as.character(marker.genes))
}


#' Convert list to data frame
#'
#' This function will convert a list to a data frame, even if they are unequal length
#'
#' @param fname list
#' @export
#' @examples
#' library(GSVAdata)
#' data(brainTxDbSets)
#' brainTxDbSets_df <- list2df(brainTxDbSets)
list2df <- function(list){
        df <- do.call(rowr::cbind.fill, c(list, fill = NA))
        names(df) = names(list)
        return(df)
}


# clean up the gene names for downstream analysis
# turn list of gene character to uniform mouse gene list format
#' @param object Seurat object
#' @param marker.genes gene names, can be one gene or vector. Must be character
#' @param unique TRUE/FALSE, output unique gene name or not
#' @export marker.genes uniform mouse gene that exsit in object@data
#' @example MouseGenes(mouse_eyes,c("Cdh5","Pecam1","Flt1"))
MouseGenes <- function(object, marker.genes, unique =F){
        # marker.genes =c("Cdh5,Pecam1,Flt1,Vwf,Plvap,Kdr") for example
        if(missing(object))
                stop("A seurat object must be provided first!")
        if(class(object) != "seurat")
                stop("A seurat object must be provided first!")
        if(missing(marker.genes))
                stop("A list of marker genes must be provided!")
        if(object@var.genes[1]==toupper(object@var.genes[1]))
                stop("This is human genome, use HumanGenes() instead!")
        
        marker.genes <- as.character(marker.genes)
        marker.genes <- unlist(strsplit(marker.genes,","))
        marker.genes <- Hmisc::capitalize(tolower(marker.genes))
        print(paste("Before filtration:",length(marker.genes)))
        marker.genes <- CaseMatch(search = marker.genes, match = rownames(x = object@raw.data))
        if(unique) marker.genes <- unique(marker.genes)
        print(paste("After filtration:",length(marker.genes)))
        return(as.character(marker.genes))
}


randomStrings <- function(n = 5000) {
        a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
        paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}


#' rename ident for certain cells only
#' One drawback, will not double check the duplicated cell names
#' @param object Seurat object
#' @param new.ident.name one character
#' @param cells.use a vector of cell names
#' @export object Seurat object with new ident
#' @example 
#' for(i in 1:length(SSC_labels_id)){
#'         SSCs <- RenameIdent.1(SSCs, new.ident.name = names(SSC_labels_id)[i],
#'                               cells.use = SSC_labels_id[[i]])
#' }
RenameIdent.1 <- function (object, new.ident.name, cells.use) 
{
        new.levels <- old.levels <- levels(x = object@ident)
        if (!(new.ident.name %in% old.levels)) {
                new.levels <- c(new.levels, new.ident.name)
        }
        
        ident.vector <- as.character(x = object@ident)
        names(x = ident.vector) <- names(object@ident)
        ident.vector[cells.use] <- new.ident.name
        object@ident <- factor(x = ident.vector, levels = new.levels)
        
        return(object)
}


SetAllIdent.1 <- function (object, id = NULL) 
{
        id <- Seurat:::SetIfNull(x = id, default = "orig.ident")
        if (id %in% colnames(x = object@meta.data)) {
                cells.use <- rownames(x = object@meta.data)
                ident.use <- factor(x = object@meta.data[, id])
                names(ident.use) = cells.use
                object@ident <- ident.use
        }
        return(object)
}


# FeaturePlot doesn't return ggplot
# SingleFeaturePlot doesn't take seurat object as input
# modified SingleFeaturePlot, take seurat object and return ggplot
SingleFeaturePlot.1 <- function (object = object, feature = feature, pt.size = 1.0,
                                 dim.1 = 1, dim.2 = 2,pch.use = 16, cols.use  = c("lightgrey","blue"),
                                 gradient.use = c("orangered", "red4"),threshold=0.1,text.size=15,
                                 cells.use = NULL,dim.codes, min.cutoff = 0, max.cutoff = Inf,
                                 use.imputed = FALSE, reduction.use = "tsne",no.axes = FALSE, no.legend = TRUE, 
                                 dark.theme = FALSE, do.return = FALSE,...) 
{
        dim.code <- GetDimReduction(object = object, reduction.type = reduction.use, 
                                    slot = "key")
        dim.codes <- paste0(dim.code, c(dim.1, dim.2))
        data.plot <- as.data.frame(GetCellEmbeddings(object = object,
                                                     reduction.type = reduction.use, 
                                                     dims.use = c(dim.1,dim.2), 
                                                     cells.use = cells.use))
        x1 <- paste0(dim.code, dim.1)
        x2 <- paste0(dim.code, dim.2)
        data.plot$x <- data.plot[, x1]
        data.plot$y <- data.plot[, x2]
        data.plot$pt.size <- pt.size
        names(x = data.plot) <- c("x", "y")
        data.use <- t(x = FetchData(object = object, vars.all = feature, 
                                    cells.use = cells.use, use.imputed = use.imputed,...))
        data.gene <- na.omit(object = data.frame(data.use[1,])) # Error in data.use[feature, ] : subscript out of bounds
        min.cutoff <- Seurat:::SetQuantile(cutoff = min.cutoff, data = data.gene)
        max.cutoff <- Seurat:::SetQuantile(cutoff = max.cutoff, data = data.gene)
        data.gene <- sapply(X = data.gene, FUN = function(x) {
                return(ifelse(test = x < min.cutoff, yes = min.cutoff, 
                              no = x))
        })
        data.gene <- sapply(X = data.gene, FUN = function(x) {
                return(ifelse(test = x > max.cutoff, yes = max.cutoff, 
                              no = x))
        })
        data.plot$gene <- data.gene
        
        if (length(x = cols.use) == 1) {
                brewer.gran <- brewer.pal.info[cols.use, ]$maxcolors
        }
        else {
                brewer.gran <- length(x = cols.use)
        }
        if (all(data.gene == 0)) {
                data.cut <- 0
        }
        else {
                data.cut <- as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.gene), 
                                                             breaks = brewer.gran)))
        }
        data.plot$col <- as.factor(x = data.cut)
        p <- ggplot(data = data.plot, mapping = aes(x = x, y = y))
        if (brewer.gran != 2) {
                if (length(x = cols.use) == 1) {
                        p <- p + geom_point(mapping = aes(color = col), 
                                            size = pt.size, shape = pch.use)
                }
                else {
                        p <- p + geom_point(mapping = aes(color = col), 
                                            size = pt.size, shape = pch.use)
                }
        }
        else {
                if (all(data.plot$gene == data.plot$gene[1])) {
                        warning(paste0("All cells have the same value of ", 
                                       feature, "."))
                        p <- p + geom_point(color = cols.use[1], size = pt.size, 
                                            shape = pch.use)
                }
                else {
                        p <- p + geom_point(mapping = aes(color = gene), 
                                            size = pt.size, shape = pch.use)
                }
        }
        if (no.axes) {
                p <- p + labs(title = feature, x = "", y = "") + theme(axis.line = element_blank(), 
                                                                       axis.text.x = element_blank(), 
                                                                       axis.text.y = element_blank(), 
                                                                       axis.ticks = element_blank(), 
                                                                       axis.title.x = element_blank(), 
                                                                       axis.title.y = element_blank())
        }
        else {
                p <- p + labs(x = dim.codes[1], y = dim.codes[2])
        }
        if (no.legend) {
                p <- p + theme(legend.position = "none")
        }
        if (dark.theme) {
                p <- p + DarkTheme()
        }
        p$data <- p$data[order(p$data$gene),]
        p1 <- ChangeColorScale(p, alpha.use = 1,
                               scaled.expression.threshold = threshold,
                               gradient.use = gradient.use)
        p1 <- p1 +ggtitle(feature)+
                theme(text = element_text(size=text.size),     #larger text including legend title							
                      axis.text.x = element_text(size=text.size*0.8),
                      axis.text.y = element_text(size=text.size*0.8),
                      plot.title = element_text(hjust = 0.5,size=text.size*1.5)) #title in middle
        return(p1)
}


# add x.lab
SingleVlnPlot.1 <- function (feature, data, cell.ident, do.sort, y.max, size.x.use, 
          size.y.use, size.title.use, adjust.use, point.size.use, x.lab,
          cols.use, gene.names, y.log, x.lab.rot, y.lab.rot, legend.position, 
          remove.legend) 
{
        feature.name <- colnames(data)
        colnames(data) <- "feature"
        feature <- "feature"
        set.seed(seed = 42)
        data$ident <- cell.ident
        if (do.sort) {
                data$ident <- factor(x = data$ident, 
                                     levels = names(x = rev(x = sort(x = tapply(X = data[,feature],
                                                                                INDEX = data$ident,
                                                                                FUN = mean)))))
        }
        if (y.log) {
                noise <- rnorm(n = length(x = data[, feature]))/200
                data[, feature] <- data[, feature] + 1
        }
        else {
                noise <- rnorm(n = length(x = data[, feature]))/1e+05
        }
        if (all(data[, feature] == data[, feature][1])) {
                warning(paste0("All cells have the same value of ", 
                               feature, "."))
        }
        else {
                data[, feature] <- data[, feature] + noise
        }
        y.max <- Seurat:::SetIfNull(x = y.max, default = max(data[, feature]))
        plot <- ggplot(data = data, mapping = aes(x = factor(x = ident), 
                                                  y = feature)) + 
                geom_violin(scale = "width", adjust = adjust.use,
                            trim = TRUE, mapping = aes(fill = factor(x = ident))) + 
                guides(fill = guide_legend(title = NULL)) + xlab(x.lab) + 
                Seurat:::NoGrid() + ggtitle(feature) + theme(plot.title = element_text(size = size.title.use, 
                                                                              face = "bold"), legend.position = legend.position, axis.title.x = element_text(face = "bold", 
                                                                                                                                                             colour = "#990000", size = size.x.use), axis.title.y = element_text(face = "bold", 
                                                                                                                                                                                                                                 colour = "#990000", size = size.y.use))
        if (point.size.use != 0) {
                plot <- plot + geom_jitter(height = 0, size = point.size.use)
        }
        plot <- plot + ggtitle(feature.name)
        if (y.log) {
                plot <- plot + scale_y_log10()
        }
        else {
                plot <- plot + ylim(min(data[, feature]), y.max)
        }
        if (feature %in% gene.names) {
                if (y.log) {
                        plot <- plot + ylab(label = "Log Expression level")
                }
                else {
                        plot <- plot + ylab(label = "Expression level")
                }
        }
        else {
                plot <- plot + ylab(label = "")
        }
        if (!is.null(x = cols.use)) {
                plot <- plot + scale_fill_manual(values = cols.use)
        }
        if (x.lab.rot) {
                plot <- plot + theme(axis.text.x = element_text(angle = 45, 
                                                                hjust = 1, size = size.x.use))
        }
        if (y.lab.rot) {
                plot <- plot + theme(axis.text.x = element_text(angle = 90, 
                                                                size = size.y.use))
        }
        if (remove.legend) {
                plot <- plot + theme(legend.position = "none")
        }
        return(plot)
}


#' Split seurat cell names by certein criteria
#' A supporting funtion to SplitSeurat
#' @param object Seurat object
#' @param split.by the criteria to split, can be gene name, or any variable in meta.data
#' @export cell.subsets list of subseted cell names by certein conditions, plus levels to split
#' @example SplitCells(mouse_eyes, split.by = "conditions")
SplitCells <- function(object = mouse_eyes, split.by = "conditions"){
        
        cell.all <- FetchData(object = object, vars.all = split.by)
        if(class(cell.all[,1]) == "numeric"){
                Levels = paste(c("Express no","Express"), split.by)
                cell.subsets[[1]] <- rownames(cell.all)[cell.all[,split.by] == 0]
                cell.subsets[[2]] <- rownames(cell.all)[cell.all[,split.by] > 0]
                cell.subsets[[3]] <- Levels # record conditions in the last return
        }
        if(class(cell.all[,1]) == "factor"){
                Levels <- levels(cell.all[,1])
                cell.subsets <- lapply(Levels, function(x)
                        rownames(cell.all)[cell.all[,split.by] == x])
                cell.subsets[[length(cell.subsets)+1]] <- Levels # record conditions in the last return
        }
        return(cell.subsets)
}


#' Split Seurat object by certein criteria
#' A supporting funtion to SplitTSNEPlot
#' @param object Seurat object
#' @param split.by the criteria to split, can be gene name, or any variable in meta.data
#' @export object.subsets list of subseted object by certein conditions,
#' plus levels to be splited
#' @example SplitCells(mouse_eyes, split.by = "conditions")
SplitSeurat <- function(object = object, split.by = "conditions"){
        
        cell.subsets <- SplitCells(object = object, split.by = split.by)
        
        object.subsets <- list()
        for(i in 1:(length(cell.subsets)-1)){
                object.subsets[[i]] <- SubsetData(object, cells.use =cell.subsets[[i]])
        }
        object.subsets[[i+1]] <- cell.subsets[[i+1]] # record conditions in the last return
        return(object.subsets)
}


#' Split Seurat by certein criteria and make tsne plot
#' @param object Seurat object
#' @param split.by the criteria to split, can be gene name, or any variable in meta.data
#' @param select.plots output order, default to NULL. If want to change,use c(2,1) for example
#' @param return.data TRUE/FASLE, return splited ojbect or not.
#' @param ... all other parameters are inherited from Seurat::TSNEPlot()
#' @export p ggplot object from TSNEplot
#' @example SplitTSNEPlot(mouse_eyes, split.by = "conditions")
#' @example SplitTSNEPlot(mouse_eyes, split.by = "Rlbp1", select.plots = c(2,1))
SplitTSNEPlot <- function(object = mouse_eyes, split.by = "conditions",
                          select.plots = NULL, return.plots = FALSE,
                          do.label = T, group.by = "ident", no.legend = TRUE,
                          pt.size = 1,label.size = 5,... ){
        
        object.subsets <- SplitSeurat(object = object, split.by = split.by)
        levels <- object.subsets[[length(object.subsets)]]
        
        p <- list()
        if(is.null(select.plots)) select.plots <- 1:length(levels)
        for(i in 1:length(select.plots)){
                p[[i]] <- TSNEPlot.1(object = object.subsets[[select.plots[i]]],
                                     do.label = do.label, group.by = group.by,
                                     do.return = T, no.legend = no.legend,
                                     pt.size = pt.size,label.size = label.size,...)+
                        ggtitle(levels[select.plots[i]])+
                        theme(text = element_text(size=20),     #larger text including legend title
                              plot.title = element_text(hjust = 0.5)) #title in middle
        }
        p <- p[lapply(p,length)>0] # remove NULL element
        if(return.plots) return(p) else print(do.call(cowplot::plot_grid, p))
}


# find and print differentially expressed genes within all major cell types
# combine SubsetData, FindAllMarkers, write.csv
#' @param object Seurat object
#' @param split.by compatible to vars.all in Seurat::FetchData. If split.by is gene name,
#' @param ... all paramethers are the same as Seurat::FindAllMarkers
#' @export gde.all data frame
#' @example FindAllMarkers.UMI(mouse_eyes)
SplitFindAllMarkers <- function(object, split.by = "conditions",test.use = "bimod",
                                write.csv = TRUE,...){
        
        object.subsets <- SplitCells(object = object, split.by = split.by)
        conditions <- object.subsets[[length(object.subsets)]] # levels of conditions
        
        object.markers <- list()
        
        for(i in 1:length(conditions)){
                object.markers[[i]] <- FindAllMarkers.UMI(object = object.subsets[[i]],
                                                          test.use = test.use,...)
                if(write.csv) write.csv(object.markers[[i]],
                                        file=paste0("./output/",
                                                    deparse(substitute(object)),
                                                    "_",conditions[i],
                                                    ".csv"))
        }
        return(object.markers)
}


SplitDotPlotGG.1 <- function (object, grouping.var, genes.plot,
                              gene.groups, cols.use = c("blue","red"),
                              col.min = -2.5, col.max = 2.5, 
                              dot.min = 0, dot.scale = 6,
                              group.by, plot.legend = FALSE,
                              do.return = FALSE, x.lab.rot = FALSE) 
{
        "Fix the mutate_impl error in orignal function"
        if (!missing(x = group.by)) {
                object <- SetAllIdent(object = object, id = group.by)
        }
        grouping.data <- FetchData(object = object, vars.all = grouping.var)[names(x = object@ident), 
                                                                             1]
        ncolor <- length(x = cols.use)
        ngroups <- length(x = unique(x = grouping.data))
        if (ncolor < ngroups) {
                stop(paste("Not enough colors supplied for number of grouping variables. Need", 
                           ngroups, "got", ncolor, "colors"))
        }
        else if (ncolor > ngroups) {
                cols.use <- cols.use[1:ngroups]
        }
        idents.old <- levels(x = object@ident)
        idents.new <- paste(object@ident, grouping.data, sep = "_")
        colorlist <- cols.use
        names(x = colorlist) <- levels(x = grouping.data)
        object@ident <- factor(x = idents.new, 
                               levels = unlist(x = lapply(X = idents.old,
                                                          FUN = function(x) {
                                                                  lvls <- list()
                                                                  for (i in seq_along(along.with = levels(x = grouping.data))) {
                                                                          lvls[[i]] <- paste(x, levels(x = grouping.data)[i], 
                                                                                             sep = "_")
                                                                  }
                                                                  return(unlist(x = lvls))
                                                          })), ordered = TRUE)
        data.to.plot <- data.frame(FetchData(object = object, vars.all = genes.plot))
        data.to.plot$cell <- rownames(x = data.to.plot)
        data.to.plot$id <- object@ident
        data.to.plot <- data.to.plot %>% tidyr::gather(key = genes.plot, 
                                                       value = expression, -c(cell, id))
        data.to.plot <- data.to.plot %>% group_by(id, genes.plot) %>% 
                summarize(avg.exp = ExpMean(x = expression), pct.exp = Seurat:::PercentAbove(x = expression, 
                                                                                             threshold = 0))
        data.to.plot <- data.to.plot %>% ungroup() %>% group_by(genes.plot) %>% 
                mutate(avg.exp = scale(x = avg.exp)) %>% 
                mutate(avg.exp.scale = as.numeric(x = cut(x = MinMax(data = avg.exp,
                                                                     max = col.max,
                                                                     min = col.min),
                                                          breaks = 20)))
        data.to.plot <- data.to.plot %>% tidyr::separate(col = id, into = c("ident1", 
                                                                            "ident2"), sep = "_") %>% rowwise() %>% mutate(palette.use = colorlist[[ident2]],
                                                                                                                           ptcolor = colorRampPalette(colors = c("grey", palette.use))(20)[avg.exp.scale]) %>% 
                tidyr::unite("id", c("ident1", "ident2"), sep = "_")
        data.to.plot$genes.plot <- factor(x = data.to.plot$genes.plot, 
                                          levels = rev(x = sub(pattern = "-", replacement = ".", 
                                                               x = genes.plot)))
        data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
        data.to.plot$id <- factor(x = data.to.plot$id, levels = levels(object@ident))
        palette.use <- unique(x = data.to.plot$palette.use)
        if (!missing(x = gene.groups)) {
                names(x = gene.groups) <- genes.plot
                data.to.plot <- data.to.plot %>% mutate(gene.groups = gene.groups[genes.plot])
        }
        p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot, 
                                                       y = id)) + geom_point(mapping = aes(size = pct.exp, 
                                                                                           color = ptcolor)) + scale_radius(range = c(0, dot.scale)) + 
                scale_color_identity() + theme(axis.title.x = element_blank(), 
                                               axis.title.y = element_blank())
        if (!missing(x = gene.groups)) {
                p <- p + facet_grid(facets = ~gene.groups, scales = "free_x", 
                                    space = "free_x", switch = "y") + theme(panel.spacing = unit(x = 1, 
                                                                                                 units = "lines"), strip.background = element_blank(), 
                                                                            strip.placement = "outside")
        }
        if (x.lab.rot) {
                p <- p + theme(axis.text.x = element_text(angle = 90, 
                                                          vjust = 0.5))
        }
        if (!plot.legend) {
                p <- p + theme(legend.position = "none")
        }
        else if (plot.legend) {
                plot.legend <- cowplot::get_legend(plot = p)
                palettes <- list()
                for (i in seq_along(along.with = colorlist)) {
                        palettes[[names(colorlist[i])]] <- colorRampPalette(colors = c("grey", 
                                                                                       colorlist[[i]]))(20)
                }
                gradient.legends <- mapply(FUN = Seurat:::GetGradientLegend, 
                                           palette = palettes, group = names(x = palettes), 
                                           SIMPLIFY = FALSE, USE.NAMES = FALSE)
                p <- p + theme(legend.position = "none")
                legends <- cowplot::plot_grid(plotlist = gradient.legends, 
                                              plot.legend, ncol = 1, rel_heights = c(1, rep.int(x = 0.5, 
                                                                                                times = length(x = gradient.legends))), scale = rep(0.5, 
                                                                                                                                                    length(gradient.legends)), align = "hv")
                p <- cowplot::plot_grid(p, legends, ncol = 2, rel_widths = c(1, 
                                                                             0.3), scale = c(1, 0.8))
        }
        suppressWarnings(print(p))
        if (do.return) {
                return(p)
        }
}

# add x.lab
VlnPlot.1 <- function (object, features.plot, ident.include = NULL, nCol = NULL, 
          do.sort = FALSE, y.max = NULL, same.y.lims = FALSE, size.x.use = 16, 
          size.y.use = 16, size.title.use = 20, adjust.use = 1, point.size.use = 1, 
          cols.use = NULL, group.by = NULL, y.log = FALSE, x.lab.rot = FALSE, 
          y.lab.rot = FALSE, legend.position = "right", single.legend = TRUE, 
          remove.legend = FALSE, do.return = FALSE, return.plotlist = FALSE, 
          x.lab = "Identity",...) 
{
        if (is.null(x = nCol)) {
                if (length(x = features.plot) > 9) {
                        nCol <- 4
                }
                else {
                        nCol <- min(length(x = features.plot), 3)
                }
        }
        data.use <- data.frame(FetchData(object = object, vars.all = features.plot, 
                                         ...), check.names = F)
        if (is.null(x = ident.include)) {
                cells.to.include <- object@cell.names
        }
        else {
                cells.to.include <- WhichCells(object = object, ident = ident.include)
        }
        data.use <- data.use[cells.to.include, , drop = FALSE]
        if (!is.null(x = group.by)) {
                ident.use <- as.factor(x = FetchData(object = object, 
                                                     vars.all = group.by)[cells.to.include, 1])
        }
        else {
                ident.use <- object@ident[cells.to.include]
        }
        gene.names <- colnames(x = data.use)[colnames(x = data.use) %in% 
                                                     rownames(x = object@data)]
        if (single.legend) {
                remove.legend <- TRUE
        }
        if (same.y.lims && is.null(x = y.max)) {
                y.max <- max(data.use)
        }
        plots <- lapply(X = features.plot, FUN = function(x) {
                return(SingleVlnPlot.1(feature = x, data = data.use[,x, drop = FALSE], cell.ident = ident.use, do.sort = do.sort, 
                                     y.max = y.max, size.x.use = size.x.use, size.y.use = size.y.use, 
                                     size.title.use = size.title.use, adjust.use = adjust.use, 
                                     point.size.use = point.size.use, cols.use = cols.use, 
                                     gene.names = gene.names, y.log = y.log, x.lab.rot = x.lab.rot, 
                                     y.lab.rot = y.lab.rot, legend.position = legend.position, 
                                     remove.legend = remove.legend, x.lab = x.lab))
        })
        if (length(x = features.plot) > 1) {
                plots.combined <- plot_grid(plotlist = plots, ncol = nCol)
                if (single.legend && !remove.legend) {
                        legend <- get_legend(plot = plots[[1]] + theme(legend.position = legend.position))
                        if (legend.position == "bottom") {
                                plots.combined <- plot_grid(plots.combined, 
                                                            legend, ncol = 1, rel_heights = c(1, 0.2))
                        }
                        else if (legend.position == "right") {
                                plots.combined <- plot_grid(plots.combined, 
                                                            legend, rel_widths = c(3, 0.3))
                        }
                        else {
                                warning("Shared legends must be at the bottom or right of the plot")
                        }
                }
        }
        else {
                plots.combined <- plots[[1]]
        }
        if (do.return) {
                if (return.plotlist) {
                        return(plots)
                }
                else {
                        return(plots.combined)
                }
        }
        else {
                if (length(x = plots.combined) > 1) {
                        plots.combined
                }
                else {
                        invisible(x = lapply(X = plots.combined, FUN = print))
                }
        }
}


# define topGOterms funtion ===
#' topGoterm incorporate the whole topgo analysis pipeline
#' @param int.genes genes of interest, extacted based DE analysis
#' @param bg.genes background genes
#' @param organism Select Human or mouse genes
#' @param getBM biomaRt::getBM result. A table with gene name and go id.
#' If no getBM is provided, function will find it automatically.
#' @param ontology Specify ontology database, select within in c("BP","CC","MF"). Default is "BP".
#' @param nodeSize an integer larger or equal to 1. Inherited from topGOdata ojbect
#' @param topNodes an integer larger or equal to 1. Show many rows in final GenTable results.
#' @export results GenTable result.
#' @example 
#' Spermatogonia.Go <-  topGOterms(int.genes = unique(Spermatogonia_markers),
#                                  bg.genes = rownames(SSCs@data),
#'                                 organism =  "Mouse",
#'                                 getBM = getBM)
topGOterms <- function(int.genes, bg.genes,organism =  "Mouse",getBM = NULL,
                       ontology = "BP",
                       nodeSize = 5,topNodes = 20){
        
        if(is.null(getBM)){
                library(biomaRt)
                Database = useMart("ensembl") %>% listDatasets() %>% as.data.frame()
                # mmusculus_gene_ensembl is at 2nd
                dataset = tail(Database[grep(organism,Database$description),"dataset"],1) 
                print(paste("Use biomaRt",dataset,"dataset"))
                # collect gene names from biomart
                mart <- biomaRt::useMart("ensembl",dataset=dataset)
                # Get ensembl gene ids and GO terms
                getBM <- biomaRt::getBM(attributes = c("external_gene_name", "go_id"),
                                        filters= "external_gene_name",
                                        values = bg.genes,
                                        mart = mart)
        }
        getBM <- getBM[getBM$go_id != '',]
        geneID2GO <- by(getBM$go_id, getBM$external_gene_name, function(x) as.character(x))
        all.genes <- sort(unique(as.character(getBM$external_gene_name)))
        print(paste("Total",length(all.genes),"genes.",length(int.genes), "gene of interest"))
        int.genes <- factor(as.integer(all.genes %in% int.genes))
        names(int.genes) = all.genes
        
        tab = as.list(ontology)
        names(tab) = ontology
        for(ont in ontology){
                go.obj <- new("topGOdata", ontology = ont,
                              allGenes = int.genes,
                              annot = annFUN.gene2GO,
                              gene2GO = geneID2GO,
                              nodeSize = nodeSize
                )
                resultFisher <- runTest(go.obj, algorithm = "elim", statistic = "fisher")
                resultKS <- runTest(go.obj, algorithm = "classic", statistic = "ks")
                resultKS.elim <- runTest(go.obj, algorithm = "elim", statistic = "ks")
                
                tab$ont <- GenTable(go.obj, classicFisher = resultFisher, 
                                    classicKS = resultKS, elimKS = resultKS.elim,topNodes = topNodes,
                                    orderBy = "elimKS", ranksOf = "classicFisher")
        }
        topGOResultsSxT <- as.data.frame(do.call(rbind,tab))
        topGOResultsSxT <- topGOResultsSxT[(length(ontology)+1):nrow(topGOResultsSxT),]
        rownames(topGOResultsSxT) = NULL
        showSigOfNodes(go.obj, score(resultFisher), firstSigNodes = 5, useInfo = 'def')
        
        return(results.tab)
}


Types2Markers <- function(df, by = "Cell_Type"){
        "Convert a list of N gene vector into a NX2 dataframe.
        Input
        df: N by 2 dataframe with gene at column1 and Cell_Type at column2.
        by: certerial to group by, chose between Cell_Type, or cluster.
        
        Output 
        List: list of N gene vector
        ---------------------------"
        # preparation
        df <- df[,c("gene","Cell_Type","cluster")]
        df[,"Cell_Type"] <- gsub(" ","_",df[,"Cell_Type"])
        df[,"cluster"] <- gsub(" ","_",df[,"cluster"])
        df <- df[gtools::mixedorder(df[,by]),]
        
        # Create a list contains marker gene vector
        x <- unique(df[,by])
        List <- list()
        List <- lapply(seq_along(x), function(y, df, i) {
                List[[i]]=unique(df[df[,by]==y[i],"gene"])},
                y=x,df=df)
        names(List) <- x
        return(List)
}


TSNEPlot.1 <- function(object, do.label = FALSE, pt.size = 1, label.size = 4, 
                       cells.use = NULL,colors.use = NULL, 
                       text.repel = FALSE, label.repel = TRUE,force= 1,...) 
{
        return(DimPlot.1(object = object, reduction.use = "tsne", 
                         cells.use = cells.use, pt.size = pt.size, do.label = do.label, 
                         label.size = label.size, cols.use = colors.use,
                         text.repel = text.repel, label.repel = label.repel,force= force,...))
}

#' Generate 3D TSNEplot
#' @param object Seurat object after performing RunTSNE(dim.embed = 3)
#' @param ... all other parameters are inherited from Seurat::TSNEPlot()
#' @export p/p3 ggplot object from TSNEplot
TSNEPlot.3D <- function (object, reduction.use = "tsne", dim.1 = 1, dim.2 = 2, dim.3 = 3,
                         cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE, 
                         cols.use = NULL, group.by = "ident", pt.shape = NULL, do.hover = FALSE, 
                         data.hover = "ident", do.identify = FALSE, do.label = FALSE, 
                         label.size = 4, no.legend = FALSE, no.axes = FALSE, dark.theme = FALSE, 
                         plot.order = NULL, plot.title = NULL, ...) 
{
        embeddings.use = GetDimReduction(object = object, reduction.type = reduction.use, 
                                         slot = "cell.embeddings")
        if (length(x = embeddings.use) == 0) {
                stop(paste(reduction.use, "has not been run for this object yet."))
        }
        if (ncol(x = embeddings.use) < 3) {
                stop(paste(reduction.use, "doesn't have the third dimension.
                           Suggest performing RunTSNE(dim.embed = 3)"))
        }
        cells.use <- Seurat:::SetIfNull(x = cells.use, default = colnames(x = object@data))
        dim.code <- GetDimReduction(object = object, reduction.type = reduction.use, 
                                    slot = "key")
        dim.codes <- paste0(dim.code, c(dim.1, dim.2))
        data.plot <- as.data.frame(x = embeddings.use)
        cells.use <- intersect(x = cells.use, y = rownames(x = data.plot))
        data.plot <- data.plot[cells.use, dim.codes]
        ident.use <- as.factor(x = object@ident[cells.use])
        if (group.by != "ident") {
                ident.use <- as.factor(x = FetchData(object = object, 
                                                     vars.all = group.by)[cells.use, 1])
        }
        data.plot$ident <- ident.use
        data.plot$x <- data.plot[, dim.codes[1]]
        data.plot$y <- data.plot[, dim.codes[2]]
        data.plot$z <- data.plot[, dim.codes[3]]
        data.plot$pt.size <- pt.size
        if (!is.null(plot.order)) {
                if (any(!plot.order %in% data.plot$ident)) {
                        stop("invalid ident in plot.order")
                }
                plot.order <- rev(c(plot.order, setdiff(unique(data.plot$ident), 
                                                        plot.order)))
                data.plot$ident <- factor(data.plot$ident, levels = plot.order)
                data.plot <- data.plot[order(data.plot$ident), ]
        }
        rgl::open3d()
        
        car::scatter3d(x = data.plot$x, y = data.plot$y, z = data.plot$z)
        
        geom_point(mapping = aes(colour = factor(x = ident)), 
                   size = pt.size)
        if (!is.null(x = pt.shape)) {
                shape.val <- FetchData(object = object, vars.all = pt.shape)[cells.use, 
                                                                             1]
                if (is.numeric(shape.val)) {
                        shape.val <- cut(x = shape.val, breaks = 5)
                }
                data.plot[, "pt.shape"] <- shape.val
                p <- ggplot(data = data.plot, mapping = aes(x = x, y = y)) + 
                        geom_point(mapping = aes(colour = factor(x = ident), 
                                                 shape = factor(x = pt.shape)), size = pt.size)
        }
        if (!is.null(x = cols.use)) {
                p <- p + scale_colour_manual(values = cols.use)
        }
        p2 <- p + xlab(label = dim.codes[[1]]) + ylab(label = dim.codes[[2]]) + 
                zlab(label = dim.codes[[3]]) +
                scale_size(range = c(pt.size, pt.size))
        p3 <- p2 + Seurat:::SetXAxisGG() + Seurat:::SetYAxisGG()  + 
                Seurat:::SetLegendPointsGG(x = 6) + 
                Seurat:::SetLegendTextGG(x = 12) + no.legend.title + theme_bw() + 
                Seurat:::NoGrid()
        p3 <- p3 + theme(legend.title = element_blank())
        if (!is.null(plot.title)) {
                p3 <- p3 + ggtitle(plot.title) + theme(plot.title = element_text(hjust = 0.5))
        }
        if (do.label) {
                centers <- data.plot %>% dplyr::group_by(ident) %>% 
                        summarize(x = median(x = x), y = median(x = y))
                p3 <- p3 + geom_point(data = centers, mapping = aes(x = x, 
                                                                    y = y), size = 0, alpha = 0) + 
                        geom_text(data = centers,
                                  mapping = aes(label = ident), size = label.size)
        }
        if (dark.theme) {
                p <- p + DarkTheme()
                p3 <- p3 + DarkTheme()
        }
        if (no.legend) {
                p3 <- p3 + theme(legend.position = "none")
        }
        if (no.axes) {
                p3 <- p3 + theme(axis.line = element_blank(), axis.text.x = element_blank(), 
                                 axis.text.y = element_blank(), axis.ticks = element_blank(), 
                                 axis.title.x = element_blank(), axis.title.y = element_blank(), 
                                 panel.background = element_blank(), panel.border = element_blank(), 
                                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                 plot.background = element_blank())
        }
        if (do.identify || do.hover) {
                if (do.bare) {
                        plot.use <- p
                }
                else {
                        plot.use <- p3
                }
                if (do.hover) {
                        if (is.null(x = data.hover)) {
                                features.info <- NULL
                        }
                        else {
                                features.info <- FetchData(object = object, 
                                                           vars.all = data.hover)
                        }
                        return(HoverLocator(plot = plot.use, data.plot = data.plot, 
                                            features.info = features.info, dark.theme = dark.theme))
                }
                else if (do.identify) {
                        return(FeatureLocator(plot = plot.use, data.plot = data.plot, 
                                              dark.theme = dark.theme, ...))
                }
        }
        if (do.return) {
                if (do.bare) {
                        return(p)
                }
                else {
                        return(p3)
                }
        }
}

#' test a matirx's max/median/min of row/column'
testMMM <- function(x,MARGIN = 2) {
        par(mfrow=c(3,1))
        Max <- apply(x,MARGIN,max)
        Median <- apply(x,MARGIN,median)
        Min <- apply(x,MARGIN,min)
        xlim = c(min(x),max(x))
        hist(Max, xlim = xlim)
        hist(Median, xlim = xlim)
        hist(Min, xlim = xlim)
        par(mfrow=c(1,1))
}