# TRUE Global Environment

CreatGeneSetsFromSingleR <- function(object = object, cell.type = "B_cell",main.type = FALSE,
                            Human =TRUE, Mouse = FALSE){
        "Find Marker genes list for one cell type"
        #=== Create Score matrix======
        if(main.type) {
                de_gene_mode = "de.genes.main"; types_mode = "main_types"
        } else {
                de_gene_mode = "de.genes"; types_mode = "types"
        }
        genes = unique(unlist(object[[de_gene_mode]][[cell.type]]))
        types = unique(object[[types_mode]])
        S_matrix = matrix(data = genes, ncol = 1,
                          dimnames = list(NULL, "gene"))
        #=== Calculate ranking score=====
        for (object_type in types){
                Gen <- object[[de_gene_mode]][[cell.type]][[object_type]]
                Gen_mat = matrix(data = c(Gen,length(Gen):1), ncol = 2,
                                 dimnames = list(NULL, c("gene",object_type)))
                S_matrix =  merge(S_matrix,Gen_mat, by="gene",all = T)
        }
        
        S_matrix = S_matrix[!is.na(S_matrix$gene),]
        # Store the gene name
        if(Human) Genes = toupper(S_matrix$gene)
        if(Mouse) Genes = Hmisc::capitalize(tolower(S_matrix$gene))
        S_matrix$gene = 0
        S_matrix = as.matrix(S_matrix)
        S_matrix[is.na(S_matrix)] <- 0
        S_matrix = apply(S_matrix,2, as.numeric)
        S_matrix[,"gene"] = rowSums(S_matrix)
        rownames(S_matrix) = Genes
        markers <-rownames(S_matrix)[order(S_matrix[,"gene"],decreasing = T)]
        
        return(markers)
}


CreatGenePanelFromSingleR <- function(object = object, main.type = TRUE, Human =TRUE, Mouse = FALSE){
        "Find Marker genes list for All cell types"
        if(main.type) {
                types = unique(object$main_types)
        } else { types = unique(object$types)
        }
        marker_list = lapply(types, function(x) {
                CreatGeneSetsFromSingleR(object = object, cell.type = x, main.type = main.type,
                                Human =Human, Mouse = Mouse)
        })
        names(marker_list) <- types
        marker_df <- list2df(marker_list)
        # arrange columns
        marker_df= marker_df[,order(colnames(marker_df))]
        return(marker_df)
}


CreateSinglerReference <- function(name, expr, types, main_types,
de.num = 200,de.main.num = 300){
    ref = list(name=name, data = expr, types=types, main_types=main_types)
    
    # if using the de method, we can predefine the variable genes
    ref$de.genes = CreateVariableGeneSet(expr,types,de.num)
    ref$de.genes.main = CreateVariableGeneSet(expr,main_types,de.main.num)
    
    # if using the sd method, we need to define an sd threshold
    sd = rowSds(expr)
    sd.thres = sort(sd, decreasing = T)[4000] # or any other threshold
    ref$sd.thres = sd.thres
    
    return(ref)
    
}


FineTune_Human <- function(x, main.type = FALSE){
        # for both main types and sub-types
        x = gsub(" ","_",x)
        x = gsub("B-cells","B_cells",x)
        x = gsub("cells","cell",x)
        x = gsub("cell","cells",x)
        x = gsub("cytes","cyte",x)
        x = gsub("cyte","cytes",x)
        x = gsub("blasts","blast",x)
        x = gsub("blast","blasts",x)
        x = gsub("Macrophages","Macrophage",x)
        x = gsub("Macrophage","Macrophages",x)
        x = gsub("Neutrophils","Neutrophil",x)
        x = gsub("Neutrophil","Neutrophils",x)
        x = gsub("Smooth_muscle_cells","Smooth_muscle",x)
        x = gsub("Smooth_muscle","Smooth_muscle_cells",x)
        if(main.type){
                x = gsub("CD4\\+_T-cells","T_cells",x)
                x = gsub("CD8\\+_T-cells","T_cells",x)
                x = gsub("HSC_-G-CSF","HSC",x)
                x = gsub("HSC_CD34\\+","HSC",x)
                x = gsub("Memory_B_cells","B_cells",x)
                x = gsub("naive B-cells","B_cells",x)
                x = gsub("Pre-B_cells_CD34-","B_cells",x)
                x = gsub("Pro-B_cells_CD34\\+","B_cells",x)
                x = gsub("Pro-Myelocytes","Myelocytes",x) 
                x = gsub("Tregs","T_cells",x)
                
        } else {
                x = gsub("CD4\\+_T-cells","T_cells:CD4\\+",x)
                x = gsub("CD8\\+_T-cells","T_cells:CD8\\+",x)
                x = gsub("CD4\\+_Tcm","T_cells:CD4\\+_central_memory",x)
                x = gsub("CD4\\+_Tem","T_cells:CD4\\+_effector_memory",x)
                x = gsub("CD8\\+_Tcm","T_cells:CD8\\+_Central_memory",x)
                x = gsub("CD8\\+_Tem","T_cells:CD8\\+_effector_memory",x)
                x = gsub("Class-switched_memory_B_cells","B_cells:Class-switched_memory",x)
                x = gsub("Macrophages_M1","Macrophages:M1",x)
                x = gsub("Macrophages_M2","Macrophages:M2",x)
                x = gsub("Memory_B_cells","B_cells:Memory",x)
                x = gsub("naive B-cells","B_cells:Naive",x)
                x = gsub("Tregs","T_cells:Tregs",x)
        }
        return(x)
}


SearchMarker <- function(df, marker){
        result = apply(df, 2, function(x) which(grepl(marker, x)))
        result_df = data.frame(sort(unlist(result)))
        colnames(result_df) = marker
        return(result_df)
}


SearchAllMarkers <- function(df, markers){
        results <- lapply(markers,function(x) SearchMarker(df,x))
        names(results) <- markers
        temp_df = results[[1]]
        for(i in 2:length(results)) {
                temp_df = merge(temp_df,results[[i]], by="row.names",all=TRUE)
                rownames(temp_df) = temp_df$Row.names
                temp_df = temp_df[,-1]
        }
        temp_df = removeNA(temp_df)
        return(temp_df)
}


SingleR.Subset.1 <- function (singler, subsetdata)
{
    s = singler
    if (!is.null(s$seurat)) {
        s$seurat = SubsetData(s$seurat, colnames(s$seurat@data)[subsetdata])
        subsetdata = unlist(lapply(s$seurat@cell.names, FUN = function(x) which(singler$singler[[1]]$SingleR.single$cell.names ==
        x)))
    }
    for (i in 1:length(s$singler)) {
        s$singler[[i]]$SingleR.single$cell.names = s$singler[[i]]$SingleR.single$cell.names[subsetdata]
        s$singler[[i]]$SingleR.clusters$cell.names = s$singler[[i]]$SingleR.clusters$cell.names[subsetdata]
        s$singler[[i]]$SingleR.single$scores = s$singler[[i]]$SingleR.single$scores[subsetdata,
        ]
        s$singler[[i]]$SingleR.single$labels = as.matrix(s$singler[[i]]$SingleR.single$labels[subsetdata,
        ])
        
        if(length(s$singler[[i]]$SingleR.single$labels1)>1) {
            s$singler[[i]]$SingleR.single$labels1 = as.matrix(s$singler[[i]]$SingleR.single$labels1[subsetdata,
            ])}
        s$singler[[i]]$SingleR.single$clusters$cl = s$singler[[i]]$SingleR.single$clusters$cl[subsetdata]
        if (!is.null(s$singler[[i]]$SingleR.single.main)) {
            s$singler[[i]]$SingleR.single.main$cell.names = s$singler[[i]]$SingleR.single.main$cell.names[subsetdata]
            s$singler[[i]]$SingleR.clusters.main$cell.names = s$singler[[i]]$SingleR.clusters.main$cell.names[subsetdata]
            s$singler[[i]]$SingleR.single.main$scores = s$singler[[i]]$SingleR.single.main$scores[subsetdata,
            ]
            s$singler[[i]]$SingleR.single.main$labels = as.matrix(s$singler[[i]]$SingleR.single.main$labels[subsetdata,
            ])
            if(length(s$singler[[i]]$SingleR.single$labels1)>1) {
                s$singler[[i]]$SingleR.single.main$labels1 = as.matrix(s$singler[[i]]$SingleR.single.main$labels1[subsetdata,
                ])}
            s$singler[[i]]$SingleR.single.main$clusters$cl = s$singler[[i]]$SingleR.single.main$clusters$cl[subsetdata]
        }
    }
    if (!is.null(s[["signatures"]])) {
        s$signatures = s$signatures[subsetdata, ]
    }
    if (!is.null(s[["other"]])) {
        s$other = s$other[subsetdata]
    }
    if (!is.null(s$meta.data)) {
        s$meta.data$orig.ident = factor(as.character(s$meta.data$orig.ident[subsetdata]))
        s$meta.data$xy = s$meta.data$xy[subsetdata, ]
        s$meta.data$clusters = factor(as.character(s$meta.data$clusters[subsetdata]))
    }
    s
}


SingleR.PlotTsne.1 <- function (SingleR, xy, labels = SingleR$labels, score.thres = 0,
                clusters = NULL, do.letters = TRUE, dot.size = 1, do.labels = FALSE,
                do.legend = TRUE, label.size = 3, title = "", colors = singler.colors,
                font.size = NULL, alpha = 0.5, label.repel = TRUE, text.repel = FALSE, force = 1)
{
    if (do.labels == TRUE)
    do.letters = FALSE
    df = data.frame(row.names = SingleR$cell.names)
    df$x = xy[, 1]
    df$y = xy[, 2]
    if (SingleR$method == "cluster") {
        df$ident = clusters.map.values(clusters, labels)
    }
    else {
        df$ident = labels
    }
    if (score.thres > 0) {
        max.score = apply(SingleR$scores, 1, max)
        df$ident[max.score < score.thres] = "X"
    }
    df$ident = factor(df$ident)
    SYMBOLS = c(LETTERS, tolower(letters), c(0:9))
    df$initIdent = SYMBOLS[as.numeric(df$ident)]
    num.levels = length(levels(df$ident))
    p = ggplot(df, aes(x = x, y = y,color = ident))
    p = p + geom_point(aes(color = ident), size = dot.size,
    alpha = alpha, stroke = 0)
    if (do.letters == TRUE) {
        symbols = SYMBOLS[1:num.levels]
        names(symbols) = lev = levels(df$ident)
        p = p + geom_point(aes(shape = ident), size = 2 * dot.size/5,
        color = "black")
        p = p + scale_shape_manual(values = symbols)
    }
    if (do.labels == TRUE) {
        centers <- df %>% dplyr::group_by(ident) %>% dplyr::summarize(x = median(x),
        y = median(y))
        p = p + geom_point(data = centers, aes(x = x, y = y),
        size = 0, alpha = 0)
        if (label.repel == TRUE) {
            p = p + ggrepel::geom_label_repel(data = centers,
            aes(label = ident),
            size = label.size,
            force = force)
        }
        else if (text.repel == TRUE){
            p = p + ggrepel::geom_text_repel(data = centers,
            aes(label = ident),
            force = force,
            size = label.size,
            color = "black")
        }
        if (label.repel == FALSE & text.repel == FALSE) {
            p = p + geom_text(data = centers,
            aes(label = ident), size = label.size, color = "black")
        }
        p = p + guides(colour = guide_legend(override.aes = list(alpha = 1)))
        x.range = layer_scales(p)$x$range$range
        add_to_x = sum(abs(x.range)) * 0.03
        p = p + xlim(x.range[1] - add_to_x, x.range[2] + add_to_x)
    }
    else {
        if (is.null(font.size)) {
            font.size = 250 * (1/num.levels)
            font.size = max(font.size, 5)
            font.size = min(font.size, 10)
        }
        if (num.levels > 35 & num.levels < 60) {
            p = p + theme(legend.position = "bottom", legend.direction = "vertical",
            legend.text = element_text(size = 6), legend.title = element_blank()) +
            guides(col = guide_legend(ncol = 5, override.aes = list(size = 2,
            alpha = 1)))
        }
        else if (num.levels > 60) {
            p = p + theme(legend.position = "bottom", legend.direction = "vertical",
            legend.text = element_text(size = 6), legend.title = element_blank()) +
            guides(col = guide_legend(ncol = 9, override.aes = list(size = 2,
            alpha = 1)))
        }
        else {
            p = p + theme(legend.text = element_text(size = font.size),
            legend.title = element_blank()) + guides(color = guide_legend(ncol = 1,
            override.aes = list(size = 3, alpha = 1)))
        }
    }
    lev = levels(df$ident)
    cols = colors[1:length(lev)]
    names(cols) = lev
    cols[names(cols) == "X"] = "black"
    p = p + scale_color_manual(values = cols)
    p = p + xlab("tSNE 1") + ylab("tSNE 2") + ggtitle(title)
    if (do.legend == FALSE) {
        p = p + theme(legend.position = "none")
    }
    p = p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))
    return(p)
}


removeNA <- function(df){
        AllNA = apply(df,2,function(x) length(which(!is.na(x))))
        df = df[,AllNA !=0]
        return(df)
}


SplitSingler <- function(singler = singler, split.by = "conditions"){
    "
    split singler object by certein criteria
    
    Inputs
    -------------------
    singler: singler object with seurat
    split.by: the criteria to split
    
    Outputs
    --------------------
    Singler.subsets: list of subseted singler object by certein conditions,
    plus levels of the conditions
    "
    if(is.null(singler$seurat))
    stop("A seurat object must be provided first, add singler$seurat = your_seurat_object")
    cell.subsets <- SplitCells(object = singler$seurat, split.by = split.by)
    
    Singler.subsets <- list()
    for(i in 1:(length(cell.subsets)-1)){
        cell.index <- which(singler$seurat@cell.names %in% cell.subsets[[i]])
        Singler.subsets[[i]] <- SingleR.Subset.1(singler=singler,
        subsetdata =cell.index)
    }
    Singler.subsets[[i+1]] <- cell.subsets[[i+1]] # record conditions in the last return
    return(Singler.subsets)
}

SplitSingleR.PlotTsne <- function(singler = singler, split.by = "conditions",
select.plots = NULL, return.plots = FALSE,
do.label= TRUE,do.letters = FALSE,main =FALSE,
show.2nd =FALSE,label.size = 4, dot.size = 3,
legend.size = NULL,... ){
    "
    split singler by certein criteria, and generate TSNE plot
    
    Inputs
    -------------------
    singler: singler object with seurat
    split.by: the criteria to split
    select.plots：change order to output, such as c(2,1)
    return.data: TRUE/FASLE, return splited ojbect or not.
    show.subtype: TRUE/FASLE, show sub cell type or not.
    
    Outputs
    --------------------
    return.plots: if return.data = TRUE
    "
    object.subsets <- SplitSingler(singler = singler, split.by = "conditions")
    levels <- object.subsets[[length(object.subsets)]]
    
    out <- list()
    p <- list()
    if(is.null(select.plots)) select.plots <- 1:length(levels)
    sp = select.plots
    if(main) main_or_sub = "SingleR.single.main" else main_or_sub = "SingleR.single"
    if(show.2nd) st =2 else st =1
    for(i in 1:length(select.plots)){
        out[[i]] <- SingleR.PlotTsne.1(SingleR= object.subsets[[sp[i]]]$singler[[st]][[main_or_sub]],
        xy = object.subsets[[sp[i]]]$meta.data$xy,
        do.label=do.label,do.letters = do.letters,
        labels = object.subsets[[sp[i]]]$singler[[st]][[main_or_sub]]$labels,
        label.size = label.size, dot.size = dot.size,...)
        out[[i]] = out[[i]] +
        ggtitle(levels[sp[i]])+
        theme(text = element_text(size=20),
        plot.title = element_text(hjust = 0.5))
        if(!is.null(legend.size))
        out[[i]] = out[[i]] + theme(legend.text = element_text(size=legend.size))
    }
    out <- out[lapply(out,length)>0] # remove NULL element
    if(return.plots) return(out) else print(do.call(plot_grid, out))
}
