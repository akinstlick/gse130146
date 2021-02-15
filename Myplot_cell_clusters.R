Myplot_cell_clusters <- function (cds, x = 1, y = 2, color_by = "Cluster", markers = NULL, 
    show_cell_names = FALSE, cell_size = 1.5, cell_name_size = 2, addValue=0, #cell_shape=16,
    ...) 
## function was based on monocle:plot_cell_clusters
## addValue is a new parameter to plot ggplot
## return p to allow more visualization
# revised from monocle::plot_cell_clusters to allow ncol(reducedDimA(cds)) >  ncol(cds)
{
    require(viridis)
    if (is.null(cds@reducedDimA) | length(pData(cds)$Cluster) == 
        0) {
        stop("Error: Clustering is not performed yet. Please call clusterCells() before calling this function.")
    }
    gene_short_name <- NULL
    sample_name <- NULL
    data_dim_1 <- NULL
    data_dim_2 <- NULL
    lib_info <- pData(cds)
    tSNE_dim_coords <- reducedDimA(cds)
    
    tSNE_dim_coords <- tSNE_dim_coords[,1:ncol(cds)]      ## by xy
    
    data_df <- data.frame(t(tSNE_dim_coords[c(x, y), ]))
    colnames(data_df) <- c("data_dim_1", "data_dim_2")
    data_df$sample_name <- colnames(cds)
    data_df <- merge(data_df, lib_info, by.x = "sample_name", 
        by.y = "row.names")
    markers_exprs <- NULL
    if (is.null(markers) == FALSE) {
        markers_fData <- subset(fData(cds), gene_short_name %in% 
            markers)
        if (nrow(markers_fData) >= 1) {
            cds_subset <- cds[row.names(markers_fData), ]
            if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", 
                "negbinomial.size")) {
                integer_expression <- TRUE
            }
            else {
                integer_expression <- FALSE
            }
            if (integer_expression) {
                cds_exprs <- exprs(cds_subset)
                if (is.null(sizeFactors(cds_subset))) {
                  stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
                }
                cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
                cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
            }
            else {
                cds_exprs <- reshape2::melt(as.matrix(exprs(cds_subset)))
            }
            markers_exprs <- cds_exprs
            colnames(markers_exprs)[1:2] <- c("feature_id", 
                "cell_id")
            markers_exprs <- merge(markers_exprs, markers_fData, 
                by.x = "feature_id", by.y = "row.names")
            markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
            markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
        }
    }
    if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 
        0) {
        data_df <- merge(data_df, markers_exprs, by.x = "sample_name", 
            by.y = "cell_id")
        g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2)) + 
            facet_wrap(~feature_label)
    }
    else {
        g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2))
    }
    if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 
        0) {
        # g <- g + geom_point(aes(color = log10(value + addValue),shape=cell_shape),
         g <- g + geom_point(aes(color = log10(value + addValue)),
            size = I(cell_size), na.rm = TRUE) 
         + scale_color_viridis(name = paste0("log10(value + ",addValue,")"), 
            ...)
    }
    else {
         #g <- g + geom_point(aes_string(color = color_by,shape=cell_shape), size = I(cell_size), 
         g <- g + geom_point(aes_string(color = color_by), size = I(cell_size), 
                                                 na.rm = TRUE )   
    }
    g <- g + monocle_theme_opts() + xlab(paste("Component", 
        x)) + ylab(paste("Component", y)) + theme(legend.position = "top", 
        legend.key.height = grid::unit(0.35, "in")) + theme(legend.key = element_blank()) + 
        theme(panel.background = element_rect(fill = "white")) + 
        theme(text = element_text(size = 15))
    g 
    return(g)
}

# https://github.com/cole-trapnell-lab/monocle-release/blob/master/R/plotting.R
monocle_theme_opts <- function()
{

    theme(strip.background = element_rect(colour = 'white', fill = 'white')) +

    theme(panel.border = element_blank()) +

    theme(axis.line.x = element_line(size=0.25, color="black")) +

    theme(axis.line.y = element_line(size=0.25, color="black")) +

    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +

    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 

    theme(panel.background = element_rect(fill='white')) +

    theme(legend.key=element_blank())

}


