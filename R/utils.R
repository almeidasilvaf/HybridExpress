
#' Get data frames of nodes and edges for expression triangle
#'
#' @return A list of data frames with the following elements:
#' \describe{
#'   \item{nodes}{Data frame with node coordinates.}
#'   \item{edges}{Data frame with edge coordinates.}
#'   \item{edges_curves}{Data frame with curved edge coordinates.}
#' }
#'
#' @noRd
get_triangle_graph <- function() {
    
    # Data frame of nodes
    triangle_nodes <- data.frame(
        node = c("P1", "P2", "F1", "Midparent"),
        x = c(0, 4, 2, 2),
        y = c(0, 0, 4, 1.5)
    )
    triangle_nodes$xmin <- triangle_nodes$x - 0.35
    triangle_nodes$xmax <- triangle_nodes$x + 0.35
    
    triangle_nodes$ymin <- triangle_nodes$y - 0.25
    triangle_nodes$ymax <- triangle_nodes$y + 0.25
    
    # Data frame of edges
    triangle_edges_lines <- data.frame(
        x = c(0, 4, 0.35, 2),
        xend = c(1.65, 2.35, 3.65, 2),
        y = c(0.25, 0.25, 0, 1.75),
        yend = c(3.75, 3.75, 0, 3.75),
        label = c("P1-F1", "P2-F1", "P1-P2", "F1-Mid")
    )
    
    triangle_edges_curves <- data.frame(
        x = c(0.35, 2.05),
        xend = c(1.95, 3.65),
        y = c(0.1, 1.25),
        yend = c(1.25, 0.1),
        label = c("P1-Mid", "P2-Mid")
    )
    
    # Store results in a list
    triangle_graph <- list(
        nodes = triangle_nodes,
        edges = triangle_edges_lines,
        edges_curves = triangle_edges_curves
    )
    
    return(triangle_graph)
}


#' Get numbers of differentially expressed genes for expression triangle
#'
#' @param deg_counts Data frame with number of differentially expressed
#' genes per contrast as returned by \code{get_deg_counts}.
#'
#' @return A list of data frames with the following elements:
#' \describe{
#'   \item{total}{Coordinates for the total number of DEGs.}
#'   \item{up}{Coordinates for the number of up-regulated genes.}
#' }
#'
#' @noRd
get_triangle_numbers <- function(deg_counts) {
    
    # Total number of DEGs per contrast
    triangle_total <- data.frame(
        contrast = c("P2_vs_P1", "F1_vs_P1", "F1_vs_P2", "F1_vs_midparent"),
        x = c(2, 0.60, 3.40, 2.25),
        y = c(-0.25, 2.5, 2.5, 2.75)
    )
    triangle_total <- merge(triangle_total, deg_counts, by = "contrast")
    triangle_total$label <- paste0(
        triangle_total$total, "\n(", triangle_total$perc_total, "%)"
    )
    
    # Number of up-regulated genes in each generation
    triangle_up <- data.frame(
        generation = c(
            "P2_vs_P1-P1", "P2_vs_P1-P2", 
            "F1_vs_P1-P1", "F1_vs_P1-F1",
            "F1_vs_P2-P2", "F1_vs_P2-F1",
            "F1_vs_midparent-midparent", "F1_vs_midparent-F1"
        ),
        n_up = c(
            deg_counts$down[deg_counts$contrast == "P2_vs_P1"],
            deg_counts$up[deg_counts$contrast == "P2_vs_P1"],
            deg_counts$down[deg_counts$contrast == "F1_vs_P1"],
            deg_counts$up[deg_counts$contrast == "F1_vs_P1"],
            deg_counts$down[deg_counts$contrast == "F1_vs_P2"],
            deg_counts$up[deg_counts$contrast == "F1_vs_P2"],
            deg_counts$down[deg_counts$contrast == "F1_vs_midparent"],
            deg_counts$up[deg_counts$contrast == "F1_vs_midparent"]
        ),
        perc_up = c(
            deg_counts$perc_down[deg_counts$contrast == "P2_vs_P1"],
            deg_counts$perc_up[deg_counts$contrast == "P2_vs_P1"],
            deg_counts$perc_down[deg_counts$contrast == "F1_vs_P1"],
            deg_counts$perc_up[deg_counts$contrast == "F1_vs_P1"],
            deg_counts$perc_down[deg_counts$contrast == "F1_vs_P2"],
            deg_counts$perc_up[deg_counts$contrast == "F1_vs_P2"],
            deg_counts$perc_down[deg_counts$contrast == "F1_vs_midparent"],
            deg_counts$perc_up[deg_counts$contrast == "F1_vs_midparent"]
        ),
        x = c(
            0.55, 3.45, # P2-P1
            -0.35, 1.40, # F1-P1
            4.35, 2.60, # F1-P2
            1.75, 1.75 # F1-midparent
        ),
        y = c(
            -0.25, -0.25, # P2-P1
            0.55, 4, # F1-P1
            0.55, 4, # F1-P2
            2.05, 3.45 # F1-midparent
        )
    )
    triangle_up$label <- paste0(
        triangle_up$n_up, "\n(", triangle_up$perc_up, "%)"
    )
    
    # Store results in a list
    triangle_numbers <- list(
        total = triangle_total,
        up = triangle_up
    )
    
    return(triangle_numbers)
}



#' Get a summary table of differential expression results per gene in all contrasts
#' 
#' @param deg_list A list of data frames with gene-wise test statistics for
#' differentially expressed genes as returned by \code{get_deg_list()}.
#' 
#' @return A data frame with the following variables:
#' \describe{
#'   \item{Gene}{Character, gene ID.}
#'   \item{P2_vs_P1}{Character, DE results. One of "up", "down", or 
#'                   "unchanged".}
#'   \item{F1_vs_P1}{Character, DE results. One of "up", "down", or 
#'                   "unchanged".}
#'   \item{F1_vs_P2}{Character, DE results. One of "up", "down", or 
#'                   "unchanged".}
#'   \item{F1_vs_midparent}{Character, DE results. One of "up", "down", or 
#'                          "unchanged".}
#' }
#' 
#' @noRd
get_deg_summary <- function(deg_list) {
    
    summary_df <- lapply(seq_along(deg_list), function(x) {
        
        contrast_name <- names(deg_list)[x]
        df <- deg_list[[x]]
        df[, contrast_name] <- ifelse(df$log2FoldChange > 0, "up", "down")
        df$Gene <- rownames(df)
        
        return(df[, c("Gene", contrast_name)])
    })
    
    final_summary <- Reduce(function(x, y) { 
        return(merge(x, y, by = "Gene", all = TRUE))
    }, summary_df)
    
    final_summary[is.na(final_summary)] <- "unchanged"
    
    return(final_summary)
}


#' Plot a scatterplot of log2 fold changes between contrasts
#' 
#' @param partition_table A data frame with genes per expression partition
#' as returned by \code{expression_partitioning()}.
#' @param palette A character vector with the color palette to use.
#' @param group_by Character indicating the name of the variable 
#' in \strong{partition_table} to use to group genes. One of "Category" or
#' "Class". Default: "Category".
#'
#' @return A ggplot object with a scatterplot.
#' 
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual theme_bw
#' theme labs element_blank geom_vline geom_hline element_text
#' @importFrom rlang .data
#' @noRd
#' 
partition_scatterplot <- function(
        partition_table, palette, group_by = "Category"
) {
    
    p_scatter <- ggplot(
        partition_table, aes(x = .data$lFC_F1_vs_P1, y = .data$lFC_F1_vs_P2)
    ) +
        geom_point(
            aes(color = .data[[group_by]]), alpha = 0.4, show.legend = FALSE
        ) +
        scale_color_manual(values = palette) +
        theme_bw() +
        labs(
            subtitle = expression(Log[2] ~ "fold change between generations"),
            x = "F1 vs P1",
            y = "F1 vs P2"
        ) +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            axis.text = element_text(size = 11)
        ) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        geom_hline(yintercept = 0, linetype = "dashed")
    
    return(p_scatter)
}


#' Plot line plots representing each partition's expression profiles
#' 
#' @param partition_table A data frame with genes per expression partition
#' as returned by \code{expression_partitioning()}.
#' @param palette A character vector with the color palette to use.
#' @param add_n Logical indicating whether to include number of genes in each
#' group or not. Default: TRUE.
#' @param group_by Character indicating the name of the variable 
#' in \strong{partition_table} to use to group genes. One of "Category" or
#' "Class". Default: "Category".
#'
#' @return A list of ggplot objects with line plots.
#' 
#' @importFrom ggplot2 ggplot aes geom_point geom_line theme_bw
#' theme labs element_blank element_text ylim facet_wrap
#' @importFrom rlang .data
#' @noRd
#' 
partition_lineplots <- function(
        partition_table, palette, add_n = TRUE, group_by = "Category"
) {
    
    # Data frame of point coordinates for each group
    pline_data <- data.frame(
        Category = factor(rep(1:12, each = 3)),
        Class = factor(
            rep(c(
                "ADD", "ELD_P1", "DOWN", "ELD_P2", "UP", "UP", 
                "DOWN", "UP", "ELD_P2", "DOWN", "ELD_P1", "ADD"
            ), each = 3),
            levels = c("UP", "DOWN", "ADD", "ELD_P1", "ELD_P2")),
        x = factor(rep(c("P1", "F1", "P2"), 12), levels = c("P1", "F1", "P2")),
        y = c(
            1, 2, 3,
            1, 3, 3,
            2, 1, 3,
            3, 3, 1,
            1, 3, 2,
            2, 3, 1,
            3, 1, 3,
            1, 3, 1,
            1, 1, 3,
            3, 1, 2,
            3, 1, 1,
            3, 2, 1
        )
    )
    
    # Split data frame by `group_by` variable
    pline_data <- split(pline_data, pline_data[[group_by]])
    
    ## Number of genes per level of `group_by`
    n <- as.numeric(table(partition_table$Category))
    
    p_line <- lapply(seq_along(pline_data), function(x) {
        
        subtitle <- paste0(group_by, " ", names(pline_data)[x])
        if(add_n) { subtitle <- paste0(subtitle, ", N = ", n[x]) }
        
        p <- ggplot(pline_data[[x]], aes(x = .data$x, y = .data$y, group = 1)) +
            geom_line(color = palette[x], linewidth = 1) +
            geom_point(size = 2, color = "gray20") +
            theme_bw() +
            labs(x = "", y = "", subtitle = subtitle) +
            theme(
                panel.grid.minor = element_blank(),
                axis.text.x = element_text(size = 11, face = "bold"),
                axis.text.y = element_blank()
            ) +
            ylim(c(0, 4))
        
        if(group_by == "Class") { p <- p + facet_wrap("Category") }
        
        return(p)
    })
    names(p_line) <- names(pline_data)
    
    return(p_line)
}


#' Process color palette for plotting functions
#'
#' @param palette Character vector with user-defined color palette, or NULL
#' if the user does not pass a custom color palette. Default: NULL.
#' @param plot Character indicating the name of the plot for which
#' a color palette will be extracted if \strong{palette = NULL}. One 
#' of 'triangle', 'part1', 'part2',   
#'
#' @return A character vector of color names.
#' @noRd
#'
ppal <- function(palette = NULL, plot = "triangle") {
    
    col <- palette
    if(is.null(col)) {
        
        if(plot == "triangle") { # Expression triangle
            col <- c("dodgerblue3", "firebrick", "purple4", "mediumorchid4")
        } else if(plot == "partition") { # Partitions
            col <- c(
                "#5050FFFF", "#CE3D32FF", "#749B58FF", "#F0E685FF", 
                "#466983FF", "#BA6338FF", "#5DB1DDFF", "#802268FF", 
                "#6BD76BFF", "#D595A7FF", "gold2", "#837B8DFF"
            )
        } else if(plot == "pca") { # PCA plot
            col <- c("#3182BDFF", "#E6550DFF", "#31A354FF", "#756BB1FF", 
                     "#636363FF", "#6BAED6FF", "#FD8D3CFF", "#74C476FF", 
                     "#9E9AC8FF", "#969696FF", "#9ECAE1FF", "#FDAE6BFF", 
                     "#A1D99BFF", "#BCBDDCFF", "#BDBDBDFF", "#C6DBEFFF", 
                     "#FDD0A2FF", "#C7E9C0FF", "#DADAEBFF", "#D9D9D9FF"
            )
        } else {
            stop("Unknown plot type.")
        }
    }
    
    return(col)
}


#' Extract row and column metadata from `SummarizedExperiment` objects
#'
#' @param se A `SummarizedExperiment` object.
#' @param rowdata_cols Columns to use from the rowData element of the
#' `SummarizedExperiment` object. It can be either a character vector
#' of column names or a numeric vector of column indices.
#' @param coldata_cols Columns to use from the colData element of the
#' `SummarizedExperiment` object. It can be either a character vector
#' of column names or a numeric vector of column indices.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{rowdata}{A data frame of row metadata containing only the selected
#'                  columns.}
#'   \item{coldata}{A data frame of column metadata containing only the
#'                  selected columns.}
#' }
#'
#' @noRd
#' @importFrom SummarizedExperiment rowData colData
se2metadata <- function(se, rowdata_cols = NULL, coldata_cols = NULL) {
    
    final_rowdata <- NA
    final_coldata <- NA
    
    # Extract row and column metadata
    rowdata <- SummarizedExperiment::rowData(se)
    coldata <- SummarizedExperiment::colData(se)
    
    # Extract user-defined columns from metadata
    ## Row metadata
    if(ncol(rowdata) > 0) {
        r_cols <- rowdata_cols
        if(is.null(rowdata_cols)) { r_cols <- seq_along(rowdata) }
        
        final_rowdata <- as.data.frame(rowdata[, r_cols, drop = FALSE])
    }
    
    if(ncol(coldata) > 0) {
        c_cols <- coldata_cols
        if(is.null(coldata_cols)) { c_cols <- seq_along(coldata) }
        final_coldata <- as.data.frame(coldata[, c_cols, drop = FALSE])
    }
    
    # Return resuls as a list
    metadata_list <- list(
        rowdata = final_rowdata,
        coldata = final_coldata
    )
    
    return(metadata_list)
}


#' Map levels of metadata variables to colors for plotting
#'
#' @param metadata A data frame with column or row metadata. If column
#' metadata is passed, row names must contain sample names. If row metadata is
#' passed, row names must contain gene names.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{metadata}{A metadata data frame as in \strong{metadata}, but
#'                   with rows sorted by levels of every column.}
#'   \item{colors}{A list of named character vectors containing the mapping
#'                 between levels of metadata variables and colors.}
#' }
#'
#' @noRd
#' @importFrom stats setNames
metadata2colors <- function(metadata) {
    
    coldata <- NA
    colors <- NA
    
    if(is.data.frame(metadata)) {
        # Get variable names in metadata
        col_names <- names(metadata)
        if(length(col_names) > 3) {
            stop("Maximum number of columns for row and sample metadata is 3.")
        }
        
        # Sort elements in all columns
        coldata <- metadata[do.call(order, metadata), , drop = FALSE]
        
        # Create a list of named vectors with variable levels and colors
        colors <- lapply(seq_along(col_names), function(x) {
            levels <- unique(coldata[, x])
            cols <- setNames(custom_palette(x)[seq_along(levels)], levels)
            return(cols)
        })
        names(colors) <- col_names
    }
    
    # Return results as a list
    results <- list(
        metadata = coldata,
        colors = colors
    )
    
    return(results)
}



#' Generate custom color palette
#'
#' @param pal Numeric specifying palette number, from 1 to 3.
#'
#' @return Character vector of custom color palette with 20 colors
#' @noRd
custom_palette <- function(pal = 1) {
    pal1 <- c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF",
              "#9467BDFF", "#8C564BFF", "#E377C2FF", "#7F7F7FFF",
              "#BCBD22FF", "#17BECFFF", "#AEC7E8FF", "#FFBB78FF",
              "#98DF8AFF", "#FF9896FF", "#C5B0D5FF", "#C49C94FF",
              "#F7B6D2FF", "#C7C7C7FF", "#DBDB8DFF", "#9EDAE5FF")
    
    pal2 <- c("#3182BDFF", "#E6550DFF", "#31A354FF", "#756BB1FF",
              "#636363FF", "#6BAED6FF", "#FD8D3CFF", "#74C476FF",
              "#9E9AC8FF", "#969696FF", "#9ECAE1FF", "#FDAE6BFF",
              "#A1D99BFF", "#BCBDDCFF", "#BDBDBDFF", "#C6DBEFFF",
              "#FDD0A2FF", "#C7E9C0FF", "#DADAEBFF", "#D9D9D9FF")
    
    pal3 <- c("#393B79FF", "#637939FF", "#8C6D31FF", "#843C39FF",
              "#7B4173FF", "#5254A3FF", "#8CA252FF", "#BD9E39FF",
              "#AD494AFF", "#A55194FF", "#6B6ECFFF", "#B5CF6BFF",
              "#E7BA52FF", "#D6616BFF", "#CE6DBDFF", "#9C9EDEFF",
              "#CEDB9CFF", "#E7CB94FF", "#E7969CFF", "#DE9ED6FF")
    
    l <- list(pal1, pal2, pal3)
    l_final <- l[[pal]]
    return(l_final)
}


