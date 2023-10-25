
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
#' in \strong{partition_table} to use to group genes. One of "Group" or
#' "Class". Default: "Group".
#'
#' @return A ggplot object with a scatterplot.
#' 
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual theme_bw
#' theme labs element_blank geom_vline geom_hline element_text
#' @importFrom rlang .data
#' @noRd
#' 
partition_scatterplot <- function(partition_table, palette, group_by = "Group") {
    
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
#' in \strong{partition_table} to use to group genes. One of "Group" or
#' "Class". Default: "Group".
#'
#' @return A list of ggplot objects with line plots.
#' 
#' @importFrom ggplot2 ggplot aes geom_point geom_line theme_bw
#' theme labs element_blank element_text ylim facet_wrap
#' @importFrom rlang .data
#' @noRd
#' 
partition_lineplots <- function(
        partition_table, palette, add_n = TRUE, group_by = "Group"
) {
    
    # Data frame of point coordinates for each group
    pline_data <- data.frame(
        Group = factor(rep(1:12, each = 3)),
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
    n <- as.numeric(table(partition_table$Group))
    
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
        
        if(group_by == "Class") { p <- p + facet_wrap("Group") }
        
        return(p)
    })
    names(p_line) <- names(pline_data)
    
    return(p_line)
}




