
#' Plot a triangle of comparisons of DEG sets among generations
#'
#' @param deg_counts Data frame with number of differentially expressed
#' genes per contrast as returned by \code{get_deg_counts}.
#' @param palette Character vector of length 4 indicating the colors of the 
#' boxes for P1, P2, F1, and midparent, respectively. If NULL, 
#' a default color palette will be used.
#' @param box_labels Character vector of length 4 indicating the labels of the
#' boxes for P1, P2, F1, and midparent, respectively. Default: NULL, which
#' will lead to labels "P1", "P2", 
#' "F1", and "Midparent", respectively.
#'
#' @return A ggplot object with an expression triangle.
#' @details
#' The expression triangle plot shows the number of differentially expressed
#' genes (DEGs) for each contrast. Numbers in the center of the lines (in bold)
#' indicate the total number of DEGs, while numbers near boxes indicate
#' the number of up-regulated genes in each generation of the triangle.
#' 
#' @importFrom ggplot2 ggplot aes geom_segment geom_curve geom_rect 
#' geom_text theme_void ylim
#' @importFrom rlang .data
#' @export
#' @rdname plot_expression_triangle
#' @examples
#' data(deg_counts)
#' plot_expression_triangle(deg_counts)
plot_expression_triangle <- function(
        deg_counts, palette = NULL, box_labels = NULL
) {
    
    pal <- ppal(palette, "triangle")
    
    # Get plot data
    tgraph <- get_triangle_graph()
    numbers <- get_triangle_numbers(deg_counts)
    
    if(!is.null(box_labels)) {
        tgraph$nodes$node <- gsub("P1", box_labels[1], tgraph$nodes$node)
        tgraph$nodes$node <- gsub("P2", box_labels[2], tgraph$nodes$node)
        tgraph$nodes$node <- gsub("F1", box_labels[3], tgraph$nodes$node)
        tgraph$nodes$node <- gsub("Midparent", box_labels[4], tgraph$nodes$node)
    }
    
    # Plot triangle
    p <- ggplot() +
        geom_segment(
            data = tgraph$edges, 
            aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend)
        ) +
        geom_curve(
            data = tgraph$edges_curves,
            aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
            curvature = 0.2, linetype = "dashed"
        ) +
        geom_rect(
            data = tgraph$nodes, 
            aes(
                xmin = .data$xmin, ymin = .data$ymin, 
                xmax = .data$xmax, ymax = .data$ymax
            ),
            color = "grey20", fill = pal
        ) +
        geom_text(
            data = tgraph$nodes, 
            aes(x = .data$x, y = .data$y, label = .data$node), color = "gray90"
        ) +
        geom_text(
            data = numbers$total, 
            aes(x = .data$x, y = .data$y, label = .data$label), 
            fontface = "bold", size = 4
        ) +
        geom_text(
            data = numbers$up,
            aes(x = .data$x, y = .data$y, label = .data$label),
            size = 4
        ) +
        ylim(-0.5, 4.5) +
        theme_void()
    
    return(p)
}


#' Plot expression partitions
#'
#' @param partition_table A data frame with genes per expression partition
#' as returned by \code{expression_partitioning()}.
#' @param group_by Character indicating the name of the variable 
#' in \strong{partition_table} to use to group genes. One of "Group" or
#' "Class". Default: "Group".
#' @param palette Character vector with color names to be used for each level
#' of the variable specified in \strong{group_by}. 
#' If \strong{group_by = "Group"}, this must be a vector of length 12.
#' If \strong{group_by = "Class"}, this must be a vector of length 5.
#' If NULL, a default color palette will be used.
#'
#' @return A ggplot object with a plot showing genes in each expression
#' partition.
#' 
#' @importFrom patchwork wrap_plots
#' @importFrom ggplot2 unit
#' @export
#' @rdname plot_expression_partitions
#' @examples
#' data(deg_list)
#' partition_table <- expression_partitioning(deg_list)
#' plot_expression_partitions(partition_table)
plot_expression_partitions <- function(
        partition_table, group_by = "Group", palette = NULL
) {
    
    pdata <- partition_table
    pdata <- pdata[!is.na(pdata$lFC_F1_vs_P1) & !is.na(pdata$lFC_F1_vs_P2), ]
    
    # Define color palette
    pal <- ppal(palette, "partition")[seq_along(levels(pdata[[group_by]]))]
    names(pal) <- levels(pdata[[group_by]])
    
    # 1) Create scatterplot
    p_scatter <- partition_scatterplot(pdata, pal, group_by = group_by)
    
    # 2) Create line plots to represent partitions
    p_line <- partition_lineplots(pdata, pal, group_by = group_by)
    
    # Combine figures
    if(group_by == "Class") {
        p_final <- wrap_plots(
            wrap_plots(p_line, ncol = 1) & 
                theme(plot.margin = unit(rep(1, 4), "pt")),
            p_scatter,
            ncol = 2,
            widths = c(1, 1.5)
        )
    } else {
        p_final <- wrap_plots(
            wrap_plots(p_line[1:4], nrow = 1), # row 1
            wrap_plots(
                wrap_plots(p_line[c(5, 7)], nrow = 2),
                p_scatter,
                wrap_plots(p_line[c(6, 8)], nrow = 2),
                ncol = 3, widths = c(1, 2, 1)
            ), # row 2
            wrap_plots(p_line[9:12], nrow = 1), # row 3
            ncol = 1,
            heights = c(1, 3, 1)
        )
    }

    return(p_final)
}



#' Plot a barplot of gene frequencies per expression partition
#'
#' @param partition_table A data frame with genes per expression partition
#' as returned by \code{expression_partitioning()}.
#' @param group_by Character indicating the name of the variable 
#' in \strong{partition_table} to use to group genes. One of "Group" or
#' "Class". Default: "Group".
#' @param palette Character vector with color names to be used for each level
#' of the variable specified in \strong{group_by}. 
#' If \strong{group_by = "Group"}, this must be a vector of length 12.
#' If \strong{group_by = "Class"}, this must be a vector of length 5.
#' If NULL, a default color palette will be used.
#' 
#' @return A ggplot object with a barplot showing gene frequencies per
#' partition next to explanatory line plots depicting each partition.
#'
#' @importFrom ggplot2 geom_bar coord_flip scale_x_discrete geom_text
#' scale_y_continuous
#' @export
#' @rdname plot_partition_frequencies
#' @examples
#' data(deg_list)
#' partition_table <- expression_partitioning(deg_list)
#' plot_partition_frequencies(partition_table)
plot_partition_frequencies <- function(
        partition_table, group_by = "Group", palette = NULL
) {
    
    # Define color palette
    lev <- levels(partition_table[[group_by]])
    pal <- ppal(palette, "partition")[seq_along(lev)]
    names(pal) <- lev

    # Get barplot data
    freqs <- table(partition_table[[group_by]])
    freq_df <- data.frame(
        Group = factor(names(freqs), levels = names(freqs)), 
        N = as.numeric(freqs)
    )
    freq_df$Perc <- paste0(round((freq_df$N / sum(freq_df$N)) * 100, 2), "%")
    ymax <- round(max(freq_df$N) + mean(freq_df$N), -2)
    
    # Create barplot
    p_bar <- ggplot(freq_df, aes(x = .data$Group, y = .data$N)) +
        geom_bar(fill = pal, color = "gray20", stat = "identity") +
        geom_text(aes(label = .data$Perc), hjust = -0.2) +
        theme_bw() +
        scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +
        theme(plot.subtitle = element_text(size = 13)) +
        labs(y = "Count", subtitle = "Frequency of genes per partition") +
        scale_x_discrete(limits = rev(levels(freq_df$Group))) +
        coord_flip() 
    
    
    # Get line plots
    p_line <- partition_lineplots(
        partition_table, pal, add_n = FALSE, group_by = group_by
    )
    
    # Combine plots
    ncols <- ifelse(group_by == "Group", 2, 1)
    p_final <- wrap_plots(
        wrap_plots(p_line, ncol = ncols) & 
            theme(plot.margin = unit(rep(1, 4), "pt")),
        p_bar,
        nrow = 1,
        widths = c(1, 2)
    )
    
    return(p_final)
}



#' Perform a principal component analysis (PCA) and plot PCs
#'
#' @param se A `SummarizedExperiment` object with a count matrix and sample
#' metadata.
#' @param PCs Numeric vector indicating which principal components to show
#' in the x-axis and y-axis, respectively. Default: \code{c(1,2)}.
#' @param color_by Character with the name of the column in \code{colData(se)}
#' to use to group samples by color. Default: NULL.
#' @param shape_by Character with the name of the column in \code{colData(se)}
#' to use to group samples by shape. Default: NULL.
#' @param add_mean Logical indicating whether to add a diamond symbol
#' with the mean value for each level of the variable indicated 
#' in \strong{color_by}. Default: FALSE
#' @param palette Character vector with colors to use for each level of the
#' variable indicated in \strong{color_by}. If NULL, a default color palette
#' will be used.
#'
#' @return A ggplot object with a PCA plot showing 2 principal components
#' in each axis along with their % of variance explained.
#'
#' @importFrom SummarizedExperiment assay colData
#' @importFrom DESeq2 vst
#' @importFrom stats prcomp
#' @importFrom ggplot2 geom_point theme scale_color_manual element_blank
#' @export
#' @rdname pca_plot
#' @examples
#' data(se_chlamy)
#' se <- add_midparent_expression(se_chlamy)
#' se$Ploidy[is.na(se$Ploidy)] <- "midparent"
#' se$Generation[is.na(se$Generation)] <- "midparent"
#' pca_plot(se, color_by = "Generation", shape_by = "Ploidy", add_mean = TRUE)
pca_plot <- function(
        se, PCs = c(1, 2), color_by = NULL, shape_by = NULL, add_mean = FALSE,
        palette = NULL
) {
    
    pc <- paste0("PC", PCs)
    pal <- ppal(palette, "pca")
    
    # Get sample metadata
    cdata <- as.data.frame(colData(se))

    # Perform PCA on vs-transformed and get data frame to plot
    pca_df <- prcomp(t(vst(assay(se))))
    
    # Create final plot data: 1) PCA + coldata; 2) % variance explained
    pdata <- merge(as.data.frame(pca_df$x), cdata, by = "row.names")
    names(pdata)[1] <- "Sample"
    
    varexp <- data.frame(
        row.names = colnames(pca_df$x), 
        Var = round(100 * pca_df$sdev^2 / sum(pca_df$sdev^2), 1)
    )
    
    # (optional) Add mean by a particular variable
    point_mean <- NULL
    if(add_mean) {
        mean_df <- lapply(split(pdata, pdata[[color_by]]), function(x) {
            cmeans <- colMeans(x[, pc])
            df <- data.frame(pcx = cmeans[1], pcy = cmeans[2], g = x[1, color_by])
            names(df) <- c(pc, color_by)
            return(df)
        })
        mean_df <- Reduce(rbind, mean_df)
        point_mean <- geom_point(
            data = mean_df, aes(color = .data[[color_by]]), size = 8, shape = 18
        )
    }
    
    # Plot PCA
    p_pca <- ggplot(pdata, aes(x = .data[[pc[1]]], y = .data[[pc[2]]])) +
        geom_point(
            aes(
                color = .data[[color_by]], shape = .data[[shape_by]]
            ), size = 3, alpha = 0.7
        ) +
        point_mean +
        scale_color_manual(values = pal) +
        labs(
            title = "PCA of samples",
            x = paste0(pc[1], " (", varexp[pc[1], ], "%)"),
            y = paste0(pc[2], " (", varexp[pc[2], ], "%)")
        ) +
        theme_bw() +
        theme(panel.grid = element_blank())
    
    return(p_pca)
}



