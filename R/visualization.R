
#' Plot a triangle of comparisons of DEG sets among generations
#'
#' @param deg_counts Data frame with number of differentially expressed
#' genes per contrast as returned by \code{get_deg_counts}.
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
plot_expression_triangle <- function(deg_counts) {
    
    # Get plot data
    tgraph <- get_triangle_graph()
    numbers <- get_triangle_numbers(deg_counts)
    
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
            color = "grey20", 
            fill = c("dodgerblue3", "firebrick", "purple4", "mediumorchid4")
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




