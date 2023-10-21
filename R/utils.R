
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
        contrast = c("P1_vs_P2", "P1_vs_F1", "P2_vs_F1", "midparent_vs_F1"),
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
            "P1_vs_P2-P1", "P1_vs_P2-P1-P2", 
            "P1_vs_F1-P1", "P1_vs_F1-F1",
            "P2_vs_F1-P2", "P2_vs_F1-F1",
            "midparent_vs_F1-midparent", "midparent_vs_F1-F1"
        ),
        n_up = c(
            deg_counts$up[deg_counts$contrast == "P1_vs_P2"],
            deg_counts$down[deg_counts$contrast == "P1_vs_P2"],
            deg_counts$up[deg_counts$contrast == "P1_vs_F1"],
            deg_counts$down[deg_counts$contrast == "P1_vs_F1"],
            deg_counts$up[deg_counts$contrast == "P2_vs_F1"],
            deg_counts$down[deg_counts$contrast == "P2_vs_F1"],
            deg_counts$up[deg_counts$contrast == "midparent_vs_F1"],
            deg_counts$down[deg_counts$contrast == "midparent_vs_F1"]
        ),
        perc_up = c(
            deg_counts$perc_up[deg_counts$contrast == "P1_vs_P2"],
            deg_counts$perc_down[deg_counts$contrast == "P1_vs_P2"],
            deg_counts$perc_up[deg_counts$contrast == "P1_vs_F1"],
            deg_counts$perc_down[deg_counts$contrast == "P1_vs_F1"],
            deg_counts$perc_up[deg_counts$contrast == "P2_vs_F1"],
            deg_counts$perc_down[deg_counts$contrast == "P2_vs_F1"],
            deg_counts$perc_up[deg_counts$contrast == "midparent_vs_F1"],
            deg_counts$perc_down[deg_counts$contrast == "midparent_vs_F1"]
        ),
        x = c(
            0.55, 3.45, # P1-P2
            -0.35, 1.40, # P1-F1
            4.35, 2.60, # P2-F1
            1.75, 1.75 # midparent-F1
        ),
        y = c(
            -0.25, -0.25, # P1-P2
            0.55, 4, # P1-F1
            0.55, 4, # P2-F1
            2.05, 3.45 # mid-parent-F1
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




