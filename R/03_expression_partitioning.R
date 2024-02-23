
#' Partition genes in groups based on their expression patterns
#' 
#' @param deg_list A list of data frames with gene-wise test statistics for
#' differentially expressed genes as returned by \code{get_deg_list()}.
#'
#' @return A data with the following variables:
#' \describe{
#'   \item{Gene}{Character, gene ID.}
#'   \item{Category}{Factor, expression group. Category names are numbers from
#'   1 to 12.}
#'   \item{Class}{Factor, expression group class. One of "UP" (transgressive
#'   up-regulation), "DOWN" (transgressive down-regulation), 
#'   "ADD" (additivity), "ELD_P1" (expression-level dominance toward
#'   the parent 1), or "ELD_P2" (expression-level dominance toward
#'   the parent 2).}
#' }
#'
#' @rdname expression_partitioning
#' @export
#' @examples
#' data(deg_list)
#' exp_partitions <- expression_partitioning(deg_list)
expression_partitioning <- function(deg_list) {
    
    # Get DEG summary table
    ds <- get_deg_summary(deg_list)
    
    # Create a data frame with genes concatenated text for contrasts
    genes_and_deg_string <- data.frame(
        Gene = ds$Gene,
        DEG = paste(ds$F1_vs_P2, ds$F1_vs_P1, ds$P2_vs_P1, sep = "-")
    )
    
    # Data frame with classification rules for groups and classes
    cl <- c(
        "ADD", "ELD_P2", "DOWN", "ELD_P1", "UP", "UP", 
        "DOWN", "UP", "ELD_P1", "DOWN", "ELD_P2", "ADD"
    )
    partition_rules <- data.frame(
        Category = factor(seq_len(12), levels = seq_len(12)),
        Class = factor(cl, levels = c("UP", "DOWN", "ADD", "ELD_P1", "ELD_P2")),
        DEG = c(
            "down-up-up", #             1
            "unchanged-up-up", #        2
            "down-down-up", #           3
            "up-unchanged-down", #      4
            "up-up-up", #               5
            "up-up-down", #             6
            "down-down-unchanged", #    7
            "up-up-unchanged", #        8
            "down-unchanged-up", #      9
            "down-down-down", #         10
            "unchanged-down-down", #    11
            "up-down-down" #            12
        )
    )

    # Classify genes in 12 groups
    classified_genes <- merge(
        genes_and_deg_string, partition_rules, all.x = TRUE
    )
    
    # Get a table with log2foldChange for contrasts `F1_vs_P1` and `F1_vs_P2`
    log2fc_df <- attributes(deg_list)$plotdata
    
    # Combine classification data frame with log2foldchange
    class_df <- merge(
        classified_genes[, c("Gene", "Category", "Class")], log2fc_df,
        all.x = TRUE
    )
    class_df <- class_df[!is.na(class_df$Category), ]
    class_df <- class_df[order(class_df$Category), ]
    
    rownames(class_df) <- NULL
    
    return(class_df)
}

