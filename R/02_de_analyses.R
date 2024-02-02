
#' Get a table of differential expression expression statistics with \strong{DESeq2}
#'
#' @param se A `SummarizedExperiment` object with a count matrix and sample
#' metadata.
#' @param coldata_column Character indicating the name of column 
#' in \code{colData(se)} where information on the generation are stored. 
#' Default: "Generation".
#' @param parent1 Character indicating which level of the 
#' variable \strong{coldata_column} represents parent 1. Default: "P1".
#' @param parent2 Character indicating which level of the 
#' variable \strong{coldata_column} represents parent 2. Default: "P2".
#' @param offspring Character indicating which level of the 
#' variable \strong{coldata_column} represents the offspring (hybrid
#' or allopolyploid). Default: "F1"
#' @param midparent Character indicating which level of the variable
#' \strong{coldata_column} represents the midparent value. Default:
#' "midparent", as returned by \code{add_midparent_expression()}.
#' @param lfcThreshold Numeric indicating the log2 fold-change threshold
#' to use to consider differentially expressed genes. Default: 0.
#' @param alpha Numeric indicating the adjusted P-value threshold to use
#' to consider differentially expressed genes. Default: 0.01.
#' @param ... Additional arguments to be passed to \code{DESeq2::results()}.
#' 
#' @return A list of data frames with DESeq2's gene-wise tests statistics
#' for each contrast. Each data frame contains the same columns as the
#' output of \code{DESeq2::results()}. Contrasts (list names) are:
#' \describe{
#'   \item{P2_vs_P1}{Parent 2 (numerator) versus parent 1 (denominator).}
#'   \item{F1_vs_P1}{Offspring (numerator) versus parent 1 (denominator).}
#'   \item{F1_vs_P2}{Offspring (numerator) versus parent 2 (denominator).}
#'   \item{F1_vs_midparent}{Offspring (numerator) versus midparent (denominator).}
#' }
#' 
#' The data frame with gene-wise test statistics in each list element contains
#' the following variables:
#' \describe{
#'   \item{baseMean}{Numeric, base mean.}
#'   \item{log2FoldChange}{Numeric, log2-transformed fold changes.}
#'   \item{lfcSE}{Numeric, standard error of the log2-transformed fold changes.}
#'   \item{stat}{Numeric, observed test statistic.}
#'   \item{pvalue}{Numeric, p-value.}
#'   \item{padj}{Numeric, P-value adjusted for multiple testing.}
#' }
#' 
#' The list contains two additional attributes named \strong{ngenes} (numeric,
#' total number of genes), and \strong{plotdata}, which is a 3-column data
#' frame with variables "gene" (character, gene ID), "lFC_F1_vs_P1" (numeric,
#' log2 fold change between F1 and P1), and "lFC_F1_vs_P2" (numeric, 
#' log2 fold change between F1 and P2).
#'
#' @export
#' @rdname get_deg_list
#' @importFrom stats as.formula
#' @importFrom DESeq2 DESeqDataSet estimateSizeFactors sizeFactors
#' @examples
#' data(se_chlamy)
#' se <- add_midparent_expression(se_chlamy)
#' se <- add_size_factors(se, spikein = TRUE)
#' deg_list <- get_deg_list(se)
get_deg_list <- function(
        se, coldata_column = "Generation", 
        parent1 = "P1", parent2 = "P2", offspring = "F1", 
        midparent = "midparent",
        lfcThreshold = 0,
        alpha = 0.01,
        ...
) {
   
    ngenes <- nrow(se)
    
    # Create DESeq object
    form <- as.formula(paste0("~", coldata_column))
    deseq <- suppressMessages(DESeq2::DESeqDataSet(se, design = form))
    
    # Perform DGE
    res <- DESeq2::DESeq(deseq)
    
    # Extract gene-wise test statistics for each contrast
    contrasts <- list(
        "P2_vs_P1" = c(parent2, parent1),
        "F1_vs_P1" = c(offspring, parent1),
        "F1_vs_P2" = c(offspring, parent2),
        "F1_vs_midparent" = c(offspring, midparent)
    )
    
    dge_list <- lapply(contrasts, function(x) {
        res_df <- as.data.frame(DESeq2::results(
            res, contrast = c(coldata_column, x[1], x[2]), ...
        ))
        
        return(res_df)
    })

    # Extract log2foldChange for contrast F1_P1 and F1_P2
    cont <- c("F1_vs_P1", "F1_vs_P2")
    dge_f1 <- dge_list[cont]
    log2fc_df <- lapply(seq_along(dge_f1), function(x) {

        stats_df <- data.frame(
            Gene = rownames(dge_f1[[x]]),
            log2FoldChange = dge_f1[[x]]$log2FoldChange
        )
        names(stats_df) <- c("Gene", paste0("lFC_", cont[x]))
        
        return(stats_df)
    })
    log2fc_df <- Reduce(
        function(x, y) merge(x, y, by = "Gene", all.x = TRUE), 
        log2fc_df
    )
    log2fc_df[is.na(log2fc_df)] <- 0
    
    # Filter DGE list based on user-defined thresholds
    dge_list <- lapply(dge_list, function(x) {
        
        res_df <- x[!is.na(x$padj), ]
        res_df <- res_df[res_df$padj < alpha &
                             abs(res_df$log2FoldChange) > lfcThreshold, ]
        return(res_df)
    })
    attributes(dge_list)$ngenes <- ngenes
    attributes(dge_list)$plotdata <- log2fc_df
    
    return(dge_list)
}


#' Get a count table of differentially expressed genes per contrast
#' 
#' @param deg_list A list of data frames with gene-wise test statistics for
#' differentially expressed genes as returned by \code{get_deg_list()}.
#'
#' @return A data frame with the following variables:
#' \describe{
#'   \item{contrast}{Character, contrast name.}
#'   \item{up}{Numeric, number of up-regulated genes.}
#'   \item{down}{Numeric, number of down-regulated genes.}
#'   \item{total}{Numeric, total number of differentially expressed genes.}
#'   \item{perc_up}{Numeric, percentage of up-regulated genes.}
#'   \item{perc_down}{Numeric, percentage of down-regulated genes.}
#'   \item{perc_total}{Numeric, percentage of diffferentially expressed genes.}
#' }
#'
#' @export
#' @rdname get_deg_counts
#' @examples
#' data(deg_list)
#' deg_counts <- get_deg_counts(deg_list)
get_deg_counts <- function(deg_list) {
    
    # Get number of up- and down-regulated genes per contrast
    count_df <- Reduce(rbind, lapply(seq_along(deg_list), function(x) {
        
        deg_df <- deg_list[[x]]
        deg_df <- deg_df[!is.na(deg_df$log2FoldChange), ]
        
        ## Number of up- and down-regulated genes
        stats_df <- data.frame(
            contrast = names(deg_list)[[x]],
            up = nrow(deg_df[deg_df$log2FoldChange > 0, ]),
            down = nrow(deg_df[deg_df$log2FoldChange < 0, ]),
            total = nrow(deg_df)
        )
        
        ## Add percentage of up- and down-regulated genes
        stats_df$perc_up <- (stats_df$up / attributes(deg_list)$ngenes) * 100
        stats_df$perc_up <- round(stats_df$perc_up, 1)

        stats_df$perc_down <- (stats_df$down / attributes(deg_list)$ngenes) * 100
        stats_df$perc_down <- round(stats_df$perc_down, 1)
        
        stats_df$perc_total <- (stats_df$total / attributes(deg_list)$ngenes) * 100
        stats_df$perc_total <- round(stats_df$perc_total, 1)

        return(stats_df)
    }))
    
    return(count_df)
}




