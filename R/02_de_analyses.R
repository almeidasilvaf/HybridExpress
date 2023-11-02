
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
#' @param spikein_norm Logical indicating whether or not to normalize data
#' using spike-ins. Default: FALSE.
#' @param spikein_pattern Character with the pattern (regex) to use
#' to identify spike-in features in the count matrix. Only valid 
#' if \code{spikein_norm = TRUE}.
#' @param lfcThreshold Numeric indicating the log2 fold-change threshold
#' to use to consider differentially expressed genes. Default: 0.
#' @param alpha Numeric indicating the adjusted P-value threshold to use
#' to consider differentially expressed genes. Default: 0.05.
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
#' @export
#' @rdname get_deg_list
#' @importFrom stats as.formula
#' @importFrom DESeq2 DESeqDataSet estimateSizeFactors sizeFactors
#' @examples
#' data(se_chlamy)
#' se <- add_midparent_expression(se_chlamy)
#' deg_list <- get_deg_list(se, spikein_norm = TRUE)
get_deg_list <- function(
        se, coldata_column = "Generation", 
        parent1 = "P1", parent2 = "P2", offspring = "F1",
        spikein_norm = FALSE, spikein_pattern = "ERCC",
        lfcThreshold = 0,
        alpha = 0.05,
        ...
) {
   
    ngenes <- nrow(se)
    
    # Create DESeq object
    form <- as.formula(paste0("~", coldata_column))
    deseq <- suppressMessages(DESeq2::DESeqDataSet(se, design = form))
    if(spikein_norm) {
        spikein_rows <- grepl(spikein_pattern, rownames(se))
        spikein_se <- se[spikein_rows, ]
        spikein_deseq <- DESeqDataSet(spikein_se, design = form)
        spikein_deseq <- estimateSizeFactors(spikein_deseq)
        
        spikein_size_factors <- sizeFactors(spikein_deseq)
        
        deseq <- deseq[!spikein_rows, ]
        DESeq2::sizeFactors(deseq) <- spikein_size_factors
        ngenes <- nrow(deseq)
    }
    
    # Perform DGE
    res <- DESeq2::DESeq(deseq)
    
    # Extract gene-wise test statistics for each contrast
    contrasts <- list(
        "P2_vs_P1" = c(parent2, parent1),
        "F1_vs_P1" = c(offspring, parent1),
        "F1_vs_P2" = c(offspring, parent2),
        "F1_vs_midparent" = c(offspring, "midparent")
    )
    
    dge_list <- lapply(contrasts, function(x) {
        res_df <- as.data.frame(DESeq2::results(
            res, contrast = c(coldata_column, x[1], x[2]), ...
        ))
        res_df <- res_df[res_df$padj < alpha & 
                             abs(res_df$log2FoldChange) >= lfcThreshold, ]
        
        return(res_df)
    })
    
    attributes(dge_list)$ngenes <- ngenes
    
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




