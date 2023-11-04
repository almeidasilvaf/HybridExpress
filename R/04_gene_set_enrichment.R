

#' Perform overrepresentation analysis for a set of genes
#'
#' @param genes Character vector containing genes for which
#' overrepresentation will be tested.
#' @param genesets Named list of character vectors, where list names
#' represent functional categories (e.g., GO, pathway, etc.), and vectors
#' represent their associated genes.
#' @param universe Character vector of genes to be used as universe.
#' @param adj Character indicating the multiple testing correction method
#' as in \code{stats::p.adjust()}.
#'
#' @return A data frame of overrepresentation results with the following
#' variables:
#' \describe{
#'   \item{term}{character, functional term ID/name.}
#'   \item{genes}{numeric, intersection length between input genes and genes
#'                in a particular functional term.}
#'   \item{all}{numeric, number of all genes in a particular functional term.}
#'   \item{pval}{numeric, P-value for the hypergeometric test.}
#'   \item{padj}{numeric, P-value adjusted for multiple comparisons using
#'               the method specified in parameter \strong{adj}.}
#' }
#'
#' @importFrom stats p.adjust phyper
#' @noRd
ora_run <- function(
        genes, genesets, universe, adj = "BH"
) {
    
    # Remove genes that are not in `universe`
    genes <- intersect(genes, universe)
    gene_sets <- lapply(genesets, function(x) intersect(x, universe))
    
    # Define universe and input genes
    n_universe <- length(universe)
    n_genes <- length(genes)
    
    # Define variables of the hypergeometric test
    x <- vapply(gene_sets, function(x) length(intersect(x, genes)), numeric(1))
    m <- vapply(gene_sets, length, numeric(1))
    n <- n_universe - m
    k <- n_genes
    
    # Get P-values from a hypergeometric test
    pvals <- phyper(x - 1, m, n, k, lower.tail = FALSE)
    
    # Store results in a data frame
    results <- data.frame(
        term = names(pvals),
        genes = x,
        all = m,
        pval = pvals,
        row.names = NULL
    )
    
    results$padj <- p.adjust(results$pval, method = adj)
    
    return(results)
}


#' Perform overrepresentation analysis for a set of genes
#'
#' @param genes Character vector containing genes for overrepresentation
#' analysis.
#' @param annotation Annotation data frame with genes in the first column and
#' functional annotation in the other columns. This data frame can be exported
#' from Biomart or similar databases.
#' @param column Column or columns of \strong{annotation} to be used for
#' enrichment.
#' Both character or numeric values with column indices can be used.
#' If users want to supply more than one column, input a character or
#' numeric vector. Default: all columns from \strong{annotation}.
#' @param background Character vector of genes to be used as background
#' for the overrepresentation analysis.
#' @param correction Multiple testing correction method. One of "holm",
#' "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr" or "none".
#' Default is "BH".
#' @param alpha Numeric indicating the adjusted P-value threshold for
#' significance. Default: 0.05.
#' @param min_setsize Numeric indicating the minimum gene set size to be
#' considered. Gene sets correspond to levels of each variable
#' in \strong{annotation}). Default: 5.
#' @param max_setsize Numeric indicating the maximum gene set size to be
#' considered. Gene sets correspond to levels of each variable
#' in \strong{annotation}). Default: 500.
#' @param bp_param BiocParallel back-end to be used.
#' Default: BiocParallel::SerialParam()
#'
#' @return A data frame of overrepresentation results with the following
#' variables:
#' \describe{
#'   \item{term}{Character, functional term ID/name.}
#'   \item{genes}{Numeric, intersection length between input genes and genes
#'                in a particular functional term.}
#'   \item{all}{Numeric, number of all genes in a particular functional term.}
#'   \item{pval}{Numeric, P-value for the hypergeometric test.}
#'   \item{padj}{Numeric, P-value adjusted for multiple comparisons using
#'               the method specified in parameter \strong{adj}.}
#'   \item{category}{Character, name of the grouping variable (i.e., column
#'                   name of \strong{annotation}).}
#' }
#'
#' @rdname ora
#' @export
#' @importFrom BiocParallel bplapply
#' @examples
#' data(se_chlamy)
#' data(go_chlamy)
#' data(deg_list)
#' 
#' # Perform ORA for up-regulated genes in contrast F1_vs_P1
#' up_genes <- deg_list$F1_vs_P1
#' up_genes <- rownames(up_genes[up_genes$log2FoldChange > 0, ])
#' 
#' background <- rownames(se_chlamy)
#' 
#' ora(up_genes, go_chlamy, background = background)
#' 
ora <- function(
        genes, annotation, column = NULL, background, 
        correction = "BH", alpha = 0.05, min_setsize = 5, max_setsize = 500,
        bp_param = BiocParallel::SerialParam()
) {
    
    names(annotation)[1] <- "Gene"
    
    # Filtered `annotation` data frame to keep only specified columns
    col <- names(annotation)
    if(!is.null(column)) {
        col <- if(is.numeric(column)) c(1, column) else c("Gene", column)
    }
    fannot <- annotation[, col]
    
    # Get a data frame of enrichment results for each column
    enrichment <- bplapply(seq_along(fannot)[-1], function(x) {
        
        ## Remove missing values (unannotated genes for a particular category)
        fannot[, x][is.na(fannot[, x])] <- "unannotated"
        
        ## Split variable into a named list
        gene_sets <- split(fannot[, 1], fannot[, x])
        gene_sets <- gene_sets[!duplicated(gene_sets)]
        gene_sets <- gene_sets[names(gene_sets) != "unannotated"]
        
        ## Keep only element with `min_setsize` <= n <= `max_setsize` elements
        set_sizes <- as.numeric(lengths(gene_sets))
        keep <- which(set_sizes >= min_setsize & set_sizes <= max_setsize)
        gene_sets <- gene_sets[keep]
        
        ## Get ORA results and add `category` column
        ora_df <- ora_run(genes, gene_sets, background, correction)
        
        if(nrow(ora_df) > 0) {
            ora_df$category <- names(fannot)[x]
            
            ## Filter non-significant observations out
            ora_df <- ora_df[ora_df$padj <= alpha, ]
        }
        
    }, BPPARAM = bp_param)
    
    enrichment <- Reduce(rbind, enrichment)
    
    return(enrichment)
}