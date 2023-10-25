

#' Expression data (in counts) for 3 Chlamydomonas lines (P1, P2, and F1)
#'
#' Two lines (referred to as parent 1 and parent 2) with different ploidy
#' levels were crossed to generate an allopolyploid (F1).
#'
#' @name se_chlamy
#' @format A `SummarizedExperiment` object with an assay (count) and
#' colData.
#' @examples
#' data(se_chlamy)
#' @usage data(se_chlamy)
"se_chlamy"


#' List of differentially expressed genes for all contrasts
#'
#' This object was obtained with \code{get_deg_list()} using the example
#' data set \strong{se_chlamy}.
#' 
#' @name deg_list
#' @format A list of data frames with gene-wise test statistics for
#' differentially expressed genes for each contrast. Contrasts are
#' "P2_vs_P1", "F1_vs_P1", "F1_vs_P2", and "F1_vs_midparent",
#' where the ID before 'vs' represents the numerator, and the ID after 'vs'
#' represents the denominator.
#' @examples 
#' data(deg_list)
#' @usage data(deg_list)
"deg_list"


#' Data frame with frequencies (absolute and relative) of DEGs per contrast
#'
#' This object was obtained with \code{get_deg_counts()} using the example
#' data set \strong{deg_list}.
#' 
#' @name deg_counts
#' @format A data frame with the frequencies (absolute and relative) of
#' up- and down-regulated genes in each contrast. Relative frequencies
#' are calculated relative to the total number of genes in the count matrix
#' used for differential expression analysis.
#' @examples 
#' data(deg_counts)
#' @usage data(deg_counts)
"deg_counts"

