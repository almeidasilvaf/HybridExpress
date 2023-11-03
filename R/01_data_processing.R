

#' Add midparent expression to `SummarizedExperiment` object
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
#' @param method Character indicating the method to use to create midparent
#' values. One of 'mean' (default), 'sum', or 'weightedmean'.
#' @param weights Numeric vector of length 2 indicating the weights to give
#' to parents 1 and 2 (respectively) if \code{method == "weightedmean"}.
#' Setting \code{method == "weightedmean"} is used sometimes when parents have
#' different ploidy levels. In such cases, the ploidy levels of parents 1 and 2
#' can be passed in a vector. Default: \code{c(1, 2)}.
#' 
#' @return A `SummarizedExperiment` object.
#'
#' @importFrom SummarizedExperiment colData assay SummarizedExperiment
#' @export
#' @rdname add_midparent_expression
#' @examples
#' data(se_chlamy)
#' new_se <- add_midparent_expression(se_chlamy)
add_midparent_expression <- function(
        se, coldata_column = "Generation", parent1 = "P1", parent2 = "P2",
        method = "mean", weights = c(1, 1)
) {

    # Create a vector with samples from each parent - randomly sampled
    cdata <- as.data.frame(colData(se))
    p1_samples <- sample(rownames(cdata[cdata[, coldata_column] == parent1, ]))
    p2_samples <- sample(rownames(cdata[cdata[, coldata_column] == parent2, ]))
    
    # Get a count matrix by getting the mean of a sample pair (2 parents)
    count <- SummarizedExperiment::assay(se)
    nreplicates <- min(c(length(p1_samples), length(p2_samples)))
    
    count_midparent <- Reduce(cbind, lapply(seq_len(nreplicates), function(x) {
        if(method == "mean") {
            midm <- rowMeans(count[, c(p1_samples[x], p2_samples[x])])
        } else if(method == "sum") {
            midm <- rowSums(count[, c(p1_samples[x], p2_samples[x])])
        } else if(method == "weightedmean") {
            midm <- rowSums(count[, c(p1_samples[x], p2_samples[x])] * weights)
            midm <- midm / sum(weights)
        } else {
            stop("Parameter 'method' must be one of 'mean', 'sum', or 'weightedmean'.")
        }
        
        return(as.matrix(midm))
    }))
    colnames(count_midparent) <- paste0("midparent", seq_len(nreplicates))
    
    # Create a new `SummarizedExperiment` object
    ## 1) Create a new count matrix
    new_counts <- round(cbind(count, count_midparent))
    
    ## 2) Create a new colData
    new_cdata <- data.frame(
        matrix(NA, nrow = ncol(count_midparent), ncol = ncol(cdata))
    )
    colnames(new_cdata) <- names(cdata)
    rownames(new_cdata) <- colnames(count_midparent)
    new_cdata[, coldata_column] <- "midparent"
    
    new_cdata <- rbind(cdata, new_cdata)
    
    ## Create `SummarizedExperiment` object
    new_se <- SummarizedExperiment::SummarizedExperiment(
        assay = list(counts = new_counts),
        colData = new_cdata
    )
    
    return(new_se)
}