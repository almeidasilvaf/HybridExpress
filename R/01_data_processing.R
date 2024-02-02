

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


#' Add size factors to normalize count data by library size or by biomass
#'
#' @param se A `SummarizedExperiment` object with a count matrix and sample
#' metadata.
#' @param spikein Logical indicating whether or not to normalize data
#' using spike-ins. If FALSE, data will be normalized by library size.
#' Default: FALSE.
#' @param spikein_pattern Character with the pattern (regex) to use
#' to identify spike-in features in the count matrix. Only valid 
#' if \code{spikein_norm = TRUE}.
#' 
#' @return A `SummarizedExperiment` object as in \strong{se}, but with an extra
#' column in the colData slot named "sizeFactor". This column contains size
#' factors that will be used by DESeq2 when performing differential expression
#' analyses.
#'
#' @importFrom DESeq2 DESeqDataSet estimateSizeFactors sizeFactors
#' @importFrom methods as
#' @export
#' @rdname add_size_factors
#' @examples
#' data(se_chlamy)
#' se_norm <- add_size_factors(se_chlamy)
add_size_factors <- function(
        se, spikein = FALSE, spikein_pattern = "ERCC"
) {
    
    # Create a DESeqDataSet with no design
    deseq <- suppressMessages(DESeq2::DESeqDataSet(se, design = ~1))
    
    # Estimate size factors
    ## Option 1: Using library size
    sf <- estimateSizeFactors(deseq)
    sf <- DESeq2::sizeFactors(sf)
    
    ## Option 2: Using spike-in controls
    if(spikein) {
        se_spikein <- se[grepl(spikein_pattern, rownames(se)), ]
        deseq_spikein <- DESeqDataSet(se_spikein, design = ~1)
        deseq_spikein <- estimateSizeFactors(deseq_spikein)
        
        sf <- sizeFactors(deseq_spikein)
        
        ## Remove rows with spike-in controls (no longer needed)
        deseq <- deseq[!grepl(spikein_pattern, rownames(se)), ]
    }
    
    # Add size factors to colData of DESeqDataSet
    DESeq2::sizeFactors(deseq) <- sf
    
    # Create SummarizedExperiment object from DESeqDataSet
    final_se <- as(deseq, "SummarizedExperiment")
    rownames(final_se) <- rownames(deseq)
    
    return(final_se)
}

