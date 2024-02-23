
#' Check if a column exists in the colData of a `SummarizedExperiment` object
#'
#' @param se A `SummarizedExperiment` object.
#' @param column Character with name of the column to check.
#' 
#' @return TRUE if the column exists, and ERROR otherwise.
#'
#' @importFrom SummarizedExperiment colData
#' @noRd
#' @examples
#' data(se_chlamy)
#' check_coldata_column(se_chlamy, "Generation")
check_coldata_column <- function(se, column) {
    
    cdata <- colData(se)
    if(!column %in% names(cdata)) {
        stop("Column '", column, "' is not present in the `colData` slot.")
    }
    
    return(TRUE)
}


#' Check if levels exist in a colData column of a `SummarizedExperiment` object
#'
#' @param se A `SummarizedExperiment` object.
#' @param column Character with name of the column where levels are.
#' @param levels Character with levels to check for presence in \strong{column}.
#' 
#' @return TRUE if the column exists, and ERROR otherwise.
#'
#' @importFrom SummarizedExperiment colData
#' @noRd
#' @examples
#' data(se_chlamy)
#' check_coldata_levels(se_chlamy, "Generation", levels = c("P1", "P2"))
check_coldata_levels <- function(se, column, levels) {
    
    col <- unique(colData(se)[[column]])
    
    if(any(levels %in% col == FALSE)) {
        stop(
            "All levels (", paste0(levels, collapse = ","), 
            ") must be in column '", column, "'."
        )
    }
    
    return(TRUE)
}

