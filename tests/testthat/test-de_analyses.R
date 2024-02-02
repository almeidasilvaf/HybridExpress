
# Load data ----
data(se_chlamy)
data(deg_list)

## SummarizedExperiment object with midparent expression
se <- add_midparent_expression(se_chlamy)
se <- add_size_factors(se, spikein = TRUE)

# Start tests ----
test_that("get_deg_list() and get_deg_counts() returns list and data frames", {
    
    # 1) Get list of DEGs
    degs <- get_deg_list(se)
    
    # 2) Get data frame of DEG counts
    deg_counts <- get_deg_counts(degs) 
    
    # 1) Tests for get_deg_list()
    expect_equal(length(degs), 4)
    expect_true(is(degs[[1]], "data.frame"))
    
    # 2) Tests for get_deg_counts()
    expect_true(is(deg_counts, "data.frame"))
    expect_equal(ncol(deg_counts), 7)
})

