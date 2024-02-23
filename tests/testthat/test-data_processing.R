
# Load data ----
data(se_chlamy)

# Start tests ----
test_that("add_midparent_expression() adds columns with midparent expression", {
    
    new_se <- add_midparent_expression(se_chlamy)
    new_se2 <- add_midparent_expression(se_chlamy, method = "sum")
    new_se3 <- add_midparent_expression(se_chlamy, method = "weightedmean")
    
    expect_s4_class(new_se, "SummarizedExperiment")
    expect_s4_class(new_se2, "SummarizedExperiment")
    expect_s4_class(new_se3, "SummarizedExperiment")
    expect_true(any(grepl("midparent", colnames(new_se))))
    
    expect_error(add_midparent_expression(se_chlamy, method = "error"))
    expect_error(add_midparent_expression(se_chlamy, coldata_column = "error"))
    expect_error(add_midparent_expression(se_chlamy, parent1 = "error"))
})

test_that("add_size_factors() adds a column named 'sizeFactor' for DESeq2", {
    
    se_norm1 <- add_size_factors(se_chlamy, spikein = FALSE)
    se_norm2 <- add_size_factors(se_chlamy, spikein = TRUE)
    
    expect_true("sizeFactor" %in% colnames(colData(se_norm1)))
    expect_true("sizeFactor" %in% colnames(colData(se_norm2)))
    expect_true(!identical(se_norm1$sizeFactor, se_norm2$sizeFactor))
})
