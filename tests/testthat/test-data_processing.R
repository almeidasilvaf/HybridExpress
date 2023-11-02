
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
})

