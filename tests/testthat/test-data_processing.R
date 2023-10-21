
# Load data ----
data(se_chlamy)


# Start tests ----
test_that("add_midparent_expression() adds columns with midparent expression", {
    
    new_se <- add_midparent_expression(se_chlamy)
    
    expect_s4_class(new_se, "SummarizedExperiment")
    expect_true(any(grepl("midparent", colnames(new_se))))
})

