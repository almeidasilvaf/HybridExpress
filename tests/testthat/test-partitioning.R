
# Load data ----
data(deg_list)


# Start tests ----
test_that("expression_partitioning() returns groups", {
    
    parts <- expression_partitioning(deg_list)
    
    expect_equal(ncol(parts), 5)
    expect_true(is(parts, "data.frame"))
})
