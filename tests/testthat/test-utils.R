
# Load data ----
data(se_chlamy)
se <- add_midparent_expression(se_chlamy)
deg_list <- get_deg_list(se, spikein_norm = TRUE)
deg_counts <- get_deg_counts(deg_list)


# Start tests ----
test_that("get_triangle_graph() returns a list of data frames", {
    
    triangle <- get_triangle_graph()
    
    expect_true(is(triangle, "list"))
    expect_equal(length(triangle), 3)
})


test_that("get_triangle_numbers() returns a list of data frames", {
    
    numbers <- get_triangle_numbers(deg_counts)
    
    expect_true(is(numbers, "list"))
    expect_true(is(numbers[[1]], "data.frame"))
    expect_equal(length(numbers), 2)
})
