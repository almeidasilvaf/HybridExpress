
# Load data ----
data(deg_counts)

# Start tests ----
test_that("plot_expression_triangle() returns a ggplot object", {
    
    p <- plot_expression_triangle(deg_counts)
    
    expect_true(is(p, "ggplot"))
})
