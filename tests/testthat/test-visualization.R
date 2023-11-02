
# Load data ----
data(deg_counts)
data(deg_list)
data(se_chlamy)

se <- add_midparent_expression(se_chlamy)
ptable <- expression_partitioning(deg_list)
se$Ploidy[is.na(se$Ploidy)] <- "midparent"
se$Generation[is.na(se$Generation)] <- "midparent"


# Start tests ----
test_that("plot_expression_triangle() returns a ggplot object", {
    
    p <- plot_expression_triangle(
        deg_counts, 
        box_labels = c("Parent 1", "Parent 2", "Hybrid", "Midparent")
    )
    
    expect_true(is(p, "ggplot"))
})


test_that("plot_expression_partitions() returns a multi-panel ggplot", {
    
    p1 <- plot_expression_partitions(ptable, group_by = "Category")
    p2 <- plot_expression_partitions(ptable, group_by = "Class")
    
    expect_true(is(p1, "ggplot"))
    expect_true(is(p2, "ggplot"))
})

test_that("plot_partition_frequencies() returns a multi-panel ggplot", {
    
    p1 <- plot_partition_frequencies(ptable, group_by = "Category")
    p2 <- plot_partition_frequencies(ptable, group_by = "Class")
    
    expect_true(is(p1, "ggplot"))
    expect_true(is(p2, "ggplot"))
})

test_that("pca_plot() returns a ggplot object", {
    
    p1 <- pca_plot(
        se, color_by = "Ploidy", shape_by = "Generation", add_mean = TRUE
    )
    
    expect_true(is(p1, "ggplot"))
})


test_that("plot_samplecor() plots a heatmap of pairwise correlations", {
    
    p <- plot_samplecor(se, ntop = 50)
    
    expect_equal(attr(class(p), "package"), "ComplexHeatmap")
})




