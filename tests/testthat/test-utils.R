
# Load data ----
data(se_chlamy)
data(deg_list)
data(deg_counts)

partition_table <- expression_partitioning(deg_list)

palette_group <- c(
    `1` = "#5050FFFF", `2` = "#CE3D32FF", `3` = "#749B58FF", `4` = "#F0E685FF", 
    `5` = "#466983FF", `6` = "#BA6338FF", `7` = "#5DB1DDFF", `8` = "#802268FF", 
    `9` = "#6BD76BFF", `10` = "#D595A7FF", `11` = "gold2", `12` = "#837B8DFF"
)

palette_class <- c(
    UP = "#5050FFFF", DOWN = "#CE3D32FF", ADD = "#749B58FF", 
    ELD_P1 = "#F0E685FF", ELD_P2 = "#466983FF"
)

## Create fake assay, colData, and rowData for testing purposes
exp <- matrix(
    rnorm(1000, mean = 10, sd = 2), nrow = 100, ncol = 100, byrow = TRUE
)
rownames(exp) <- paste0("Gene", seq_len(nrow(exp)))
colnames(exp) <- paste0("Sample", seq_len(ncol(exp)))

col_metadata <- data.frame(
    row.names = colnames(exp),
    Class = paste0("Class", sample(1:5, size = ncol(exp), replace = TRUE)),
    Weight = stats::rnorm(ncol(exp), 50, 15)
)

row_metadata <- data.frame(
    row.names = rownames(exp),
    Pathway = paste0("Pathway", sample(1:5, size = nrow(exp), replace = TRUE))
)


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


test_that("get_deg_summary() returns a data frame with DE summary", {
    
    summary_df <- get_deg_summary(deg_list)
    
    expect_equal(ncol(summary_df), 5)
    expect_true(is(summary_df, "data.frame"))
})


test_that("partition_scatterplot() creates a scatterplot", {
    
    p1 <- partition_scatterplot(partition_table, palette_class, "Class")
    p2 <- partition_scatterplot(partition_table, palette_group, "Category")
    
    expect_true(is(p1, "ggplot"))
})


test_that("partition_lineplots() creates a list of line plots", {
    
    p1 <- partition_lineplots(
        partition_table, palette_class, group_by = "Class"
    )
    
    expect_true(is.list(p1))
    expect_true(is(p1[[1]], "ggplot"))
    expect_equal(length(p1), 5)
})


test_that("ppal() generates a vector of color palettes", {
    
    pal1 <- ppal(c("red", "green", "blue"))
    pal2 <- ppal(NULL, "triangle")
    pal3 <- ppal(NULL, "pca")
    pal4 <- ppal(NULL, "partition")
    
    expect_true(is.character(pal1))
    expect_true(is.character(pal2))
    expect_true(is.character(pal3))
    expect_true(is.character(pal4))
    
    expect_error(ppal(NULL, "error"))
})


test_that("metadata2colors() returns a list of metadata and named vectors", {
    
    cols <- metadata2colors(col_metadata)
    
    expect_equal(names(cols), c("metadata", "colors"))
    expect_equal(length(cols), 2)
    
    metadata_4columns <- cbind(
        col_metadata, col_metadata, col_metadata, col_metadata
    )
    expect_error(metadata2colors(metadata_4columns))
})


test_that("se2metadata() returns a list of row and coldata", {
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(exp = exp),
        colData = col_metadata,
        rowData = row_metadata
    )
    m <- se2metadata(se)
    
    expect_equal(length(m), 2)
    expect_equal(names(m), c("rowdata", "coldata"))
})
