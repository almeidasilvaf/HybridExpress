
# Load data ----
data(se_chlamy)
data(go_chlamy)
data(deg_list)
 
# Perform ORA for up-regulated genes in contrast F1_vs_P1
up_genes <- deg_list$F1_vs_P1
up_genes <- rownames(up_genes[up_genes$log2FoldChange > 0, ])

background <- rownames(se_chlamy)

# Start tests ----
test_that("ora() performs overrepresentation analysis", {
    
    ora_df <- ora(up_genes, go_chlamy, background = background, column = 2)
    
    expect_true(is(ora_df, "data.frame"))
})
