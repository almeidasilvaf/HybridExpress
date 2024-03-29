Data acquisition
================

# Overview

Here, we describe all code we used to obtain the example data sets in
this package.

## `data/se_chlamy.rda`

``` r
library(tidyverse)

coldata_path <- "~/meta_data.csv"
counts_path <- "~/TS_count.txt"

# Create colData
coldata <- read_csv(coldata_path, show_col_types = FALSE) |>
    tibble::column_to_rownames("sample") |>
    mutate(ploidy_gen = str_c(ploidy, generation, sep = "_")) |>
    filter(ploidy_gen %in% c("haploid_Anc", "diploid_Anc", "triploid_Anc")) |>
    select(Ploidy = ploidy) |>
    mutate(Generation = case_when(
        Ploidy == "diploid" ~ "P1",
        Ploidy == "haploid" ~ "P2",
        Ploidy == "triploid" ~ "F1"
    )) |>
    mutate(Generation = factor(Generation, levels = c("F1", "P1", "P2")))

# Create assay
counts <- read_tsv(counts_path, show_col_types = FALSE) |>
    select(-transcript_id) |>
    filter(
        !gene_id %in% c(
            "__alignment_not_unique", "__ambiguous",
            "__no_feature", "__not_aligned", "__too_low_aQual"
        )
    ) |>
    mutate(gene_id = str_replace_all(gene_id, "_4532.v6.1", "")) |>
    tibble::column_to_rownames("gene_id") |>
    as.matrix()


# Filter the expression matrix
counts <- counts[, rownames(coldata)]

# Rename samples
rownames(coldata) <- paste0("S", seq_len(ncol(counts)))
colnames(counts) <- paste0("S", seq_len(ncol(counts)))

# Keep only genes with rowSums > 10 + spike-ins
spikeins <- grep("ERCC", rownames(counts))
expressed <- which(rowSums(counts) > 10)

fcount <- counts[sort(unique(c(expressed, spikeins))), ]

# Create `SummarizedExperiment` object
se_chlamy <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = fcount),
    colData = coldata
)

# Save object
usethis::use_data(se_chlamy, compress = "xz", overwrite = TRUE)
```

## `data/deg_list.rda`

This object contains a list of differentially expressed genes as
returned by `get_deg_list()`.

``` r
# Load data and add midparent expression
data(se_chlamy)
se <- add_midparent_expression(se_chlamy)

# Get DEG list
deg_list <- get_deg_list(se, spikein_norm = TRUE)

# Save object
usethis::use_data(deg_list, compress = "xz", overwrite = TRUE)
```

## `data/deg_counts.rda`

This object contains a list of DEG counts for all contrasts as returned
by `get_deg_counts()`.

``` r
# Load data and add midparent expression
deg_counts <- get_deg_counts(deg_list)

# Save data
usethis::use_data(deg_counts, compress = "xz", overwrite = TRUE)
```

## `data/`

``` r
# Read table with functional annotation (from Phytozome)
annotation <- readr::read_tsv("~/functional_annotation.txt.gz")

# Get a table with genes and GO IDs
go_table <- annotation |>
    select(gene = locusName, GO) |>
    mutate(gene = str_replace_all(gene, "_4532", "")) |>
    separate_longer_delim(GO, delim = " ") |>
    dplyr::distinct()

# Replace GO IDs with names
go_plaza <- readr::read_tsv(
    "https://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_05/GO/go.cre.csv.gz",
    skip = 8, show_col_types = FALSE
) |>
    select(GO = go, description) |>
    distinct()

go_chlamy <- go_table |>
    left_join(go_plaza) |>
    select(gene, GO = description) |>
    as.data.frame()

# Save data set
usethis::use_data(go_chlamy, compress = "xz")
```
