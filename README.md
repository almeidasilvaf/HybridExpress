
<!-- README.md is generated from README.Rmd. Please edit that file -->

# HybridExpress

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/almeidasilvaf/HybridExpress)](https://github.com/almeidasilvaf/HybridExpress/issues)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![check-bioc](https://github.com/almeidasilvaf/HybridExpress/actions/workflows/rworkflows.devel.yml/badge.svg)](https://github.com/almeidasilvaf/HybridExpress/actions/workflows/rworkflows.devel.yml)
[![Codecov test
coverage](https://codecov.io/gh/almeidasilvaf/HybridExpress/branch/devel/graph/badge.svg)](https://app.codecov.io/gh/almeidasilvaf/HybridExpress?branch=devel)
<!-- badges: end -->

The goal of `HybridExpress` is to perform comparative transcriptomic
analyses of hybrids relative to their progenitor species (a.k.a.
*experimental trios*). The package features:

1.  **Calculation of midparent expression values**, which are *in
    silico* samples obtained from the mean, sum, or weighted mean of
    random sample pairs from each parent;
2.  **Exploratory analyses of sample grouping** with PCA plots and
    heatmaps of hierarchically-clustered pairwise sample correlations;
3.  **Identification of differentially expressed genes** between hybrids
    and their progenitor species, hybrids and midparent values, and the
    two parents. If spike-in standards are available, `HybridExpress`
    uses them to normalize the count data by transcriptome size;
4.  **Classification of genes in expression-based categories and
    classes** based on [Rapp et
    al. (2009)](https://doi.org/10.1186/1741-7007-7-18). The 12
    expression categories proposed by [Rapp et
    al. (2009)](https://doi.org/10.1186/1741-7007-7-18) are grouped into
    5 major classes (transgressive up-regulation, transgressive
    down-regulation, additivity, expression-level dominance toward
    parent 1, and expression-level dominance toward parent 2);
5.  **Functional analyses** through the identification of
    overrepresented functional terms for gene sets of interest.

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `HybridExpress` from
[Bioconductor](http://bioconductor.org/) using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("HybridExpress")
```

And the development version from
[GitHub](https://github.com/almeidasilvaf/HybridExpress) with:

``` r
BiocManager::install("almeidasilvaf/HybridExpress")
```

## Citation

Below is the citation output from using `citation('HybridExpress')` in
R. Please run this yourself to check for any updates on how to cite
**HybridExpress**.

``` r
print(citation('HybridExpress'), bibtex = TRUE)
#> To cite package 'HybridExpress' in publications use:
#> 
#>   Almeida-Silva F, Prost-Boxoen L, Van de Peer Y (2023).
#>   _HybridExpress: Comparative analysis of RNA-seq data for hybrids and
#>   their progenitors_. R package version 0.99.0,
#>   <https://github.com/almeidasilvaf/HybridExpress>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {HybridExpress: Comparative analysis of RNA-seq data for hybrids and their progenitors},
#>     author = {Fabricio Almeida-Silva and Lucas Prost-Boxoen and Yves {Van de Peer}},
#>     year = {2023},
#>     note = {R package version 0.99.0},
#>     url = {https://github.com/almeidasilvaf/HybridExpress},
#>   }
```

Please note that the `HybridExpress` was only made possible thanks to
many other R and bioinformatics software authors, which are cited either
in the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `HybridExpress` project is released with a
[Contributor Code of
Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

- Continuous code testing is possible thanks to [GitHub
  actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
  through *[usethis](https://CRAN.R-project.org/package=usethis)*,
  *[remotes](https://CRAN.R-project.org/package=remotes)*, and
  *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)* customized
  to use [Bioconductor’s docker
  containers](https://www.bioconductor.org/help/docker/) and
  *[BiocCheck](https://bioconductor.org/packages/3.17/BiocCheck)*.
- Code coverage assessment is possible thanks to
  [codecov](https://codecov.io/gh) and
  *[covr](https://CRAN.R-project.org/package=covr)*.
- The [documentation
  website](http://almeidasilvaf.github.io/HybridExpress) is
  automatically updated thanks to
  *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
- The code is styled automatically thanks to
  *[styler](https://CRAN.R-project.org/package=styler)*.
- The documentation is formatted thanks to
  *[devtools](https://CRAN.R-project.org/package=devtools)* and
  *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.17/biocthis)*.
