
<div align="center">

<img src="man/figures/logo.jpg" alt="" width="35%" />

# CellSpaceViz

**Creating plots commonly used in all single-cell and/or spatial omics
studies**

</div>

------------------------------------------------------------------------

# Overview

I created **CellSpaceViz** while generating figures for an article
arising from my PhD project.

This work leveraged a **spatially resolved single-nucleus multiomics
approach** to delineate stress-associated molecular changes in the
**prefrontal cortex** of psychiatric patients.

Throughout the project, I repeatedly found myself writing the same
visualization code for:

- dimensionality reduction plots
- patterns of marker gene expressions
- cell-type composition summaries
- differential expression analyses
- pathway enrichment analyses

Eventually, these functions became reusable enough that I decided to
package them into **CellSpaceViz**.

The goal of the package is to simplify the generation of
publication-ready visualizations commonly used in:

- single-cell RNA-seq
- single-nucleus RNA-seq
- spatial transcriptomics
- multiomics studies
- differential expression analyses

------------------------------------------------------------------------

# Installation

The package can be installed directly from GitHub.

First, install **remotes** only if it is not already available:

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
```

Then install the package:

``` r
remotes::install_github(
  "tamim-ahsan/CellSpaceViz"
)
```

Load the package:

``` r
library(CellSpaceViz)
```

------------------------------------------------------------------------

# Example Datasets Included in the Package

| Dataset           | Purpose                        |
|-------------------|--------------------------------|
| `reduced_dim_df`  | Dimensionality reduction plots |
| `marker_df`       | Marker-gene dot plots          |
| `nuclei_count_df` | Cell or nuclei count plots     |
| `deg_summary_df`  | DEG count summaries            |
| `deg_results_df`  | DEA metric plots               |
| `dea_volcano_df`  | Volcano plots                  |
| `dea_heatmap_df`  | LFC heatmaps                   |
| `gsea_results_df` | GSEA pathway enrichment plots  |

------------------------------------------------------------------------

# Learning More

Every function and dataset is documented.

Use:

``` r
?function_name
```

or

``` r
?dataset_name
```

to access the help pages.

Examples:

``` r
?make_marker_dotplot
```

``` r
?make_dea_volcano_plot
```

``` r
?dea_heatmap_df
```

``` r
?gsea_results_df
```

------------------------------------------------------------------------

# Included Functions and Example Datasets

Each plotting function is accompanied by a small example dataset
designed to demonstrate its functionality and simplify testing.

For additional details, users can access the documentation using:

``` r
?function_name
```

or

``` r
?dataset_name
```

For example:

``` r
?make_dea_volcano_plot
```

or

``` r
?dea_volcano_df
```

------------------------------------------------------------------------

### 1. `make_reduced_dim_plots()`

Creates dimensionality reduction visualizations such as:

- UMAP
- t-SNE
- PCA

colored and faceted by grouping variables.

**Example dataset**

``` r
reduced_dim_df
```

Documentation:

``` r
?make_reduced_dim_plots
?reduced_dim_df
```

------------------------------------------------------------------------

### 2. `make_marker_dotplot()`

Creates marker-gene dot plots displaying:

- average expression
- fraction of expressing cells

across cell types.

**Example dataset**

``` r
marker_df
```

Documentation:

``` r
?make_marker_dotplot
?marker_df
```

------------------------------------------------------------------------

### 3. `make_nuclei_count_barplot()`

Creates horizontal bar plots summarizing:

- cell counts
- nuclei counts
- spot counts

across cell types.

**Example dataset**

``` r
nuclei_count_df
```

Documentation:

``` r
?make_nuclei_count_barplot
?nuclei_count_df
```

------------------------------------------------------------------------

### 4. `make_deg_count_plot()`

Creates stacked bar plots showing:

- up-regulated genes
- down-regulated genes
- total DEGs

for each cell type.

**Example dataset**

``` r
deg_summary_df
```

Documentation:

``` r
?make_deg_count_plot
?deg_summary_df
```

------------------------------------------------------------------------

### 5. `make_dea_metric_plot()`

Creates differential expression metric bubble plots where:

- x-axis represents log2 fold-change
- bubble size represents significance

using:

``` text
-log10(adjusted p-value)
```

**Example dataset**

``` r
deg_results_df
```

Documentation:

``` r
?make_dea_metric_plot
?deg_results_df
```

------------------------------------------------------------------------

### 6. `make_dea_volcano_plot()`

Creates publication-style volcano plots with:

- significance color coding
- automatic gene labeling
- customizable significance thresholds

**Example dataset**

``` r
dea_volcano_df
```

Documentation:

``` r
?make_dea_volcano_plot
?dea_volcano_df
```

------------------------------------------------------------------------

### 7. `make_dea_lfc_heatmap()`

Creates heatmaps displaying:

- log2 fold-change values
- significance annotations

across cell types or experimental conditions.

Significance is represented by:

``` text
*
**
***
```

annotations.

**Example dataset**

``` r
dea_heatmap_df
```

Documentation:

``` r
?make_dea_lfc_heatmap
?dea_heatmap_df
```

------------------------------------------------------------------------

### 8. `make_gsea_pathway_plot()`

Creates pathway enrichment visualizations from Gene Set Enrichment
Analysis (GSEA).

Features include:

- pathway ranking by NES
- positive and negative enrichment visualization
- ontology filtering
- significance-based color scaling

**Example dataset**

``` r
gsea_results_df
```

Documentation:

``` r
?make_gsea_pathway_plot
?gsea_results_df
```

------------------------------------------------------------------------

# Future Plans

This README provides only a brief introduction to the package.

I plan to create a more detailed documentation/tutorial soon showing:

- how to use the functions and
- what the plots look like

Stay tuned!

------------------------------------------------------------------------

# A Final (and Honest) Note 😄

This is my **first attempt at writing an R package**.

CellSpaceViz came together through a combination of:

- reading a [book](https://r-pkgs.org/) on R package development
- watching a lot of YouTube tutorials
- substantial trial and error
- many conversations with AI assistants, particularly **Microsoft
  Copilot** and **Google Gemini**

It has been a surprisingly fun learning experience, and I hope the
package proves useful to others working with similar datasets.

If CellSpaceViz saves somebody from repeatedly rewriting the same
plotting code during manuscript preparation, then it has already
achieved its goal 🚀.
