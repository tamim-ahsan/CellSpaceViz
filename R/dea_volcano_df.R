#' Example differential expression analysis dataset for volcano plots
#'
#' A simulated differential expression analysis dataset containing
#' log2 fold-change and adjusted p-value statistics for 500 genes
#' from a single cell type.
#'
#' The dataset is intended for demonstrating and testing
#' `make_dea_volcano_plot()`.
#'
#' The data were simulated to produce a realistic volcano plot with
#' a mixture of up-regulated, down-regulated, significant, and
#' non-significant genes.
#'
#' @format A data frame with 500 rows and 4 variables:
#' \describe{
#'   \item{cell_type}{
#'   Cell type annotation. All observations belong to
#'   `"Cell-type A"`.
#'   }
#'   \item{gene}{
#'   Simulated gene identifier.
#'   }
#'   \item{log2FoldChange}{
#'   Simulated log2 fold-change value.
#'   }
#'   \item{padj}{
#'   Simulated adjusted p-value.
#'   }
#' }
#'
#' @details
#' The dataset contains 500 genes:
#' \itemize{
#'   \item 100 strongly up-regulated genes with positive log2 fold-changes
#'   and highly significant adjusted p-values.
#'   \item 100 strongly down-regulated genes with negative log2 fold-changes
#'   and moderately significant adjusted p-values.
#'   \item 300 genes with small effect sizes and non-significant
#'   adjusted p-values.
#' }
#'
#' The adjusted p-values were generated to provide a mixture of:
#' \itemize{
#'   \item highly significant genes (`padj < 0.01`)
#'   \item moderately significant genes (`0.01 <= padj < 0.05`)
#'   \item non-significant genes (`padj >= 0.05`)
#' }
#'
#' This structure creates a characteristic volcano plot shape with
#' distinct up-regulated and down-regulated gene populations and a
#' realistic spread of statistical significance values.
#'
#' @source Simulated data.
#'
#' @examples
#' head(dea_volcano_df)
#'
"dea_volcano_df"
