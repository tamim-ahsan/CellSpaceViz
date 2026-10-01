#' Example differential expression results dataset
#'
#' A simulated differential expression dataset containing
#' log2 fold-change and adjusted p-value statistics for
#' 400 genes across four cell types.
#'
#' The dataset is intended for demonstrating and testing
#' `make_deg_metric_plot()`.
#'
#' Each cell type contains 100 genes with a mixture of:
#' \itemize{
#'   \item strongly significant genes (`padj < 0.01`)
#'   \item moderately significant genes (`0.01 <= padj < 0.05`)
#'   \item non-significant genes (`padj >= 0.05`)
#' }
#'
#' This produces realistic variation in bubble sizes when visualising
#' `-log10(adjusted p-value)` and a broad distribution of positive,
#' negative, and near-zero log2 fold-change values.
#'
#' @format A data frame with 400 rows and 4 variables:
#' \describe{
#'   \item{cell_type}{
#'   Cell type annotation with four levels:
#'   `"Cell-type A"`, `"Cell-type B"`,
#'   `"Cell-type C"`, and `"Cell-type D"`.
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
#' The dataset contains:
#' \itemize{
#'   \item 100 genes assigned to Cell-type A
#'   \item 100 genes assigned to Cell-type B
#'   \item 100 genes assigned to Cell-type C
#'   \item 100 genes assigned to Cell-type D
#' }
#'
#' Within each cell type:
#' \itemize{
#'   \item 40 genes are highly significant (`padj < 0.01`)
#'   \item 30 genes are moderately significant (`0.01 <= padj < 0.05`)
#'   \item 30 genes are non-significant (`padj >= 0.05`)
#' }
#'
#' The values are simulated and intended solely for package
#' documentation, examples, and automated testing.
#'
#' @source Simulated data.
#'
#' @examples
#' head(deg_results_df)
#'
"deg_results_df"
