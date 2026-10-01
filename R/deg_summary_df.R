#' Example differential expression summary dataset
#'
#' A simulated differential expression summary dataset containing
#' counts of up-regulated, down-regulated, and total differentially
#' expressed genes (DEGs) across four cell types.
#'
#' The dataset is intended for demonstrating and testing
#' `make_deg_count_plot()`.
#'
#' @format A data frame with 4 rows and 4 variables:
#' \describe{
#'   \item{cell_type}{
#'   Cell type annotation with four levels:
#'   `"Cell-type A"`, `"Cell-type B"`,
#'   `"Cell-type C"`, and `"Cell-type D"`.
#'   }
#'   \item{up_regulated}{
#'   Number of significantly up-regulated genes.
#'   }
#'   \item{down_regulated}{
#'   Number of significantly down-regulated genes.
#'   }
#'   \item{total_degs}{
#'   Total number of differentially expressed genes
#'   (`up_regulated + down_regulated`).
#'   }
#' }
#'
#' @details
#' Cell-type A exhibits the largest number of differentially
#' expressed genes, while Cell-type D exhibits the smallest.
#' The dataset was generated specifically for package examples,
#' documentation, and automated testing of DEG visualization
#' functions.
#'
#' The counts are:
#' \itemize{
#'   \item Cell-type A: 450 up-regulated and 250 down-regulated genes (700 total DEGs)
#'   \item Cell-type B: 320 up-regulated and 280 down-regulated genes (600 total DEGs)
#'   \item Cell-type C: 180 up-regulated and 120 down-regulated genes (300 total DEGs)
#'   \item Cell-type D: 90 up-regulated and 60 down-regulated genes (150 total DEGs)
#' }
#'
#' @source Simulated data.
#'
#' @examples
#' head(deg_summary_df)
#'
"deg_summary_df"
