#' Example cell-type count dataset
#'
#' A simulated metadata dataset containing four cell types with
#' differing numbers of observations. The dataset is designed for
#' demonstrating and testing count-based visualizations such as
#' `make_nuclei_count_barplot()`.
#'
#' @format A data frame with 3,350 rows and 1 variable:
#' \describe{
#'   \item{cell_type}{Cell type annotation with four levels:
#'   `"Cell-type A"`, `"Cell-type B"`, `"Cell-type C"`,
#'   and `"Cell-type D"`.}
#' }
#'
#' @details
#' The dataset contains:
#' \itemize{
#'   \item 1,500 observations assigned to Cell-type A
#'   \item 1,000 observations assigned to Cell-type B
#'   \item 600 observations assigned to Cell-type C
#'   \item 250 observations assigned to Cell-type D
#' }
#'
#' @source Simulated data.
#'
#' @examples
#' head(nuclei_count_df)
#'
"nuclei_count_df"
