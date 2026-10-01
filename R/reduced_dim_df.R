#' Sample dataset with reduced dimensions
#'
#' A simulated dataset containing UMAP coordinates, cell type labels, and multi-omics
#' modality information representing 200 single nuclei/cells for testing and demonstration purposes.
#'
#' @format A data frame with 200 rows and 4 variables:
#' \describe{
#'   \item{UMAP1}{Numeric, first UMAP coordinate}
#'   \item{UMAP2}{Numeric, second UMAP coordinate}
#'   \item{cell_type}{Factor, cell type classification (T cell, B cell, NK cell, Monocyte)}
#'   \item{modality}{Character, assay modality}
#' }
#' @source Simulated synthetic data for CellSpaceViz.
#' @examples
#' data(reduced_dim_df)
#'
"reduced_dim_df"
