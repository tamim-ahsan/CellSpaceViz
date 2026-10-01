#' Example GSEA enrichment results dataset
#'
#' A simulated Gene Set Enrichment Analysis (GSEA) results dataset
#' containing pathway enrichment statistics across multiple Gene
#' Ontology (GO) categories.
#'
#' The dataset is intended for demonstrating and testing
#' `make_gsea_pathway_plot()`.
#'
#' The dataset includes positively enriched pathways, negatively
#' enriched pathways, significant pathways, and non-significant
#' pathways to demonstrate enrichment filtering, pathway ranking,
#' ontology filtering, and colour-scale mapping.
#'
#' @format A data frame with 18 rows and 4 variables:
#' \describe{
#'   \item{Description}{
#'   Gene set or pathway description.
#'   }
#'   \item{NES}{
#'   Normalized enrichment score. Positive values indicate
#'   positively enriched pathways and negative values indicate
#'   negatively enriched pathways.
#'   }
#'   \item{padj}{
#'   Adjusted p-value.
#'   }
#'   \item{Ontology}{
#'   Gene Ontology category annotation.
#'   }
#' }
#'
#' @details
#' The dataset contains pathways from three Gene Ontology (GO)
#' categories:
#' \itemize{
#'   \item GO_BP (Biological Process)
#'   \item GO_MF (Molecular Function)
#'   \item GO_CC (Cellular Component)
#' }
#'
#' Each ontology contains a mixture of:
#' \itemize{
#'   \item positively enriched pathways (NES > 0)
#'   \item negatively enriched pathways (NES < 0)
#'   \item statistically significant pathways (padj < 0.05)
#'   \item non-significant pathways (padj >= 0.05)
#' }
#'
#' The dataset was specifically designed to demonstrate:
#' \itemize{
#'   \item selection of top positively enriched pathways
#'   \item selection of top negatively enriched pathways
#'   \item filtering by significance threshold
#'   \item ontology-specific pathway visualization
#'   \item colour mapping using `-log10(adjusted p-value)`
#' }
#'
#' Example ontology categories include:
#' \itemize{
#'   \item GO_BP pathways such as Oxidative Phosphorylation,
#'   Cell Cycle, and Apoptosis
#'   \item GO_MF pathways such as ATP Binding,
#'   DNA Binding, and Protein Kinase Activity
#'   \item GO_CC pathways such as Mitochondrial Matrix,
#'   Ribosome, and Synapse
#' }
#'
#' The values are simulated and intended solely for package
#' examples, documentation, and automated testing.
#'
#' @source Simulated data.
#'
#' @examples
#' head(gsea_results_df)
#'
"gsea_results_df"
