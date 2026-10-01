test_that("example dataset is available and correctly formatted", {

  expect_true(exists("gsea_results_df"))

  required_cols <- c(
    "Description",
    "NES",
    "padj",
    "Ontology"
  )

  expect_true(
    all(required_cols %in% colnames(gsea_results_df))
  )

})

test_that("returns a ggplot object", {

  p <- make_gsea_pathway_plot(
    df = gsea_results_df,
    enrichment_score_col = "NES",
    padj_col = "padj",
    description_col = "Description"
  )

  expect_s3_class(p, "ggplot")

})

test_that("plot title works", {

  p <- make_gsea_pathway_plot(
    df = gsea_results_df,
    enrichment_score_col = "NES",
    padj_col = "padj",
    description_col = "Description",
    plot_title = "Example GSEA Plot"
  )

  expect_s3_class(p, "ggplot")

})

test_that("ontology filtering works for GO_BP", {

  p <- make_gsea_pathway_plot(
    df = gsea_results_df,
    enrichment_score_col = "NES",
    padj_col = "padj",
    description_col = "Description",
    ontology_col = "Ontology",
    ontology = "GO_BP"
  )

  expect_s3_class(p, "ggplot")

})

test_that("ontology filtering works for GO_MF", {

  p <- make_gsea_pathway_plot(
    df = gsea_results_df,
    enrichment_score_col = "NES",
    padj_col = "padj",
    description_col = "Description",
    ontology_col = "Ontology",
    ontology = "GO_MF"
  )

  expect_s3_class(p, "ggplot")

})

test_that("ontology filtering works for GO_CC", {

  p <- make_gsea_pathway_plot(
    df = gsea_results_df,
    enrichment_score_col = "NES",
    padj_col = "padj",
    description_col = "Description",
    ontology_col = "Ontology",
    ontology = "GO_CC"
  )

  expect_s3_class(p, "ggplot")

})

test_that("all supported viridis palettes work", {

  palettes <- c(
    "magma",
    "inferno",
    "plasma",
    "viridis",
    "cividis",
    "rocket",
    "mako",
    "turbo"
  )

  for (pal in palettes) {

    p <- make_gsea_pathway_plot(
      df = gsea_results_df,
      enrichment_score_col = "NES",
      padj_col = "padj",
      description_col = "Description",
      fill_scale_option = pal
    )

    expect_s3_class(p, "ggplot")

  }

})

test_that("missing enrichment score column throws error", {

  expect_error(

    make_gsea_pathway_plot(
      df = gsea_results_df,
      enrichment_score_col = "bad_column",
      padj_col = "padj",
      description_col = "Description"
    )

  )

})

test_that("missing padj column throws error", {

  expect_error(

    make_gsea_pathway_plot(
      df = gsea_results_df,
      enrichment_score_col = "NES",
      padj_col = "bad_column",
      description_col = "Description"
    )

  )

})

test_that("missing description column throws error", {

  expect_error(

    make_gsea_pathway_plot(
      df = gsea_results_df,
      enrichment_score_col = "NES",
      padj_col = "padj",
      description_col = "bad_column"
    )

  )

})

test_that("non-numeric NES column throws error", {

  test_df <- gsea_results_df

  test_df$NES <- as.character(
    test_df$NES
  )

  expect_error(

    make_gsea_pathway_plot(
      df = test_df,
      enrichment_score_col = "NES",
      padj_col = "padj",
      description_col = "Description"
    )

  )

})

test_that("non-numeric padj column throws error", {

  test_df <- gsea_results_df

  test_df$padj <- as.character(
    test_df$padj
  )

  expect_error(

    make_gsea_pathway_plot(
      df = test_df,
      enrichment_score_col = "NES",
      padj_col = "padj",
      description_col = "Description"
    )

  )

})

test_that("negative adjusted p-values throw error", {

  test_df <- gsea_results_df

  test_df$padj[1] <- -0.01

  expect_error(

    make_gsea_pathway_plot(
      df = test_df,
      enrichment_score_col = "NES",
      padj_col = "padj",
      description_col = "Description"
    )

  )

})

test_that("invalid ontology specification throws error", {

  expect_error(

    make_gsea_pathway_plot(
      df = gsea_results_df,
      enrichment_score_col = "NES",
      padj_col = "padj",
      description_col = "Description",
      ontology = "GO_BP"
    )

  )

})

test_that("invalid ontology column throws error", {

  expect_error(

    make_gsea_pathway_plot(
      df = gsea_results_df,
      enrichment_score_col = "NES",
      padj_col = "padj",
      description_col = "Description",
      ontology_col = "bad_column",
      ontology = "GO_BP"
    )

  )

})

test_that("invalid viridis palette throws error", {

  expect_error(

    make_gsea_pathway_plot(
      df = gsea_results_df,
      enrichment_score_col = "NES",
      padj_col = "padj",
      description_col = "Description",
      fill_scale_option = "rainbow"
    )

  )

})

test_that("invalid padj threshold throws error", {

  expect_error(

    make_gsea_pathway_plot(
      df = gsea_results_df,
      enrichment_score_col = "NES",
      padj_col = "padj",
      description_col = "Description",
      padj_threshold = -0.01
    )

  )

  expect_error(

    make_gsea_pathway_plot(
      df = gsea_results_df,
      enrichment_score_col = "NES",
      padj_col = "padj",
      description_col = "Description",
      padj_threshold = 1.5
    )

  )

})

test_that("invalid pathway count throws error", {

  expect_error(

    make_gsea_pathway_plot(
      df = gsea_results_df,
      enrichment_score_col = "NES",
      padj_col = "padj",
      description_col = "Description",
      n_pathways_to_show = 0
    )

  )

})

test_that("strict significance threshold throws error when no pathways pass", {

  expect_error(

    make_gsea_pathway_plot(
      df = gsea_results_df,
      enrichment_score_col = "NES",
      padj_col = "padj",
      description_col = "Description",
      padj_threshold = 1e-20
    )

  )

})

test_that("bar layer is generated", {

  p <- make_gsea_pathway_plot(
    df = gsea_results_df,
    enrichment_score_col = "NES",
    padj_col = "padj",
    description_col = "Description"
  )

  bar_layer <- ggplot2::ggplot_build(p)$data[[1]]

  expect_true(
    nrow(bar_layer) > 0
  )

})

test_that("all displayed pathways are plotted", {

  p <- make_gsea_pathway_plot(
    df = gsea_results_df,
    enrichment_score_col = "NES",
    padj_col = "padj",
    description_col = "Description",
    n_pathways_to_show = 5
  )

  bar_layer <- ggplot2::ggplot_build(p)$data[[1]]

  expect_lte(
    nrow(bar_layer),
    10
  )

})
