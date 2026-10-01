test_that("example dataset is available and correctly formatted", {

  expect_true(exists("marker_df"))

  required_cols <- c(
    "cell_type",
    "geneA1", "geneA2", "geneA3",
    "geneB1", "geneB2", "geneB3",
    "geneC1", "geneC2", "geneC3",
    "geneD1", "geneD2", "geneD3"
  )

  expect_true(
    all(required_cols %in% colnames(marker_df))
  )

})

marker_genes <- c(
  "geneA1", "geneA2", "geneA3",
  "geneB1", "geneB2", "geneB3",
  "geneC1", "geneC2", "geneC3",
  "geneD1", "geneD2", "geneD3"
)

gene_groups <- c(
  rep("Cell-type A", 3),
  rep("Cell-type B", 3),
  rep("Cell-type C", 3),
  rep("Cell-type D", 3)
)

test_that("returns a ggplot object", {

  p <- make_marker_dotplot(
    df = marker_df,
    marker_genes = marker_genes,
    ct_col = "cell_type",
    ct_groups = gene_groups
  )

  expect_s3_class(p, "ggplot")

})

test_that("works with all supported viridis palettes", {

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

    p <- make_marker_dotplot(
      df = marker_df,
      marker_genes = marker_genes,
      ct_col = "cell_type",
      ct_groups = gene_groups,
      color_scale_option = pal
    )

    expect_s3_class(p, "ggplot")
  }

})

test_that("custom plot title works", {

  p <- make_marker_dotplot(
    df = marker_df,
    marker_genes = marker_genes,
    ct_col = "cell_type",
    ct_groups = gene_groups,
    plot_title = "Test Dot Plot"
  )

  expect_s3_class(p, "ggplot")

})

test_that("custom colour scale title works", {

  p <- make_marker_dotplot(
    df = marker_df,
    marker_genes = marker_genes,
    ct_col = "cell_type",
    ct_groups = gene_groups,
    color_scale_name = "Average Expression"
  )

  expect_s3_class(p, "ggplot")

})

test_that("cells_or_spots accepts cells", {

  p <- make_marker_dotplot(
    df = marker_df,
    marker_genes = marker_genes,
    ct_col = "cell_type",
    ct_groups = gene_groups,
    cells_or_spots = "cells"
  )

  expect_s3_class(p, "ggplot")

})

test_that("cells_or_spots accepts spots", {

  p <- make_marker_dotplot(
    df = marker_df,
    marker_genes = marker_genes,
    ct_col = "cell_type",
    ct_groups = gene_groups,
    cells_or_spots = "spots"
  )

  expect_s3_class(p, "ggplot")

})

test_that("invalid cells_or_spots throws error", {

  expect_error(
    make_marker_dotplot(
      df = marker_df,
      marker_genes = marker_genes,
      ct_col = "cell_type",
      ct_groups = gene_groups,
      cells_or_spots = "bananas"
    )
  )

})

test_that("missing cell type column throws error", {

  expect_error(
    make_marker_dotplot(
      df = marker_df,
      marker_genes = marker_genes,
      ct_col = "missing_column",
      ct_groups = gene_groups
    )
  )

})

test_that("missing marker gene column throws error", {

  expect_error(
    make_marker_dotplot(
      df = marker_df,
      marker_genes = c(marker_genes, "fake_gene"),
      ct_col = "cell_type",
      ct_groups = c(gene_groups, "Fake")
    )
  )

})

test_that("marker gene and group lengths must match", {

  expect_error(
    make_marker_dotplot(
      df = marker_df,
      marker_genes = marker_genes,
      ct_col = "cell_type",
      ct_groups = gene_groups[-1]
    )
  )

})

test_that("duplicate marker genes throw error", {

  expect_error(
    make_marker_dotplot(
      df = marker_df,
      marker_genes = c(marker_genes, "geneA1"),
      ct_col = "cell_type",
      ct_groups = c(gene_groups, "Cell-type A")
    )
  )

})

test_that("missing gene group entries throw error", {

  bad_groups <- gene_groups
  bad_groups[1] <- NA

  expect_error(
    make_marker_dotplot(
      df = marker_df,
      marker_genes = marker_genes,
      ct_col = "cell_type",
      ct_groups = bad_groups
    )
  )

})

test_that("invalid viridis option throws error", {

  expect_error(
    make_marker_dotplot(
      df = marker_df,
      marker_genes = marker_genes,
      ct_col = "cell_type",
      ct_groups = gene_groups,
      color_scale_option = "rainbow"
    )
  )

})

test_that("factor levels are preserved", {

  test_df <- marker_df

  test_df$cell_type <- factor(
    test_df$cell_type,
    levels = c(
      "Cell-type D",
      "Cell-type C",
      "Cell-type B",
      "Cell-type A"
    )
  )

  p <- make_marker_dotplot(
    df = test_df,
    marker_genes = marker_genes,
    ct_col = "cell_type",
    ct_groups = gene_groups
  )

  expect_s3_class(p, "ggplot")

})

test_that("non-numeric marker columns throw error", {

  test_df <- marker_df

  test_df$geneA1 <- as.character(test_df$geneA1)

  expect_error(
    make_marker_dotplot(
      df = test_df,
      marker_genes = marker_genes,
      ct_col = "cell_type",
      ct_groups = gene_groups
    )
  )

})

test_that("plot contains facet structure", {

  p <- make_marker_dotplot(
    df = marker_df,
    marker_genes = marker_genes,
    ct_col = "cell_type",
    ct_groups = gene_groups
  )

  expect_false(is.null(p$facet))

})
