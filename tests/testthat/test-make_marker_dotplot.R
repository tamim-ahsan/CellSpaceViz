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

test_that("marker dataset is available", {

  expect_true(exists("marker_df"))

})

test_that("returns ggplot object", {

  p <- make_marker_dotplot(
    df = marker_df,
    marker_genes = marker_genes,
    ct_col = "cell_type",
    gene_groups = gene_groups
  )

  expect_s3_class(p, "ggplot")

})

test_that("all viridis palettes work", {

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
      gene_groups = gene_groups,
      color_scale_option = pal
    )

    expect_s3_class(p, "ggplot")

  }

})

test_that("invalid palette throws error", {

  expect_error(

    make_marker_dotplot(
      df = marker_df,
      marker_genes = marker_genes,
      ct_col = "cell_type",
      gene_groups = gene_groups,
      color_scale_option = "rainbow"
    )

  )

})

test_that("duplicate markers throw error", {

  expect_error(

    make_marker_dotplot(
      df = marker_df,
      marker_genes = c(
        marker_genes,
        "geneA1"
      ),
      ct_col = "cell_type",
      gene_groups = c(
        gene_groups,
        "Cell-type A"
      )
    )

  )

})

test_that("group length mismatch throws error", {

  expect_error(

    make_marker_dotplot(
      df = marker_df,
      marker_genes = marker_genes,
      ct_col = "cell_type",
      gene_groups = gene_groups[-1]
    )

  )

})

test_that("missing marker throws error", {

  expect_error(

    make_marker_dotplot(
      df = marker_df,
      marker_genes = c(
        marker_genes,
        "FakeGene"
      ),
      ct_col = "cell_type",
      gene_groups = c(
        gene_groups,
        "FakeGroup"
      )
    )

  )

})

test_that("non numeric marker columns throw error", {

  tmp <- marker_df

  tmp$geneA1 <- as.character(
    tmp$geneA1
  )

  expect_error(

    make_marker_dotplot(
      df = tmp,
      marker_genes = marker_genes,
      ct_col = "cell_type",
      gene_groups = gene_groups
    )

  )

})
