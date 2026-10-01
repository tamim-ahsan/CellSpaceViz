test_that("example dataset is available and correctly formatted", {

  expect_true(exists("dea_heatmap_df"))

  required_cols <- c(
    "cell_type",
    "gene",
    "log2FoldChange",
    "padj"
  )

  expect_true(
    all(required_cols %in% colnames(dea_heatmap_df))
  )

})

test_that("returns a ggplot object", {

  p <- make_dea_lfc_heatmap(
    df = dea_heatmap_df,
    lfc_col = "log2FoldChange",
    padj_col = "padj",
    gene_name_col = "gene",
    genes_to_show = c(
      "gene1",
      "gene2",
      "gene3",
      "gene4"
    ),
    condition_col = "cell_type"
  )

  expect_s3_class(p, "ggplot")

})

test_that("plot title works", {

  p <- make_dea_lfc_heatmap(
    df = dea_heatmap_df,
    lfc_col = "log2FoldChange",
    padj_col = "padj",
    gene_name_col = "gene",
    genes_to_show = c(
      "gene1",
      "gene2",
      "gene3",
      "gene4"
    ),
    condition_col = "cell_type",
    plot_title = "Example Heatmap"
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

    p <- make_dea_lfc_heatmap(
      df = dea_heatmap_df,
      lfc_col = "log2FoldChange",
      padj_col = "padj",
      gene_name_col = "gene",
      genes_to_show = c(
        "gene1",
        "gene2",
        "gene3",
        "gene4"
      ),
      condition_col = "cell_type",
      scale_fill_option = pal
    )

    expect_s3_class(p, "ggplot")

  }

})

test_that("custom significance levels work", {

  p <- make_dea_lfc_heatmap(
    df = dea_heatmap_df,
    lfc_col = "log2FoldChange",
    padj_col = "padj",
    gene_name_col = "gene",
    genes_to_show = c(
      "gene1",
      "gene2",
      "gene3",
      "gene4"
    ),
    condition_col = "cell_type",
    sig_levels = c(
      0.001,
      0.01,
      0.05
    )
  )

  expect_s3_class(p, "ggplot")

})

test_that("missing log fold-change column throws error", {

  expect_error(

    make_dea_lfc_heatmap(
      df = dea_heatmap_df,
      lfc_col = "bad_column",
      padj_col = "padj",
      gene_name_col = "gene",
      genes_to_show = "gene1",
      condition_col = "cell_type"
    )

  )

})

test_that("missing padj column throws error", {

  expect_error(

    make_dea_lfc_heatmap(
      df = dea_heatmap_df,
      lfc_col = "log2FoldChange",
      padj_col = "bad_column",
      gene_name_col = "gene",
      genes_to_show = "gene1",
      condition_col = "cell_type"
    )

  )

})

test_that("missing gene column throws error", {

  expect_error(

    make_dea_lfc_heatmap(
      df = dea_heatmap_df,
      lfc_col = "log2FoldChange",
      padj_col = "padj",
      gene_name_col = "bad_column",
      genes_to_show = "gene1",
      condition_col = "cell_type"
    )

  )

})

test_that("missing condition column throws error", {

  expect_error(

    make_dea_lfc_heatmap(
      df = dea_heatmap_df,
      lfc_col = "log2FoldChange",
      padj_col = "padj",
      gene_name_col = "gene",
      genes_to_show = "gene1",
      condition_col = "bad_column"
    )

  )

})

test_that("non numeric lfc column throws error", {

  test_df <- dea_heatmap_df

  test_df$log2FoldChange <- as.character(
    test_df$log2FoldChange
  )

  expect_error(

    make_dea_lfc_heatmap(
      df = test_df,
      lfc_col = "log2FoldChange",
      padj_col = "padj",
      gene_name_col = "gene",
      genes_to_show = "gene1",
      condition_col = "cell_type"
    )

  )

})

test_that("non numeric padj column throws error", {

  test_df <- dea_heatmap_df

  test_df$padj <- as.character(
    test_df$padj
  )

  expect_error(

    make_dea_lfc_heatmap(
      df = test_df,
      lfc_col = "log2FoldChange",
      padj_col = "padj",
      gene_name_col = "gene",
      genes_to_show = "gene1",
      condition_col = "cell_type"
    )

  )

})

test_that("empty genes_to_show throws error", {

  expect_error(

    make_dea_lfc_heatmap(
      df = dea_heatmap_df,
      lfc_col = "log2FoldChange",
      padj_col = "padj",
      gene_name_col = "gene",
      genes_to_show = character(0),
      condition_col = "cell_type"
    )

  )

})

test_that("unknown genes throw error", {

  expect_error(

    make_dea_lfc_heatmap(
      df = dea_heatmap_df,
      lfc_col = "log2FoldChange",
      padj_col = "padj",
      gene_name_col = "gene",
      genes_to_show = c(
        "fake_gene_1",
        "fake_gene_2"
      ),
      condition_col = "cell_type"
    )

  )

})

test_that("invalid viridis palette throws error", {

  expect_error(

    make_dea_lfc_heatmap(
      df = dea_heatmap_df,
      lfc_col = "log2FoldChange",
      padj_col = "padj",
      gene_name_col = "gene",
      genes_to_show = "gene1",
      condition_col = "cell_type",
      scale_fill_option = "rainbow"
    )

  )

})

test_that("invalid significance levels throw error", {

  expect_error(

    make_dea_lfc_heatmap(
      df = dea_heatmap_df,
      lfc_col = "log2FoldChange",
      padj_col = "padj",
      gene_name_col = "gene",
      genes_to_show = "gene1",
      condition_col = "cell_type",
      sig_levels = character(0)
    )

  )

})

test_that("heatmap tile layer is generated", {

  p <- make_dea_lfc_heatmap(
    df = dea_heatmap_df,
    lfc_col = "log2FoldChange",
    padj_col = "padj",
    gene_name_col = "gene",
    genes_to_show = c(
      "gene1",
      "gene2",
      "gene3",
      "gene4"
    ),
    condition_col = "cell_type"
  )

  tile_layer <- ggplot2::ggplot_build(p)$data[[1]]

  expect_equal(
    nrow(tile_layer),
    16
  )

})

test_that("significance labels are generated", {

  p <- make_dea_lfc_heatmap(
    df = dea_heatmap_df,
    lfc_col = "log2FoldChange",
    padj_col = "padj",
    gene_name_col = "gene",
    genes_to_show = c(
      "gene1",
      "gene2",
      "gene3",
      "gene4"
    ),
    condition_col = "cell_type"
  )

  text_layer <- ggplot2::ggplot_build(p)$data[[2]]

  expect_equal(
    nrow(text_layer),
    16
  )

})

test_that("all selected genes are plotted", {

  p <- make_dea_lfc_heatmap(
    df = dea_heatmap_df,
    lfc_col = "log2FoldChange",
    padj_col = "padj",
    gene_name_col = "gene",
    genes_to_show = c(
      "gene1",
      "gene2",
      "gene3",
      "gene4"
    ),
    condition_col = "cell_type"
  )

  expect_s3_class(
    p,
    "ggplot"
  )

})
