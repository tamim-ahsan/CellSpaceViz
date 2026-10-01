test_that("example dataset is available and correctly formatted", {

  expect_true(exists("dea_volcano_df"))

  required_cols <- c(
    "cell_type",
    "gene",
    "log2FoldChange",
    "padj"
  )

  expect_true(
    all(required_cols %in% colnames(dea_volcano_df))
  )

})

test_that("returns a ggplot object", {

  p <- make_dea_volcano_plot(
    df = dea_volcano_df,
    lfc_col = "log2FoldChange",
    padj_col = "padj",
    gene_label_col = "gene"
  )

  expect_s3_class(p, "ggplot")

})

test_that("default colours work", {

  p <- make_dea_volcano_plot(
    df = dea_volcano_df,
    lfc_col = "log2FoldChange",
    padj_col = "padj",
    gene_label_col = "gene"
  )

  expect_s3_class(p, "ggplot")

})

test_that("custom colour vector works", {

  cols <- c(
    Upregulated = "#E41A1C",
    Downregulated = "#377EB8",
    Insignificant = "grey70"
  )

  p <- make_dea_volcano_plot(
    df = dea_volcano_df,
    lfc_col = "log2FoldChange",
    padj_col = "padj",
    gene_label_col = "gene",
    color_vec = cols
  )

  expect_s3_class(p, "ggplot")

})

test_that("plot title works", {

  p <- make_dea_volcano_plot(
    df = dea_volcano_df,
    lfc_col = "log2FoldChange",
    padj_col = "padj",
    gene_label_col = "gene",
    plot_title = "Volcano Plot"
  )

  expect_s3_class(p, "ggplot")

})

test_that("custom labels work", {

  p <- make_dea_volcano_plot(
    df = dea_volcano_df,
    lfc_col = "log2FoldChange",
    padj_col = "padj",
    gene_label_col = "gene",
    x_lab = "LogFC",
    y_lab = "Significance",
    color_lab = "DEA Category"
  )

  expect_s3_class(p, "ggplot")

})

test_that("missing log fold-change column throws error", {

  expect_error(

    make_dea_volcano_plot(
      df = dea_volcano_df,
      lfc_col = "bad_column",
      padj_col = "padj",
      gene_label_col = "gene"
    )

  )

})

test_that("missing adjusted p-value column throws error", {

  expect_error(

    make_dea_volcano_plot(
      df = dea_volcano_df,
      lfc_col = "log2FoldChange",
      padj_col = "bad_column",
      gene_label_col = "gene"
    )

  )

})

test_that("missing gene-label column throws error", {

  expect_error(

    make_dea_volcano_plot(
      df = dea_volcano_df,
      lfc_col = "log2FoldChange",
      padj_col = "padj",
      gene_label_col = "bad_column"
    )

  )

})

test_that("non-numeric log fold-change column throws error", {

  test_df <- dea_volcano_df

  test_df$log2FoldChange <- as.character(
    test_df$log2FoldChange
  )

  expect_error(

    make_dea_volcano_plot(
      df = test_df,
      lfc_col = "log2FoldChange",
      padj_col = "padj",
      gene_label_col = "gene"
    )

  )

})

test_that("non-numeric adjusted p-value column throws error", {

  test_df <- dea_volcano_df

  test_df$padj <- as.character(
    test_df$padj
  )

  expect_error(

    make_dea_volcano_plot(
      df = test_df,
      lfc_col = "log2FoldChange",
      padj_col = "padj",
      gene_label_col = "gene"
    )

  )

})

test_that("negative adjusted p-values throw error", {

  test_df <- dea_volcano_df

  test_df$padj[1] <- -0.01

  expect_error(

    make_dea_volcano_plot(
      df = test_df,
      lfc_col = "log2FoldChange",
      padj_col = "padj",
      gene_label_col = "gene"
    )

  )

})

test_that("zero adjusted p-values are handled", {

  test_df <- dea_volcano_df

  test_df$padj[1] <- 0

  p <- make_dea_volcano_plot(
    df = test_df,
    lfc_col = "log2FoldChange",
    padj_col = "padj",
    gene_label_col = "gene"
  )

  expect_s3_class(p, "ggplot")

})

test_that("invalid padj threshold throws error", {

  expect_error(
    make_dea_volcano_plot(
      df = dea_volcano_df,
      lfc_col = "log2FoldChange",
      padj_col = "padj",
      gene_label_col = "gene",
      padj_thresh = -1
    )
  )

  expect_error(
    make_dea_volcano_plot(
      df = dea_volcano_df,
      lfc_col = "log2FoldChange",
      padj_col = "padj",
      gene_label_col = "gene",
      padj_thresh = 1.5
    )
  )

})

test_that("invalid sig_gene_to_show throws error", {

  expect_error(

    make_dea_volcano_plot(
      df = dea_volcano_df,
      lfc_col = "log2FoldChange",
      padj_col = "padj",
      gene_label_col = "gene",
      sig_gene_to_show = -1
    )

  )

})

test_that("incorrect colour vector length throws error", {

  expect_error(

    make_dea_volcano_plot(
      df = dea_volcano_df,
      lfc_col = "log2FoldChange",
      padj_col = "padj",
      gene_label_col = "gene",
      color_vec = c(
        "red",
        "blue"
      )
    )

  )

})

test_that("named colours missing categories throw error", {

  cols <- c(
    Upregulated = "#E41A1C",
    Downregulated = "#377EB8"
  )

  expect_error(

    make_dea_volcano_plot(
      df = dea_volcano_df,
      lfc_col = "log2FoldChange",
      padj_col = "padj",
      gene_label_col = "gene",
      color_vec = cols
    )

  )

})

test_that("volcano point layer is generated", {

  p <- make_dea_volcano_plot(
    df = dea_volcano_df,
    lfc_col = "log2FoldChange",
    padj_col = "padj",
    gene_label_col = "gene"
  )

  point_layer <- ggplot2::ggplot_build(p)$data[[1]]

  expect_true(
    nrow(point_layer) > 0
  )

})

test_that("all genes are plotted", {

  p <- make_dea_volcano_plot(
    df = dea_volcano_df,
    lfc_col = "log2FoldChange",
    padj_col = "padj",
    gene_label_col = "gene"
  )

  point_layer <- ggplot2::ggplot_build(p)$data[[1]]

  expect_equal(
    nrow(point_layer),
    nrow(dea_volcano_df)
  )

})

test_that("gene labels are generated", {

  p <- make_dea_volcano_plot(
    df = dea_volcano_df,
    lfc_col = "log2FoldChange",
    padj_col = "padj",
    gene_label_col = "gene",
    sig_gene_to_show = 5
  )

  label_layer <- ggplot2::ggplot_build(p)$data[[3]]

  expect_true(
    nrow(label_layer) > 0
  )

})
