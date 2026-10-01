test_that("example dataset is available and correctly formatted", {

  expect_true(exists("deg_results_df"))

  required_cols <- c(
    "cell_type",
    "gene",
    "log2FoldChange",
    "padj"
  )

  expect_true(
    all(required_cols %in% colnames(deg_results_df))
  )

})

test_that("returns a ggplot object", {

  p <- make_deg_metric_plot(
    df = deg_results_df,
    ct_col = "cell_type",
    lfc_col = "log2FoldChange",
    padj_col = "padj"
  )

  expect_s3_class(p, "ggplot")

})

test_that("default colours work", {

  p <- make_deg_metric_plot(
    df = deg_results_df,
    ct_col = "cell_type",
    lfc_col = "log2FoldChange",
    padj_col = "padj"
  )

  expect_s3_class(p, "ggplot")

})

test_that("custom colour vector works", {

  cols <- c(
    "Cell-type A" = "#E41A1C",
    "Cell-type B" = "#377EB8",
    "Cell-type C" = "#4DAF4A",
    "Cell-type D" = "#984EA3"
  )

  p <- make_deg_metric_plot(
    df = deg_results_df,
    ct_col = "cell_type",
    lfc_col = "log2FoldChange",
    padj_col = "padj",
    color_vec = cols
  )

  expect_s3_class(p, "ggplot")

})

test_that("custom legend labels work", {

  p <- make_deg_metric_plot(
    df = deg_results_df,
    ct_col = "cell_type",
    lfc_col = "log2FoldChange",
    padj_col = "padj",
    color_lab = "Cell Type",
    size_lab = "Significance"
  )

  expect_s3_class(p, "ggplot")

})

test_that("plot title works", {

  p <- make_deg_metric_plot(
    df = deg_results_df,
    ct_col = "cell_type",
    lfc_col = "log2FoldChange",
    padj_col = "padj",
    plot_title = "Example DEG Metric Plot"
  )

  expect_s3_class(p, "ggplot")

})

test_that("custom axis labels work", {

  p <- make_deg_metric_plot(
    df = deg_results_df,
    ct_col = "cell_type",
    lfc_col = "log2FoldChange",
    padj_col = "padj",
    x_lab = "LogFC",
    y_lab = "Cell Type"
  )

  expect_s3_class(p, "ggplot")

})

test_that("missing cell type column throws error", {

  expect_error(

    make_deg_metric_plot(
      df = deg_results_df,
      ct_col = "bad_column",
      lfc_col = "log2FoldChange",
      padj_col = "padj"
    )

  )

})

test_that("missing log fold-change column throws error", {

  expect_error(

    make_deg_metric_plot(
      df = deg_results_df,
      ct_col = "cell_type",
      lfc_col = "bad_column",
      padj_col = "padj"
    )

  )

})

test_that("missing adjusted p-value column throws error", {

  expect_error(

    make_deg_metric_plot(
      df = deg_results_df,
      ct_col = "cell_type",
      lfc_col = "log2FoldChange",
      padj_col = "bad_column"
    )

  )

})

test_that("non-numeric log fold-change column throws error", {

  test_df <- deg_results_df

  test_df$log2FoldChange <- as.character(
    test_df$log2FoldChange
  )

  expect_error(

    make_deg_metric_plot(
      df = test_df,
      ct_col = "cell_type",
      lfc_col = "log2FoldChange",
      padj_col = "padj"
    )

  )

})

test_that("non-numeric adjusted p-value column throws error", {

  test_df <- deg_results_df

  test_df$padj <- as.character(
    test_df$padj
  )

  expect_error(

    make_deg_metric_plot(
      df = test_df,
      ct_col = "cell_type",
      lfc_col = "log2FoldChange",
      padj_col = "padj"
    )

  )

})

test_that("negative adjusted p-values throw error", {

  test_df <- deg_results_df

  test_df$padj[1] <- -0.01

  expect_error(

    make_deg_metric_plot(
      df = test_df,
      ct_col = "cell_type",
      lfc_col = "log2FoldChange",
      padj_col = "padj"
    )

  )

})

test_that("insufficient colours throw error", {

  expect_error(

    make_deg_metric_plot(
      df = deg_results_df,
      ct_col = "cell_type",
      lfc_col = "log2FoldChange",
      padj_col = "padj",
      color_vec = c(
        "red",
        "blue"
      )
    )

  )

})

test_that("named colours missing cell types throw error", {

  cols <- c(
    "Cell-type A" = "#E41A1C",
    "Cell-type B" = "#377EB8"
  )

  expect_error(

    make_deg_metric_plot(
      df = deg_results_df,
      ct_col = "cell_type",
      lfc_col = "log2FoldChange",
      padj_col = "padj",
      color_vec = cols
    )

  )

})

test_that("factor level ordering is preserved", {

  test_df <- deg_results_df

  test_df$cell_type <- factor(
    test_df$cell_type,
    levels = c(
      "Cell-type D",
      "Cell-type C",
      "Cell-type B",
      "Cell-type A"
    )
  )

  p <- make_deg_metric_plot(
    df = test_df,
    ct_col = "cell_type",
    lfc_col = "log2FoldChange",
    padj_col = "padj"
  )

  expect_s3_class(p, "ggplot")

})

test_that("zero adjusted p-values are handled", {

  test_df <- deg_results_df

  test_df$padj[1] <- 0

  p <- make_deg_metric_plot(
    df = test_df,
    ct_col = "cell_type",
    lfc_col = "log2FoldChange",
    padj_col = "padj"
  )

  expect_s3_class(p, "ggplot")

})

test_that("bubble layer is generated", {

  p <- make_deg_metric_plot(
    df = deg_results_df,
    ct_col = "cell_type",
    lfc_col = "log2FoldChange",
    padj_col = "padj"
  )

  bubble_layer <- ggplot2::ggplot_build(p)$data[[1]]

  expect_true(
    nrow(bubble_layer) > 0
  )

})

test_that("all observations are plotted", {

  p <- make_deg_metric_plot(
    df = deg_results_df,
    ct_col = "cell_type",
    lfc_col = "log2FoldChange",
    padj_col = "padj"
  )

  bubble_layer <- ggplot2::ggplot_build(p)$data[[1]]

  expect_equal(
    nrow(bubble_layer),
    nrow(deg_results_df)
  )

})
