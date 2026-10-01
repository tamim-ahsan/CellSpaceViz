test_that("example dataset is available and correctly formatted", {

  expect_true(exists("deg_summary_df"))

  required_cols <- c(
    "cell_type",
    "up_regulated",
    "down_regulated",
    "total_degs"
  )

  expect_true(
    all(required_cols %in% colnames(deg_summary_df))
  )

})

test_that("returns a ggplot object", {

  p <- make_deg_count_plot(
    df = deg_summary_df,
    ct_col = "cell_type",
    deg_num_col = "total_degs",
    up_col = "up_regulated",
    down_col = "down_regulated"
  )

  expect_s3_class(p, "ggplot")

})

test_that("default colours work", {

  p <- make_deg_count_plot(
    df = deg_summary_df,
    ct_col = "cell_type",
    deg_num_col = "total_degs",
    up_col = "up_regulated",
    down_col = "down_regulated"
  )

  expect_s3_class(p, "ggplot")

})

test_that("custom colours work", {

  cols <- c(
    up_regulated = "#E41A1C",
    down_regulated = "#377EB8"
  )

  p <- make_deg_count_plot(
    df = deg_summary_df,
    ct_col = "cell_type",
    deg_num_col = "total_degs",
    up_col = "up_regulated",
    down_col = "down_regulated",
    color_vec = cols
  )

  expect_s3_class(p, "ggplot")

})

test_that("plot title works", {

  p <- make_deg_count_plot(
    df = deg_summary_df,
    ct_col = "cell_type",
    deg_num_col = "total_degs",
    up_col = "up_regulated",
    down_col = "down_regulated",
    plot_title = "Example DEG Plot"
  )

  expect_s3_class(p, "ggplot")

})

test_that("custom axis labels work", {

  p <- make_deg_count_plot(
    df = deg_summary_df,
    ct_col = "cell_type",
    deg_num_col = "total_degs",
    up_col = "up_regulated",
    down_col = "down_regulated",
    x_lab = "DEG Count",
    y_lab = "Cell Type"
  )

  expect_s3_class(p, "ggplot")

})

test_that("missing cell-type column throws error", {

  expect_error(

    make_deg_count_plot(
      df = deg_summary_df,
      ct_col = "missing_column",
      deg_num_col = "total_degs",
      up_col = "up_regulated",
      down_col = "down_regulated"
    )

  )

})

test_that("missing DEG count column throws error", {

  expect_error(

    make_deg_count_plot(
      df = deg_summary_df,
      ct_col = "cell_type",
      deg_num_col = "missing_column",
      up_col = "up_regulated",
      down_col = "down_regulated"
    )

  )

})

test_that("missing up-regulated column throws error", {

  expect_error(

    make_deg_count_plot(
      df = deg_summary_df,
      ct_col = "cell_type",
      deg_num_col = "total_degs",
      up_col = "missing_column",
      down_col = "down_regulated"
    )

  )

})

test_that("missing down-regulated column throws error", {

  expect_error(

    make_deg_count_plot(
      df = deg_summary_df,
      ct_col = "cell_type",
      deg_num_col = "total_degs",
      up_col = "up_regulated",
      down_col = "missing_column"
    )

  )

})

test_that("incorrect colour vector length throws error", {

  expect_error(

    make_deg_count_plot(
      df = deg_summary_df,
      ct_col = "cell_type",
      deg_num_col = "total_degs",
      up_col = "up_regulated",
      down_col = "down_regulated",
      color_vec = c("red")
    )

  )

})

test_that("named colours missing regulation category throws error", {

  cols <- c(
    up_regulated = "#E41A1C"
  )

  expect_error(

    make_deg_count_plot(
      df = deg_summary_df,
      ct_col = "cell_type",
      deg_num_col = "total_degs",
      up_col = "up_regulated",
      down_col = "down_regulated",
      color_vec = cols
    )

  )

})

test_that("factor ordering is preserved", {

  test_df <- deg_summary_df

  test_df$cell_type <- factor(
    test_df$cell_type,
    levels = c(
      "Cell-type D",
      "Cell-type C",
      "Cell-type B",
      "Cell-type A"
    )
  )

  p <- make_deg_count_plot(
    df = test_df,
    ct_col = "cell_type",
    deg_num_col = "total_degs",
    up_col = "up_regulated",
    down_col = "down_regulated"
  )

  expect_s3_class(p, "ggplot")

})

test_that("stacked bar layer is generated", {

  p <- make_deg_count_plot(
    df = deg_summary_df,
    ct_col = "cell_type",
    deg_num_col = "total_degs",
    up_col = "up_regulated",
    down_col = "down_regulated"
  )

  plot_data <- ggplot2::ggplot_build(p)$data[[1]]

  expect_true(
    nrow(plot_data) > 0
  )

})

test_that("segment labels are generated", {

  p <- make_deg_count_plot(
    df = deg_summary_df,
    ct_col = "cell_type",
    deg_num_col = "total_degs",
    up_col = "up_regulated",
    down_col = "down_regulated"
  )

  label_layer <- ggplot2::ggplot_build(p)$data[[2]]

  expect_true(
    nrow(label_layer) > 0
  )

})

test_that("total DEG labels are generated", {

  p <- make_deg_count_plot(
    df = deg_summary_df,
    ct_col = "cell_type",
    deg_num_col = "total_degs",
    up_col = "up_regulated",
    down_col = "down_regulated"
  )

  total_layer <- ggplot2::ggplot_build(p)$data[[3]]

  expect_equal(
    nrow(total_layer),
    nrow(deg_summary_df)
  )

})
