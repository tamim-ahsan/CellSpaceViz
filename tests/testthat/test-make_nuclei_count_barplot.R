test_that("example dataset is available and correctly formatted", {

  expect_true(exists("nuclei_count_df"))

  expect_true(
    "cell_type" %in% colnames(nuclei_count_df)
  )

})

test_that("returns a ggplot object", {

  p <- make_nuclei_count_barplot(
    df = nuclei_count_df,
    ct_col = "cell_type"
  )

  expect_s3_class(p, "ggplot")

})

test_that("accepts custom colour vector", {

  cols <- c(
    "Cell-type A" = "#E41A1C",
    "Cell-type B" = "#377EB8",
    "Cell-type C" = "#4DAF4A",
    "Cell-type D" = "#984EA3"
  )

  p <- make_nuclei_count_barplot(
    df = nuclei_count_df,
    ct_col = "cell_type",
    color_vec = cols
  )

  expect_s3_class(p, "ggplot")

})

test_that("accepts Cell, Nucleus and Spot labels", {

  p1 <- make_nuclei_count_barplot(
    df = nuclei_count_df,
    ct_col = "cell_type",
    cell_or_spot = "Cell"
  )

  p2 <- make_nuclei_count_barplot(
    df = nuclei_count_df,
    ct_col = "cell_type",
    cell_or_spot = "Nucleus"
  )

  p3 <- make_nuclei_count_barplot(
    df = nuclei_count_df,
    ct_col = "cell_type",
    cell_or_spot = "Spot"
  )

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")

})

test_that("invalid cell_or_spot throws error", {

  expect_error(

    make_nuclei_count_barplot(
      df = nuclei_count_df,
      ct_col = "cell_type",
      cell_or_spot = "Banana"
    )

  )

})

test_that("missing cell type column throws error", {

  expect_error(

    make_nuclei_count_barplot(
      df = nuclei_count_df,
      ct_col = "missing_column"
    )

  )

})

test_that("insufficient colours throw error", {

  expect_error(

    make_nuclei_count_barplot(
      df = nuclei_count_df,
      ct_col = "cell_type",
      color_vec = c(
        "#FF0000",
        "#0000FF"
      )
    )

  )

})

test_that("named colour vector missing groups throws error", {

  cols <- c(
    "Cell-type A" = "#E41A1C",
    "Cell-type B" = "#377EB8"
  )

  expect_error(

    make_nuclei_count_barplot(
      df = nuclei_count_df,
      ct_col = "cell_type",
      color_vec = cols
    )

  )

})

test_that("plot title works", {

  p <- make_nuclei_count_barplot(
    df = nuclei_count_df,
    ct_col = "cell_type",
    plot_title = "Example Plot"
  )

  expect_s3_class(p, "ggplot")

})

test_that("factor level ordering is preserved", {

  df <- nuclei_count_df

  df$cell_type <- factor(
    df$cell_type,
    levels = c(
      "Cell-type D",
      "Cell-type C",
      "Cell-type B",
      "Cell-type A"
    )
  )

  p <- make_nuclei_count_barplot(
    df = df,
    ct_col = "cell_type"
  )

  expect_s3_class(p, "ggplot")

})

test_that("bar counts are correct", {

  p <- make_nuclei_count_barplot(
    df = nuclei_count_df,
    ct_col = "cell_type"
  )

  plot_data <- ggplot2::ggplot_build(p)$data[[1]]

  expect_equal(
    nrow(plot_data),
    4
  )

})

test_that("labels are generated", {

  p <- make_nuclei_count_barplot(
    df = nuclei_count_df,
    ct_col = "cell_type"
  )

  text_layer <- ggplot2::ggplot_build(p)$data[[2]]

  expect_equal(
    nrow(text_layer),
    4
  )

})
