test_that("example dataset is available and correctly formatted", {

  expect_true(exists("reduced_dim_df"))

  expect_true(
    all(
      c("UMAP1", "UMAP2", "cell_type", "modality") %in%
        colnames(reduced_dim_df)
    )
  )

})

test_that("returns a ggplot object", {

  p <- make_reduced_dim_plots(
    df = reduced_dim_df,
    ct_col = "cell_type",
    dataset_col = "modality"
  )

  expect_s3_class(p, "ggplot")

})

test_that("default colours are generated successfully", {

  p <- make_reduced_dim_plots(
    df = reduced_dim_df,
    ct_col = "cell_type",
    dataset_col = "modality"
  )

  expect_s3_class(p, "ggplot")

})

test_that("custom named colours are accepted", {

  cols <- c(
    "B cell"   = "#E41A1C",
    "Monocyte" = "#377EB8",
    "NK cell"  = "#4DAF4A",
    "T cell"   = "#984EA3"
  )

  p <- make_reduced_dim_plots(
    df = reduced_dim_df,
    ct_col = "cell_type",
    dataset_col = "modality",
    color_vec = cols
  )

  expect_s3_class(p, "ggplot")

})

test_that("missing cell-type column throws error", {

  expect_error(
    make_reduced_dim_plots(
      df = reduced_dim_df,
      ct_col = "missing_column",
      dataset_col = "modality"
    )
  )

})

test_that("missing dataset column throws error", {

  expect_error(
    make_reduced_dim_plots(
      df = reduced_dim_df,
      ct_col = "cell_type",
      dataset_col = "missing_column"
    )
  )

})

test_that("missing dim_1 column throws error", {

  expect_error(
    make_reduced_dim_plots(
      df = reduced_dim_df,
      dim_1 = "missing_dim",
      ct_col = "cell_type",
      dataset_col = "modality"
    )
  )

})

test_that("missing dim_2 column throws error", {

  expect_error(
    make_reduced_dim_plots(
      df = reduced_dim_df,
      dim_2 = "missing_dim",
      ct_col = "cell_type",
      dataset_col = "modality"
    )
  )

})

test_that("insufficient unnamed colour vector throws error", {

  expect_error(
    make_reduced_dim_plots(
      df = reduced_dim_df,
      ct_col = "cell_type",
      dataset_col = "modality",
      color_vec = c("red", "blue")
    )
  )

})

test_that("incomplete named colour vector throws error", {

  cols <- c(
    "B cell" = "red",
    "T cell" = "blue"
  )

  expect_error(
    make_reduced_dim_plots(
      df = reduced_dim_df,
      ct_col = "cell_type",
      dataset_col = "modality",
      color_vec = cols
    )
  )

})

test_that("factor level ordering is preserved", {

  test_df <- reduced_dim_df

  new_levels <- c(
    "NK cell",
    "T cell",
    "B cell",
    "Monocyte"
  )

  test_df$cell_type <- factor(
    test_df$cell_type,
    levels = new_levels
  )

  p <- make_reduced_dim_plots(
    df = test_df,
    ct_col = "cell_type",
    dataset_col = "modality"
  )

  expect_s3_class(p, "ggplot")

})

test_that("custom dimension names work", {

  test_df <- reduced_dim_df

  test_df$PC1 <- test_df$UMAP1
  test_df$PC2 <- test_df$UMAP2

  p <- make_reduced_dim_plots(
    df = test_df,
    dim_1 = "PC1",
    dim_2 = "PC2",
    ct_col = "cell_type",
    dataset_col = "modality"
  )

  expect_s3_class(p, "ggplot")

})

test_that("custom axis labels work", {

  p <- make_reduced_dim_plots(
    df = reduced_dim_df,
    ct_col = "cell_type",
    dataset_col = "modality",
    axis_label_base = "tSNE"
  )

  expect_s3_class(p, "ggplot")

})

test_that("custom point size and alpha work", {

  p <- make_reduced_dim_plots(
    df = reduced_dim_df,
    ct_col = "cell_type",
    dataset_col = "modality",
    point_size = 2,
    point_alpha = 0.5
  )

  expect_s3_class(p, "ggplot")

})

test_that("works with tibble input", {

  skip_if_not_installed("tibble")

  test_df <- tibble::as_tibble(reduced_dim_df)

  p <- make_reduced_dim_plots(
    df = test_df,
    ct_col = "cell_type",
    dataset_col = "modality"
  )

  expect_s3_class(p, "ggplot")

})

test_that("returns object with facetting", {

  p <- make_reduced_dim_plots(
    df = reduced_dim_df,
    ct_col = "cell_type",
    dataset_col = "modality"
  )

  expect_false(is.null(p$facet))

})
