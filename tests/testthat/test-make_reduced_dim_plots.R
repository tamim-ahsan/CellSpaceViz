test_that("example dataset is available", {

  expect_true(exists("reduced_dim_df"))

})

test_that("required columns exist", {

  expect_true(
    all(
      c(
        "UMAP1",
        "UMAP2",
        "cell_type",
        "modality"
      ) %in% colnames(reduced_dim_df)
    )
  )

})

test_that("returns ggplot object", {

  p <- make_reduced_dim_plots(
    df = reduced_dim_df,
    ct_col = "cell_type",
    dataset_col = "modality"
  )

  expect_s3_class(p, "ggplot")

})

test_that("custom color vector works", {

  levels <- unique(
    as.character(
      reduced_dim_df$cell_type
    )
  )

  cols <- setNames(
    scales::hue_pal()(length(levels)),
    levels
  )

  p <- make_reduced_dim_plots(
    df = reduced_dim_df,
    ct_col = "cell_type",
    dataset_col = "modality",
    color_vec = cols
  )

  expect_s3_class(p, "ggplot")

})

test_that("custom dimensions work", {

  tmp <- reduced_dim_df

  tmp$PC1 <- tmp$UMAP1
  tmp$PC2 <- tmp$UMAP2

  p <- make_reduced_dim_plots(
    df = tmp,
    dim_1 = "PC1",
    dim_2 = "PC2",
    ct_col = "cell_type",
    dataset_col = "modality"
  )

  expect_s3_class(p, "ggplot")

})

test_that("missing celltype column throws error", {

  expect_error(

    make_reduced_dim_plots(
      df = reduced_dim_df,
      ct_col = "bad_column",
      dataset_col = "modality"
    )

  )

})

test_that("missing dataset column throws error", {

  expect_error(

    make_reduced_dim_plots(
      df = reduced_dim_df,
      ct_col = "cell_type",
      dataset_col = "bad_column"
    )

  )

})

test_that("insufficient colors throw error", {

  expect_error(

    make_reduced_dim_plots(
      df = reduced_dim_df,
      ct_col = "cell_type",
      dataset_col = "modality",
      color_vec = c("red", "blue")
    )

  )

})
