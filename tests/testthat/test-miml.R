pheno <- c("diarrhea_grade", "sampling_point", "slide_print_batch",
           "trg_score", "ypn_stage", "ypt_stage")

test_that("both panels have the documented structure", {
  data(RFS_data, package = "MIML")
  data(MFS_data, package = "MIML")

  for (d in list(RFS_data, MFS_data)) {
    expect_equal(dim(d), c(288L, 516L))
    expect_equal(length(unique(d$Patient_ID)), 80L)
    expect_true(all(pheno %in% names(d)))
  }
  expect_true("RFS" %in% names(RFS_data))
  expect_true("MFS" %in% names(MFS_data))
  expect_identical(RFS_data$Patient_ID, MFS_data$Patient_ID)

  mk <- setdiff(names(RFS_data), c("Patient_ID", "RFS", "event", pheno))
  expect_equal(length(mk), 507L)
  expect_true(all(vapply(RFS_data[, mk], is.numeric, logical(1))))
})

test_that("patient-level event counts are as documented", {
  data(RFS_data, package = "MIML")
  data(MFS_data, package = "MIML")
  expect_equal(sum(RFS_data$event[!duplicated(RFS_data$Patient_ID)]), 5L)
  expect_equal(sum(MFS_data$event[!duplicated(MFS_data$Patient_ID)]), 26L)
})

test_that("miml_rfs fills every missing cell and leaves observed values alone", {
  skip_if_not_installed("lightgbm")
  data(RFS_data, package = "MIML")
  mk <- names(RFS_data)[4:6]
  d  <- RFS_data[1:60, c("Patient_ID", "RFS", "event", mk)]

  res <- miml_rfs(d, m = 1, maxit = 1, nrounds = 3,
                  n_clusters = 0, verbose = FALSE)

  expect_s3_class(res, "miml_analysis")
  expect_equal(res$endpoint, "RFS")
  expect_equal(sum(is.na(res$completed[, mk])), 0L)
  obs <- !is.na(d[[mk[1]]])
  expect_equal(res$completed[[mk[1]]][obs], d[[mk[1]]][obs])
  expect_output(print(res), "MIML analysis -- RFS")
})

test_that("the same seed gives the same completion twice", {
  skip_if_not_installed("lightgbm")
  data(RFS_data, package = "MIML")
  mk <- names(RFS_data)[4:6]
  d  <- RFS_data[1:60, c("Patient_ID", "RFS", "event", mk)]

  a <- miml_rfs(d, m = 1, maxit = 1, nrounds = 5,
                n_clusters = 0, verbose = FALSE)
  b <- miml_rfs(d, m = 1, maxit = 1, nrounds = 5,
                n_clusters = 0, verbose = FALSE)
  expect_equal(a$completed[, mk], b$completed[, mk])
})

test_that("miml_mfs reads the metastasis endpoint", {
  skip_if_not_installed("lightgbm")
  data(MFS_data, package = "MIML")
  mk <- names(MFS_data)[4:6]
  d  <- MFS_data[1:60, c("Patient_ID", "MFS", "event", mk)]

  res <- miml_mfs(d, m = 1, maxit = 1, nrounds = 3,
                  n_clusters = 0, verbose = FALSE)
  expect_equal(res$endpoint, "MFS")
  expect_equal(res$settings$time_var, "MFS")
})

test_that("the caller's RNG stream is left untouched", {
  skip_if_not_installed("lightgbm")
  data(RFS_data, package = "MIML")
  d <- RFS_data[1:60, c("Patient_ID", "RFS", "event", names(RFS_data)[4:6])]

  set.seed(99); before <- .Random.seed
  invisible(miml_rfs(d, m = 1, maxit = 1, nrounds = 3,
                     n_clusters = 0, verbose = FALSE))
  expect_identical(.Random.seed, before)
})

test_that("predictors = 'all' keeps the outcome out of the pool", {
  skip_if_not_installed("lightgbm")
  data(RFS_data, package = "MIML")
  d <- RFS_data[1:60, c("Patient_ID", "RFS", "event", names(RFS_data)[4:6])]

  res <- miml_rfs(d, m = 1, maxit = 1, nrounds = 3, n_clusters = 0,
                  predictors = "all", verbose = FALSE)
  expect_equal(res$settings$predictors, "all")
  expect_equal(sum(is.na(res$completed[, names(RFS_data)[4:6]])), 0L)
})

test_that("a missing endpoint column is caught before any work", {
  data(RFS_data, package = "MIML")
  d <- RFS_data[1:20, c("Patient_ID", "event", names(RFS_data)[4:6])]
  expect_error(miml_rfs(d, verbose = FALSE), "Columns not found")
})
