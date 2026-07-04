test_that("KaguModel$new() assigns GPMechanism by default", {
  model <- KaguModel$new(dag = CHAIN_DAG)
  for (node in names(CHAIN_DAG)) {
    expect_true(inherits(model$mechanisms[[node]], "GPMechanism"))
  }
})

test_that("KaguModel$new() errors on invalid DAG", {
  expect_error(KaguModel$new(dag = list(x = c("y"))), "unknown parent")
})

test_that("print shows 'unfitted' before fit()", {
  model <- KaguModel$new(dag = CHAIN_DAG)
  out   <- capture.output(print(model))
  expect_match(paste(out, collapse = "\n"), "unfitted")
})

test_that("print shows 'fitted' after fit()", {
  skip_on_cran()
  data  <- make_chain_data()
  model <- KaguModel$new(dag = CHAIN_DAG)
  do.call(model$fit, c(list(data = data), SAMPLE_KWARGS))
  out <- capture.output(print(model))
  expect_match(paste(out, collapse = "\n"), "fitted")
})

test_that("summary() returns a tibble with node and mean columns", {
  skip_on_cran()
  data  <- make_chain_data()
  model <- KaguModel$new(dag = CHAIN_DAG)
  do.call(model$fit, c(list(data = data), SAMPLE_KWARGS))
  df <- model$summary()
  expect_true(tibble::is_tibble(df) || is.data.frame(df))
  expect_true("node" %in% names(df))
  expect_true("mean" %in% names(df))
})

test_that("save and load round-trip preserves DAG and enables effects()", {
  skip_on_cran()
  data  <- make_chain_data()
  model <- KaguModel$new(dag = CHAIN_DAG)
  do.call(model$fit, c(list(data = data), SAMPLE_KWARGS))

  path <- withr::local_tempfile(fileext = ".rds")
  model$save(path)

  loaded <- KaguModel$load(path)
  expect_true(loaded$.fitted)
  expect_equal(loaded$dag, model$dag)
  expect_false(is.null(loaded$data))

  effect <- loaded$effects("x", "y")
  expect_true(effect$samples@.Data |> length() > 0 || length(effect$samples) > 0)
})

test_that("save with include_data=FALSE stores NULL data", {
  skip_on_cran()
  data  <- make_chain_data()
  model <- KaguModel$new(dag = CHAIN_DAG)
  do.call(model$fit, c(list(data = data), SAMPLE_KWARGS))

  path <- withr::local_tempfile(fileext = ".rds")
  model$save(path, include_data = FALSE)

  loaded <- KaguModel$load(path)
  expect_null(loaded$data)
})
