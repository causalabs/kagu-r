test_that("validate_dag accepts a valid DAG", {
  dag <- list(a = c(), b = c("a"), c = c("a"), d = c("b", "c"))
  expect_silent(validate_dag(dag))
})

test_that("validate_dag errors on unknown parent", {
  expect_error(validate_dag(list(x = c("y"))), "unknown parent")
})

test_that("validate_dag errors on cycle", {
  expect_error(validate_dag(list(a = c("b"), b = c("a"))), "cycle")
})

test_that("topological_sort respects parent-before-child ordering", {
  dag   <- list(a = c(), b = c("a"), c = c("a"), d = c("b", "c"))
  order <- topological_sort(dag)

  expect_true(which(order == "a") < which(order == "b"))
  expect_true(which(order == "a") < which(order == "c"))
  expect_true(which(order == "b") < which(order == "d"))
  expect_true(which(order == "c") < which(order == "d"))
})

test_that("topological_sort returns all nodes", {
  dag   <- list(a = c(), b = c("a"), c = c("a"), d = c("b", "c"))
  order <- topological_sort(dag)
  expect_setequal(order, names(dag))
})

test_that("descendants are correct", {
  dag <- list(a = c(), b = c("a"), c = c("a"), d = c("b", "c"))
  expect_setequal(descendants(dag, "a"), c("b", "c", "d"))
  expect_setequal(descendants(dag, "b"), c("d"))
  expect_length(descendants(dag, "d"), 0)
})

test_that("ancestors are correct", {
  dag <- list(a = c(), b = c("a"), c = c("a"), d = c("b", "c"))
  expect_setequal(ancestors(dag, "d"), c("a", "b", "c"))
  expect_setequal(ancestors(dag, "b"), c("a"))
  expect_length(ancestors(dag, "a"), 0)
})

test_that("is_ancestor returns correct booleans", {
  dag <- list(a = c(), b = c("a"), c = c("a"), d = c("b", "c"))
  expect_true(is_ancestor(dag, "a", "d"))
  expect_false(is_ancestor(dag, "d", "a"))
  expect_false(is_ancestor(dag, "b", "c"))
})

test_that("node_depth assigns correct depths", {
  dag   <- list(a = c(), b = c("a"), c = c("a"), d = c("b", "c"))
  depth <- node_depth(dag)
  expect_equal(depth[["a"]], 0L)
  expect_equal(depth[["b"]], 1L)
  expect_equal(depth[["c"]], 1L)
  expect_equal(depth[["d"]], 2L)
})
