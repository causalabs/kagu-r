#' Save and load fitted KaguModel objects
#'
#' @name io
NULL

#' Save a fitted model to disk
#'
#' @param model A fitted `KaguModel`.
#' @param path File path to write to (e.g. `"model.rds"`).
#' @param include_data Logical — whether to save the training data alongside
#'   the model. Defaults to `TRUE` so that `$effects()` works immediately
#'   after loading without needing to reattach the data.
#' @return Invisible `NULL`.
#' @export
kagu_save <- function(model, path, include_data = TRUE) {
  payload <- list(
    version    = 1L,
    dag        = model$dag,
    mechanisms = model$mechanisms,
    traces     = model$traces,
    data       = if (include_data) model$data else NULL
  )
  saveRDS(payload, file = path)
  invisible(NULL)
}

#' Load a model saved with [kagu_save()]
#'
#' @param path File path to read from.
#' @return A `KaguModel` with `dag`, `mechanisms`, `traces`, and (if saved)
#'   `data` restored. The model is marked as fitted.
#' @export
kagu_load <- function(path) {
  payload <- readRDS(path)

  version <- payload[["version"]] %||% 0L
  if (version != 1L) {
    stop(sprintf(
      "Saved model version %d is not compatible with the current version 1.",
      version
    ))
  }

  model              <- KaguModel$new(payload[["dag"]], payload[["mechanisms"]])
  model$traces       <- payload[["traces"]]
  model$data         <- payload[["data"]]  # NULL if saved with include_data = FALSE
  model$.fitted      <- TRUE

  model
}

# Null-coalescing operator (backport for R < 4.4)
`%||%` <- function(x, y) if (!is.null(x)) x else y
