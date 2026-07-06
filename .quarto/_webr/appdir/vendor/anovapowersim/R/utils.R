#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x


#' Assert that a number is in the open interval (0, 1)
#'
#' @param x Numeric.
#' @param arg Argument name used in the error message.
#' @keywords internal
#' @noRd
assert_unit_interval <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= 0 || x >= 1) {
    stop("`", arg, "` must be a single number in (0, 1). Got: ",
         paste(x, collapse = ", "), call. = FALSE)
  }
  invisible(x)
}


#' Coerce and assert that a value is in the open interval (0, 1)
#'
#' Accepts numeric scalars and numeric-looking character scalars. This keeps
#' low-level helpers usable with formatted effect-size output such as ".310"
#' while preserving strict interval validation.
#'
#' @param x Numeric or character scalar.
#' @param arg Argument name used in the error message.
#' @keywords internal
#' @noRd
as_unit_interval <- function(x, arg) {
  x_num <- if (is.character(x) && length(x) == 1L) {
    suppressWarnings(as.numeric(trimws(x)))
  } else {
    x
  }
  assert_unit_interval(x_num, arg)
  x_num
}


#' Create internal, collision-proof cell labels
#'
#' @keywords internal
#' @noRd
make_cell_labels <- function(prefix, n) {
  sprintf(".anovapowersim_%s_%03d", prefix, seq_len(n))
}

