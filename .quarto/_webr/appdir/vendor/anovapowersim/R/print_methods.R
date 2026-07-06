#' Print an anovapowersim power curve
#'
#' Compact one-screen summary: target, term, effective effect size,
#' estimated per-cell and total sample sizes, and the first and last rows of
#' the power curve.
#'
#' @param x An `anovapowersim_curve` object.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.anovapowersim_curve <- function(x, ...) {
  cat("<anovapowersim_curve>\n")
  cat("  term:          '", x$term, "'\n", sep = "")
  cat("  target power:  ",
      if (is.null(x$power) || !is.finite(x$power)) "<not specified>"
      else sprintf("%.3f", x$power),
      "\n", sep = "")
  cat("  alpha:         ", format(x$alpha), "\n", sep = "")
  cat("  effect size:   pes = ", format(round(x$target_pes, 4)), sep = "")
  if (!is.null(x$scale_factor) && is.finite(x$scale_factor) &&
      !isTRUE(all.equal(x$scale_factor, 1))) {
    cat("  (rescaled, k = ", format(round(x$scale_factor, 3)), ")", sep = "")
  }
  cat("\n")
  if (inherits(x$design, "anovapowersim_unbalanced_design_spec")) {
    cat("  n values:      explicit unbalanced cell counts\n", sep = "")
  } else {
    cat("  n values:      ", nrow(x$results), " per-cell sample sizes visited\n", sep = "")
  }
  cat("  sims per cell size: ", x$n_sims, "\n", sep = "")
  if (isTRUE(x$gpower)) {
    cat("  G*Power convention: TRUE\n", sep = "")
  }
  if (!is.null(x$ss_type)) {
    cat("  SS type:       ", x$ss_type, "\n", sep = "")
  }

  if (inherits(x$design, "anovapowersim_unbalanced_design_spec")) {
    cat("  design:        unbalanced between-subject cells\n", sep = "")
  } else {
    cat("  n needed for between-subjects cell: ",
        if (is.na(x$n_needed)) "<not reached>" else x$n_needed, "\n",
        sep = "")
    cat("  total N needed: ",
        if (is.na(x$total_n_needed)) "<not reached>" else x$total_n_needed, "\n",
        sep = "")
  }
  cat("\n")
  print(format_power_results(x$results), row.names = FALSE)
  invisible(x)
}


#' Summarise an anovapowersim power curve
#'
#' Returns the full `$results` tibble along with a small header containing
#' the target, effective effect size, and estimated `n_needed`.
#'
#' @param object An `anovapowersim_curve` object.
#' @param ... Unused.
#'
#' @return A list with elements `header` (named character) and `curve`
#'   (tibble), invisibly; printed to console as well.
#' @export
summary.anovapowersim_curve <- function(object, ...) {
  header <- c(
    term         = object$term,
    target_power = if (is.null(object$power) || !is.finite(object$power)) {
      "<not specified>"
    } else sprintf("%.3f", object$power),
    alpha        = sprintf("%.3f", object$alpha),
    target_pes   = sprintf("%.4f", object$target_pes),
    scale_factor = if (is.null(object$scale_factor) ||
                       !is.finite(object$scale_factor)) {
      "<not recorded>"
    } else sprintf("%.3f", object$scale_factor),
    n_sims       = as.character(object$n_sims),
    n_needed_between_subjects_cell =
      if (is.na(object$n_needed)) "<not reached>"
      else as.character(object$n_needed),
    total_n_needed = if (is.na(object$total_n_needed)) "<not reached>"
                     else as.character(object$total_n_needed)
  )
  if (isTRUE(object$gpower)) {
    header <- append(header, c(gpower_convention = "TRUE"), after = 4L)
  }
  if (!is.null(object$ss_type)) {
    header <- append(header, c(ss_type = object$ss_type), after = 4L)
  }
  out <- list(header = header, curve = object$results)
  cat("anovapowersim power simulation summary\n")
  cat("----------------------------------\n")
  for (nm in names(header)) {
    cat(sprintf("  %-13s %s\n", paste0(nm, ":"), header[[nm]]))
  }
  cat("\nPower curve:\n")
  print(format_power_results(object$results), row.names = FALSE)
  invisible(out)
}


#' @keywords internal
#' @noRd
format_power_results <- function(x) {
  x <- as.data.frame(x)
  numeric_3dp <- intersect(c("ncp", "power_calc", "power_sim"), names(x))
  for (nm in numeric_3dp) {
    x[[nm]] <- ifelse(is.na(x[[nm]]), NA_character_, sprintf("%.3f", x[[nm]]))
  }
  x
}

