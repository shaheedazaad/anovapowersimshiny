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
  calc_only <- is_calculation_only_curve(x)
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
  if (calc_only) {
    cat("  calculation:   calculated power only\n", sep = "")
  } else {
    cat("  sims per cell size: ", x$n_sims, "\n", sep = "")
  }
  if (!calc_only && !is.null(x$sim_correction_resolved)) {
    cat("  simulated test: ", simulated_test_label(
      x$sim_correction, x$sim_correction_resolved
    ), "\n", sep = "")
  }
  if (!calc_only && !is.null(x$custom_means_pattern)) {
    cat("  means pattern: ",
        if (isTRUE(x$custom_means_pattern)) "custom" else
          "default linear/Kronecker",
        "\n", sep = "")
  }
  if (isTRUE(x$gpower)) {
    cat("  G*Power convention: TRUE\n", sep = "")
  }
  if (!is.null(x$epsilon) && is.finite(x$epsilon) && x$epsilon < 1) {
    cat("  epsilon:       ", format(x$epsilon), "\n", sep = "")
  }
  if (isTRUE(x$custom_covariance)) {
    cat(
      "  covariance:    custom ", nrow(x$covariance), " x ",
      ncol(x$covariance), " within-subject matrix\n",
      sep = ""
    )
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
  calc_only <- is_calculation_only_curve(object)
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
    calculation  = if (calc_only) "calculated power only" else "simulation",
    n_sims       = if (calc_only) "<not run>" else as.character(object$n_sims),
    n_needed_between_subjects_cell =
      if (is.na(object$n_needed)) "<not reached>"
      else as.character(object$n_needed),
    total_n_needed = if (is.na(object$total_n_needed)) "<not reached>"
                     else as.character(object$total_n_needed)
  )
  if (isTRUE(object$gpower)) {
    header <- append(header, c(gpower_convention = "TRUE"), after = 4L)
  }
  if (!calc_only && !is.null(object$custom_means_pattern)) {
    header <- append(
      header,
      c("means pattern" = if (isTRUE(object$custom_means_pattern)) {
        "custom"
      } else {
        "default linear/Kronecker"
      }),
      after = 4L
    )
  }
  if (!calc_only && !is.null(object$sim_correction_resolved)) {
    header <- append(
      header,
      c(simulated_test = simulated_test_label(
        object$sim_correction, object$sim_correction_resolved
      )),
      after = 4L
    )
  }
  if (!is.null(object$epsilon) && is.finite(object$epsilon) &&
      object$epsilon < 1) {
    header <- append(
      header,
      c(epsilon = sprintf("%.4f", object$epsilon)),
      after = 4L
    )
  }
  if (isTRUE(object$custom_covariance)) {
    header <- append(
      header,
      c(covariance = paste0(
        "custom ", nrow(object$covariance), " x ", ncol(object$covariance),
        " within-subject matrix"
      )),
      after = 4L
    )
  }
  if (!is.null(object$ss_type)) {
    header <- append(header, c(ss_type = object$ss_type), after = 4L)
  }
  out <- list(header = header, curve = object$results)
  cat(if (calc_only) {
    "anovapowersim calculated-power summary\n"
  } else {
    "anovapowersim power simulation summary\n"
  })
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
  if ("target_pes" %in% names(x)) {
    x$target_pes <- ifelse(
      is.na(x$target_pes), NA_character_, sprintf("%.6f", x$target_pes)
    )
  }
  numeric_3dp <- intersect(
    c("ncp", "power_calc", "power_sim"),
    names(x)
  )
  for (nm in numeric_3dp) {
    x[[nm]] <- ifelse(is.na(x[[nm]]), NA_character_, sprintf("%.3f", x[[nm]]))
  }
  x
}


#' @keywords internal
#' @noRd
is_calculation_only_curve <- function(x) {
  is.null(x$n_sims) || length(x$n_sims) != 1L || is.na(x$n_sims)
}


#' @keywords internal
#' @noRd
simulated_test_label <- function(sim_correction, sim_correction_resolved) {
  test <- if (identical(sim_correction_resolved, "GG")) {
    "GG-corrected"
  } else {
    "uncorrected"
  }
  source <- if (identical(sim_correction, "auto")) {
    "auto"
  } else if (identical(sim_correction, "GG") &&
             identical(sim_correction_resolved, "none")) {
    "GG requested; no applicable multi-df within component"
  } else {
    "user-set"
  }
  paste0(test, " (", source, ")")
}
