# Vendored, dependency-free adaptation of anovapowersim's calculation-only
# achieved-power and sensitivity functions. The upstream package returns
# tibbles; this Shinylive build uses base data frames so it only needs Shiny.

#' Calculate achieved ANOVA power at a fixed sample size
#'
#' @param between Named integer vector of between-subject factor level counts.
#' @param within Named integer vector of within-subject factor level counts.
#' @param term Character scalar naming the ANOVA term to test.
#' @param target_pes Target partial eta squared for `term`.
#' @param n Number of subjects per between-subject cell. For a purely
#'   within-subject design, this is the total sample size.
#' @param alpha Significance threshold.
#' @param gpower Whether to use the G*Power-style noncentrality convention.
#' @param epsilon Population nonsphericity correction.
#'
#' @return An `anovapowersim_achieved_power` object.
power_achieved_calc <- function(between = NULL,
                                within = NULL,
                                term,
                                target_pes,
                                n,
                                alpha = 0.05,
                                gpower = FALSE,
                                epsilon = 1) {
  setup <- prepare_fixed_calc_inputs(
    between = between,
    within = within,
    term = term,
    n = n,
    alpha = alpha,
    gpower = gpower,
    epsilon = epsilon
  )
  validate_target_pes(target_pes)

  row <- analytic_power_row(
    n = setup$n,
    spec = setup$spec,
    term = setup$term,
    target_pes = target_pes,
    alpha = alpha,
    gpower = setup$gpower,
    epsilon = setup$epsilon
  )
  calculated_power <- row$power_calc[[1L]]

  structure(
    list(
      results = row,
      term = setup$term,
      alpha = alpha,
      target_pes = target_pes,
      n = setup$n,
      total_n = row$total_n[[1L]],
      n_sims = NA_integer_,
      achieved_power = calculated_power,
      calculated_power = calculated_power,
      calculation_only = TRUE,
      gpower = setup$gpower,
      epsilon = setup$epsilon,
      covariance = NULL,
      custom_covariance = FALSE,
      ss_type = NULL,
      design = setup$spec,
      call = match.call()
    ),
    class = "anovapowersim_achieved_power"
  )
}


#' Calculate ANOVA effect-size sensitivity at a fixed sample size
#'
#' @inheritParams power_achieved_calc
#' @param power Desired target power.
#' @param pes_min Lower bound of the partial eta-squared search interval.
#' @param pes_max Upper bound of the partial eta-squared search interval.
#' @param pes_tol Maximum width of the final calculated partial eta-squared
#'   bracket.
#'
#' @return An `anovapowersim_sensitivity` object.
power_sensitivity_calc <- function(between = NULL,
                                   within = NULL,
                                   term,
                                   n,
                                   power = 0.90,
                                   alpha = 0.05,
                                   pes_min = 1e-6,
                                   pes_max = 0.99,
                                   pes_tol = 0.001,
                                   gpower = FALSE,
                                   epsilon = 1) {
  setup <- prepare_fixed_calc_inputs(
    between = between,
    within = within,
    term = term,
    n = n,
    alpha = alpha,
    gpower = gpower,
    epsilon = epsilon
  )
  assert_unit_interval(power, "power")
  if (power < 0.90) {
    warning(
      "Power greater than or equal to .90 is recommended.",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  bounds <- validate_pes_search_controls(pes_min, pes_max, pes_tol)

  search <- analytic_effect_search(
    spec = setup$spec,
    term = setup$term,
    n = setup$n,
    target_power = power,
    alpha = alpha,
    gpower = setup$gpower,
    epsilon = setup$epsilon,
    pes_min = bounds$pes_min,
    pes_max = bounds$pes_max,
    pes_tol = bounds$pes_tol
  )

  if (!search$reached) {
    warning(
      sprintf(
        paste(
          "Target power %.3f was not reached by `pes_max = %.6g`.",
          "Increase `pes_max` and rerun the search."
        ),
        power,
        bounds$pes_max
      ),
      call. = FALSE
    )
  } else if (!search$converged) {
    warning(
      sprintf(
        paste(
          "The calculated-power sensitivity search reached its iteration cap before",
          "`pes_tol = %.6g`; reporting the best calculated upper bracket",
          "(pes = %.6g)."
        ),
        bounds$pes_tol,
        search$pes_upper
      ),
      call. = FALSE
    )
  }

  structure(
    list(
      results = search$results,
      term = setup$term,
      power = power,
      alpha = alpha,
      n = setup$n,
      total_n = as.integer(setup$n * max(1L, setup$spec$n_between_cells)),
      n_sims = NA_integer_,
      pes_needed = if (search$reached) search$pes_upper else NA_real_,
      pes_lower = search$pes_lower,
      pes_upper = search$pes_upper,
      pes_min = bounds$pes_min,
      pes_max = bounds$pes_max,
      pes_tol = bounds$pes_tol,
      bracket_width = search$bracket_width,
      converged = search$converged,
      iterations = search$iterations,
      initial_pes = NA_real_,
      calculation_only = TRUE,
      gpower = setup$gpower,
      epsilon = setup$epsilon,
      covariance = NULL,
      custom_covariance = FALSE,
      ss_type = NULL,
      design = setup$spec,
      call = match.call()
    ),
    class = "anovapowersim_sensitivity"
  )
}


validate_target_pes <- function(target_pes) {
  assert_unit_interval(target_pes, "target_pes")
  if (target_pes == 0.06) {
    warning(
      paste(
        "It looks like you are using a rule-of-thumb \"medium\" effect size.",
        "This might overestimate the true effect size, rendering your study",
        "underpowered. Consider basing your power calculations on previous",
        "research or empirically-derived guidelines."
      ),
      call. = FALSE,
      immediate. = TRUE
    )
  }
  invisible(target_pes)
}


prepare_fixed_calc_inputs <- function(between, within, term, n, alpha,
                                      gpower, epsilon) {
  spec <- balanced_anova_design(between = between, within = within)
  term <- resolve_design_term(term, spec)
  assert_unit_interval(alpha, "alpha")
  if (!is.logical(gpower) || length(gpower) != 1L || is.na(gpower)) {
    stop("`gpower` must be TRUE or FALSE.", call. = FALSE)
  }
  warn_gpower_within_term_df(gpower = gpower)
  epsilon <- validate_analytic_epsilon(epsilon, spec = spec, term = term)
  n <- validate_fixed_analytic_n(n, spec)

  list(
    spec = spec,
    term = term,
    n = n,
    gpower = gpower,
    epsilon = epsilon
  )
}


validate_fixed_analytic_n <- function(n, spec) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 1 || n != as.integer(n)) {
    stop("`n` must be a single positive integer.", call. = FALSE)
  }
  n <- as.integer(n)
  n_min <- minimum_analytic_n(spec)
  if (n < n_min) {
    stop(
      "`n` is too small for calculated power. Use `n >= ",
      n_min,
      "` so the denominator degrees of freedom are positive.",
      call. = FALSE
    )
  }
  n
}


validate_pes_search_controls <- function(pes_min, pes_max, pes_tol) {
  assert_unit_interval(pes_min, "pes_min")
  assert_unit_interval(pes_max, "pes_max")
  if (pes_min >= pes_max) {
    stop("`pes_min` must be smaller than `pes_max`.", call. = FALSE)
  }
  if (!is.numeric(pes_tol) || length(pes_tol) != 1L ||
      !is.finite(pes_tol) || pes_tol <= 0) {
    stop("`pes_tol` must be a single positive finite number.", call. = FALSE)
  }
  list(
    pes_min = as.numeric(pes_min),
    pes_max = as.numeric(pes_max),
    pes_tol = as.numeric(pes_tol)
  )
}


analytic_effect_search <- function(spec, term, n, target_power, alpha, gpower,
                                   epsilon, pes_min, pes_max, pes_tol,
                                   max_iter = 100L) {
  visited <- new.env(parent = emptyenv())

  run_one <- function(pes) {
    key <- sprintf("%.17g", pes)
    if (exists(key, envir = visited, inherits = FALSE)) {
      return(get(key, envir = visited, inherits = FALSE))
    }
    row <- analytic_power_row(
      n = n,
      spec = spec,
      term = term,
      target_pes = pes,
      alpha = alpha,
      gpower = gpower,
      epsilon = epsilon
    )
    row <- cbind(target_pes = as.numeric(pes), row)
    assign(key, row, envir = visited)
    row
  }

  lower <- run_one(pes_min)
  lower_power <- lower$power_calc[[1L]]
  if (is.finite(lower_power) && lower_power >= target_power) {
    return(analytic_effect_search_result(
      visited,
      reached = TRUE,
      converged = TRUE,
      pes_lower = NA_real_,
      pes_upper = pes_min,
      bracket_width = 0,
      iterations = 0L
    ))
  }

  upper <- run_one(pes_max)
  upper_power <- upper$power_calc[[1L]]
  if (!is.finite(upper_power) || upper_power < target_power) {
    return(analytic_effect_search_result(
      visited,
      reached = FALSE,
      converged = FALSE,
      pes_lower = pes_max,
      pes_upper = NA_real_,
      bracket_width = NA_real_,
      iterations = 0L
    ))
  }

  lo <- pes_min
  hi <- pes_max
  iter <- 0L
  while ((hi - lo) > pes_tol && iter < max_iter) {
    mid <- lo + (hi - lo) / 2
    if (mid <= lo || mid >= hi) {
      break
    }
    row <- run_one(mid)
    if (is.finite(row$power_calc[[1L]]) &&
        row$power_calc[[1L]] >= target_power) {
      hi <- mid
    } else {
      lo <- mid
    }
    iter <- iter + 1L
  }
  width <- hi - lo

  analytic_effect_search_result(
    visited,
    reached = TRUE,
    converged = width <= pes_tol,
    pes_lower = lo,
    pes_upper = hi,
    bracket_width = width,
    iterations = iter
  )
}


analytic_effect_search_result <- function(visited, reached, converged,
                                          pes_lower, pes_upper, bracket_width,
                                          iterations) {
  keys <- ls(visited, all.names = TRUE)
  results <- do.call(rbind, lapply(keys, get, envir = visited))
  results <- results[order(results$target_pes), , drop = FALSE]
  results <- results[!duplicated(results$target_pes), , drop = FALSE]
  row.names(results) <- NULL
  list(
    results = results,
    reached = reached,
    converged = converged,
    pes_lower = pes_lower,
    pes_upper = pes_upper,
    bracket_width = bracket_width,
    iterations = iterations
  )
}
