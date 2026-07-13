#' Calculate the sample size needed for target ANOVA power
#'
#' Calculation-only search for the per-between-cell sample size needed to reach
#' a requested power for a balanced factorial ANOVA design. Unlike
#' [power_n()], this function does not run simulations, fit ANOVA models, or
#' call `car`; numerator degrees of freedom, denominator degrees of freedom,
#' noncentrality, and power are computed analytically from the balanced design.
#'
#' @section Lifecycle:
#' \ifelse{html}{\out{<a href="https://lifecycle.r-lib.org/articles/stages.html#experimental"><img src="https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg" alt="[Experimental]"></a>}}{\strong{Experimental}}
#'
#' `power_n_calc()` is experimental while the analytic search API and reporting
#' format are refined.
#'
#' @param between Named integer vector of between-subject factor level counts,
#'   e.g. `c(group = 2)`. Use `NULL` for no between-subject factors.
#' @param within Named integer vector of within-subject factor level counts,
#'   e.g. `c(time = 3, condition = 4)`. Use `NULL` for no within-subject
#'   factors.
#' @param term Character scalar naming the ANOVA term to test, e.g.
#'   `"group:time"`. Interaction terms are order-insensitive; `"time:group"`
#'   resolves to `"group:time"` when that is the design's factor order.
#' @param target_pes Target partial eta squared for `term`.
#' @param power Desired target power.
#' @param alpha Significance threshold.
#' @param n_start Starting sample size per between-subject cell. If `NULL`,
#'   starts from the smallest analytically valid value.
#' @param n_max Maximum sample size per between-subject cell.
#' @param gpower Logical; if `TRUE`, use the G*Power-style noncentrality
#'   convention `lambda = total_n * f^2`. The default `FALSE` uses
#'   `lambda = den_df * f^2`.
#'
#' @return An `anovapowersim_curve` object with `n_needed` and
#'   `total_n_needed`. The `$results` tibble contains `n_per_cell`, `total_n`,
#'   `n_sims`, numerator and denominator degrees of freedom (`num_df`,
#'   `den_df`), the noncentrality parameter (`ncp`), calculated power
#'   (`power_calc`), and simulated power (`power_sim`). For `power_n_calc()`,
#'   `n_sims` and `power_sim` are always `NA`.
#'
#' @examples
#' power_n_calc(
#'   between = c(cond = 2),
#'   within = c(stim = 4),
#'   term = "cond:stim",
#'   target_pes = 0.14,
#'   power = 0.90
#' )
#'
#' @export
power_n_calc <- function(between = NULL,
                         within = NULL,
                         term,
                         target_pes,
                         power = 0.90,
                         alpha = 0.05,
                         n_start = NULL,
                         n_max = 5000,
                         gpower = FALSE) {
  setup <- prepare_power_n_calc_inputs(
    between = between,
    within = within,
    term = term,
    target_pes = target_pes,
    power = power,
    alpha = alpha,
    n_start = n_start,
    n_max = n_max,
    gpower = gpower
  )

  curve <- analytic_power_search(
    spec = setup$spec,
    term = setup$term,
    target_pes = target_pes,
    target_power = power,
    alpha = alpha,
    n_start = setup$n_start,
    n_max = setup$n_max,
    gpower = setup$gpower
  )

  n_needed <- estimate_calc_n_needed(curve, target = power)
  total_n_needed <- if (is.na(n_needed)) {
    NA_integer_
  } else {
    as.integer(n_needed * max(1L, setup$spec$n_between_cells))
  }

  structure(
    list(
      results = curve,
      term = setup$term,
      power = power,
      alpha = alpha,
      target_pes = target_pes,
      scale_factor = NA_real_,
      n_sims = NA_integer_,
      n_needed = n_needed,
      total_n_needed = total_n_needed,
      gpower = setup$gpower,
      ss_type = NULL,
      design = setup$spec,
      call = match.call()
    ),
    class = "anovapowersim_curve"
  )
}


#' @keywords internal
#' @noRd
prepare_power_n_calc_inputs <- function(between, within, term, target_pes,
                                        power, alpha, n_start, n_max,
                                        gpower) {
  spec <- balanced_anova_design(between = between, within = within)
  term <- resolve_design_term(term, spec)
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
  assert_unit_interval(power, "power")
  if (is.finite(power) && power < 0.90) {
    warning("Power greater than or equal to .90 is recommended.",
            call. = FALSE,
            immediate. = TRUE)
  }
  assert_unit_interval(alpha, "alpha")
  if (!is.logical(gpower) || length(gpower) != 1L || is.na(gpower)) {
    stop("`gpower` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(n_max) || length(n_max) != 1L ||
      n_max < 1 || n_max != as.integer(n_max)) {
    stop("`n_max` must be a single positive integer.", call. = FALSE)
  }

  n_min <- minimum_analytic_n(spec)
  if (is.null(n_start)) {
    n_start <- n_min
  } else if (!is.numeric(n_start) || length(n_start) != 1L ||
             n_start < 1 || n_start != as.integer(n_start)) {
    stop("`n_start` must be a single positive integer.", call. = FALSE)
  }
  n_start <- as.integer(n_start)
  n_max <- as.integer(n_max)
  if (n_start < n_min) {
    stop(
      "`n_start` is too small for analytic power calculation. ",
      "Use `n_start >= ", n_min, "` so the denominator degrees of freedom ",
      "are positive.",
      call. = FALSE
    )
  }
  if (n_max < n_start) {
    stop("`n_max` must be greater than or equal to `n_start`.",
         call. = FALSE)
  }

  list(
    spec = spec,
    term = term,
    n_start = n_start,
    n_max = n_max,
    gpower = gpower
  )
}


#' @keywords internal
#' @noRd
minimum_analytic_n <- function(spec) {
  2L
}


#' @keywords internal
#' @noRd
analytic_power_row <- function(n, spec, term, target_pes, alpha, gpower) {
  dfs <- analytic_term_dfs(spec = spec, term = term, n = n)
  total_n <- as.integer(n * max(1L, spec$n_between_cells))
  ncp <- ncp_from_pes(
    pes = target_pes,
    total_n = total_n,
    den_df = dfs$den_df,
    gpower = gpower
  )
  power_calc <- stats::pf(
    stats::qf(1 - alpha, dfs$num_df, dfs$den_df),
    dfs$num_df,
    dfs$den_df,
    ncp = ncp,
    lower.tail = FALSE
  )

  data.frame(
    n_per_cell = as.integer(n),
    total_n = total_n,
    n_sims = NA_integer_,
    num_df = dfs$num_df,
    den_df = dfs$den_df,
    ncp = ncp,
    power_calc = power_calc,
    power_sim = NA_real_,
    check.names = FALSE
  )
}


#' @keywords internal
#' @noRd
analytic_term_dfs <- function(spec, term, n) {
  term_factors <- strsplit(term, ":", fixed = TRUE)[[1L]]
  num_df <- prod(spec$level_counts[term_factors] - 1L)
  total_n <- as.integer(n * max(1L, spec$n_between_cells))

  within_factors <- intersect(term_factors, spec$within)
  if (length(within_factors)) {
    within_term_df <- prod(spec$level_counts[within_factors] - 1L)
    den_df <- (total_n - spec$n_between_cells) * within_term_df
  } else {
    den_df <- total_n - spec$n_between_cells
  }

  list(
    num_df = as.numeric(num_df),
    den_df = as.numeric(den_df)
  )
}


#' @keywords internal
#' @noRd
estimate_calc_n_needed <- function(curve, target) {
  curve <- curve[order(curve$n_per_cell), , drop = FALSE]
  above <- which(curve$power_calc >= target)
  if (length(above) == 0L) return(NA_integer_)
  as.integer(curve$n_per_cell[[above[1L]]])
}

#' @keywords internal
#' @noRd
analytic_power_search <- function(spec, term, target_pes, target_power, alpha,
                                  n_start, n_max, gpower) {
  visited <- list()

  run_one <- function(n) {
    key <- as.character(n)
    if (!is.null(visited[[key]])) return(visited[[key]])
    row <- analytic_power_row(
      n = n,
      spec = spec,
      term = term,
      target_pes = target_pes,
      alpha = alpha,
      gpower = gpower
    )
    visited[[key]] <<- row
    row
  }

  lo <- NA_integer_
  hi <- NA_integer_
  n <- as.integer(n_start)
  repeat {
    row <- run_one(n)
    if (is.finite(row$power_calc) && row$power_calc >= target_power) {
      hi <- n
      break
    }
    lo <- n
    if (n >= n_max) break
    n <- min(n_max, max(n + 1L, n * 2L))
  }

  if (!is.na(hi) && !is.na(lo)) {
    while (hi > lo + 1L) {
      mid <- as.integer(floor((lo + hi) / 2L))
      row <- run_one(mid)
      if (is.finite(row$power_calc) && row$power_calc >= target_power) {
        hi <- mid
      } else {
        lo <- mid
      }
    }
  }

  out <- do.call(rbind, visited)
  out <- out[order(out$n_per_cell), , drop = FALSE]
  out <- out[!duplicated(out$n_per_cell), , drop = FALSE]
  row.names(out) <- NULL
  out
}
