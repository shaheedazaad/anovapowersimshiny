#' Simulate ANOVA power from a balanced factorial design
#'
#' Simulation-based power estimation for balanced factorial designs under
#' sphericity. Users specify the between- and within-subject factors, the ANOVA
#' term to test, a target partial eta squared, and explicit sample sizes. The
#' function creates a default contrast pattern for the target term, scales it to
#' the requested partial eta squared, simulates datasets, refits
#' `stats::aov()`, and estimates power by counting `p < alpha`.
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
#' @param n_range Integer vector of sample sizes per between-subject cell. For
#'   pure within-subject designs, this is the total sample size.
#' @param n_sims Number of simulated datasets per sample size.
#' @param alpha Significance threshold.
#' @param ss_type Sums-of-squares type for the tested ANOVA term. `"III"` is
#'   the default for order-invariant tests in unbalanced designs. Use `"I"` to
#'   reproduce sequential `stats::aov()` tests.
#' @param gpower Logical; if `TRUE`, calibrate means to the G*Power-style
#'   noncentrality convention `lambda = total_n * f^2`. The default `FALSE`
#'   calibrates the empirical reference dataset to `target_pes`, equivalent to
#'   `lambda = den_df * f^2` for the fitted ANOVA.
#' @param progress Logical; if `TRUE`, show a text progress bar.
#' @param parallel Logical; if `TRUE`, run simulations for each sample size via
#'   the `future` ecosystem.
#' @param cores Optional positive integer number of cores to use when
#'   `parallel = TRUE`. If `NULL`, uses one fewer than the number of available
#'   cores, with a minimum of one.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return An `anovapowersim_curve` object. The `$results` tibble contains
#'   `n_per_cell`, `total_n`, `n_sims`, numerator and denominator degrees of
#'   freedom (`num_df`, `den_df`), the noncentrality parameter (`ncp`),
#'   calculated power (`power_calc`), and simulated power (`power_sim`).
#'
#' @section Examples:
#' ```{r, eval = FALSE}
#' power_curve(
#'   between = c(cond = 2),
#'   within = c(stim = 2),
#'   term = "cond:stim",
#'   target_pes = 0.14,
#'   n_range = c(16, 20, 23, 28), # n per between-subject cell
#'   n_sims = 1000,
#'   seed = 123
#' )
#'
#' power_curve(
#'   between = c(group = 2),
#'   within = c(time = 2),
#'   term = "group:time",
#'   target_pes = 0.14,
#'   n_range = c(12, 16, 20),
#'   n_sims = 5000,
#'   parallel = TRUE,
#'   cores = 4,
#'   seed = 123
#' )
#' ```
#'
#' @export
power_curve <- function(between = NULL,
                        within = NULL,
                        term,
                        target_pes,
                        n_range,
                        n_sims = 10000,
                        alpha = 0.05,
                        ss_type = "III",
                        gpower = FALSE,
                        progress = interactive(),
                        parallel = FALSE,
                        cores = NULL,
                        seed = NULL) {
  sd <- 1
  r <- 0.5
  setup <- prepare_power_curve_inputs(
    between = between,
    within = within,
    term = term,
    target_pes = target_pes,
    n_sims = n_sims,
    alpha = alpha,
    ss_type = ss_type,
    sd = sd,
    r = r,
    gpower = gpower,
    progress = progress,
    parallel = parallel,
    cores = cores
  )
  message_long_serial_run(setup$n_sims, setup$parallel)

  if (!is.null(seed)) set.seed(seed)

  ns <- normalize_n_range(n_range)
  validate_calibration_n(ns, setup$spec, "n_range")
  progress_bar <- make_progress_bar(
    enabled = setup$progress,
    total = if (setup$parallel) length(ns) else length(ns) * setup$n_sims,
    label = "Simulating power"
  )
  on.exit(close_progress_bar(progress_bar), add = TRUE)

  run_one <- function(n) {
    row <- run_design_power_at_n(
      spec = setup$spec,
      term = setup$term,
      target_pes = target_pes,
      n = n,
      n_sims = setup$n_sims,
      alpha = alpha,
      ss_type = setup$ss_type,
      sd = sd,
      r = r,
      gpower = setup$gpower,
      progress_bar = if (setup$parallel) NULL else progress_bar,
      parallel = setup$parallel,
      cores = setup$cores
    )
    if (setup$parallel) tick_progress_bar(progress_bar)
    row
  }

  curve <- purrr::map_dfr(ns, run_one)
  warn_power_disagreement(curve, setup$n_sims)

  structure(
    list(
      results = curve,
      term = setup$term,
      power = NA_real_,
      alpha = alpha,
      target_pes = target_pes,
      scale_factor = NA_real_,
      n_sims = setup$n_sims,
      n_needed = NA_integer_,
      total_n_needed = NA_integer_,
      gpower = setup$gpower,
      ss_type = setup$ss_type,
      design = setup$spec,
      call = match.call()
    ),
    class = "anovapowersim_curve"
  )
}


#' Search for the sample size needed for target ANOVA power
#'
#' Adaptive simulation search for the per-between-cell sample size needed to
#' reach a requested power for a balanced factorial ANOVA design. The search
#' doubles upward from `n_start` until it brackets the target or reaches
#' `n_max`, then refines the bracket using interpolation with midpoint
#' bisection as a fallback.
#'
#' @inheritParams power_curve
#' @param power Desired target power.
#' @param n_start Starting sample size per between-subject cell. If `NULL`,
#'   starts at the smallest value that can support empirical calibration for
#'   the requested design.
#' @param n_max Maximum sample size per between-subject cell.
#' @param tol Acceptable precision above target power. If no simulated value at
#'   or above `power` is also no more than `power + tol`, `power_n()` warns that
#'   the requested precision band was not reached.
#'
#' @return An `anovapowersim_curve` object with `n_needed` and
#'   `total_n_needed`. For `power_n()`, `n_needed` is always an explicitly
#'   simulated `n_per_cell` value, never an interpolated sample size. If the
#'   search reaches target power but no simulated value lands inside
#'   `[power, power + tol]`, `power_n()` reports the smallest explicitly
#'   simulated value at or above target power and warns that the requested
#'   precision band was not reached.
#'
#' @section Examples:
#' ```{r, eval = FALSE}
#' power_n(
#'   between = c(cond = 2),
#'   within = c(stim = 4),
#'   term = "cond:stim",
#'   target_pes = 0.14,
#'   alpha = 0.05,
#'   power = 0.90,
#'   n_sims = 1000, # use 5000+ for a more precise estimate
#'   seed = 123 # for reproducibility
#' )
#' ```
#'
#' @export
power_n <- function(between = NULL,
                    within = NULL,
                    term,
                    target_pes,
                    power = 0.90,
                    n_sims = 10000,
                    alpha = 0.05,
                    ss_type = "III",
                    n_start = NULL,
                    n_max = 1000,
                    tol = 0.03,
                    gpower = FALSE,
                    progress = interactive(),
                    parallel = FALSE,
                    cores = NULL,
                    seed = NULL) {
  sd <- 1
  r <- 0.5
  setup <- prepare_power_curve_inputs(
    between = between,
    within = within,
    term = term,
    target_pes = target_pes,
    n_sims = n_sims,
    alpha = alpha,
    ss_type = ss_type,
    sd = sd,
    r = r,
    gpower = gpower,
    progress = progress,
    parallel = parallel,
    cores = cores
  )
  assert_unit_interval(power, "power")
  if (is.finite(power) && power < 0.90) {
    warning("Power greater than or equal to .90 is recommended.",
            call. = FALSE,
            immediate. = TRUE)
  }
  message_long_serial_run(setup$n_sims, setup$parallel)
  if (!is.numeric(n_max) || length(n_max) != 1L ||
      n_max < 1 || n_max != as.integer(n_max)) {
    stop("`n_max` must be a single positive integer.", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("`tol` must be a positive finite number.", call. = FALSE)
  }

  min_n <- minimum_calibration_n(setup$spec)
  if (is.null(n_start)) {
    ncp_n <- estimate_ncp_n_needed(
      spec = setup$spec,
      term = setup$term,
      target_pes = target_pes,
      target_power = power,
      alpha = alpha,
      ss_type = setup$ss_type,
      sd = sd,
      r = r,
      gpower = setup$gpower,
      n_min = min_n,
      n_max = as.integer(n_max)
    )
    n_start <- if (is.na(ncp_n)) {
      min_n
    } else {
      max(min_n, as.integer(floor(ncp_n * 0.8)))
    }
  } else if (!is.numeric(n_start) || length(n_start) != 1L ||
             n_start < 1 || n_start != as.integer(n_start)) {
    stop("`n_start` must be a single positive integer.", call. = FALSE)
  }
  n_start <- as.integer(n_start)
  n_max <- as.integer(n_max)
  validate_calibration_n(n_start, setup$spec, "n_start")
  if (n_max < min_n) {
    stop(
      "`n_max` is too small for this design. This design has ",
      setup$spec$n_within_cells, " within-subject cell",
      if (setup$spec$n_within_cells == 1L) "" else "s",
      ", so the adaptive search needs `n_max >= ", min_n,
      "` per between-subject cell.",
      call. = FALSE
    )
  }

  if (!is.null(seed)) set.seed(seed)

  progress_bar <- make_progress_bar(
    enabled = setup$progress,
    total = estimate_adaptive_progress_total(n_start, n_max),
    label = "Searching sample size"
  )
  on.exit(close_progress_bar(progress_bar), add = TRUE)

  run_one <- function(n) {
    run_design_power_at_n(
      spec = setup$spec,
      term = setup$term,
      target_pes = target_pes,
      n = n,
      n_sims = setup$n_sims,
      alpha = alpha,
      ss_type = setup$ss_type,
      sd = sd,
      r = r,
      gpower = setup$gpower,
      progress_bar = NULL,
      parallel = setup$parallel,
      cores = setup$cores
    )
  }

  curve <- adaptive_design_search(
    run_one = run_one,
    target = power,
    n_start = n_start,
    n_max = n_max,
    tol = tol,
    progress_bar = progress_bar
  )
  warn_power_disagreement(curve, setup$n_sims)

  n_needed <- estimate_design_n_needed(curve, target = power)
  warn_precision_band_not_reached(
    curve = curve,
    target = power,
    tol = tol,
    n_needed = n_needed
  )
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
      n_sims = setup$n_sims,
      n_needed = n_needed,
      total_n_needed = total_n_needed,
      gpower = setup$gpower,
      ss_type = setup$ss_type,
      design = setup$spec,
      call = match.call()
    ),
    class = "anovapowersim_curve"
  )
}


#' @keywords internal
#' @noRd
cell_counts <- function(...) {
  dots <- list(...)
  nms <- names(dots)
  if (!length(dots)) {
    stop("Enter at least one cell, ending each cell with `n = count`.",
         call. = FALSE)
  }
  if (is.null(nms) || any(!nzchar(nms))) {
    stop("Every value in `...` must be named.", call. = FALSE)
  }

  rows <- list()
  current <- list()
  factor_names <- NULL
  for (i in seq_along(dots)) {
    nm <- nms[[i]]
    value <- dots[[i]]
    if (length(value) != 1L || is.null(value) || is.na(value)) {
      stop("`", nm, "` values must be single non-missing values.",
           call. = FALSE)
    }

    if (identical(nm, "n")) {
      if (!length(current)) {
        stop("Each `n` must follow one or more factor values.", call. = FALSE)
      }
      current$n <- value
      row_factor_names <- names(current)[names(current) != "n"]
      if (is.null(factor_names)) {
        factor_names <- row_factor_names
      } else if (!identical(row_factor_names, factor_names)) {
        stop("Every cell must use the same factor names in the same order.",
             call. = FALSE)
      }
      rows[[length(rows) + 1L]] <- current
      current <- list()
    } else {
      if (nm %in% names(current)) {
        stop("Factor `", nm, "` appears more than once before the next `n`.",
             call. = FALSE)
      }
      current[[nm]] <- value
    }
  }

  if (length(current)) {
    stop("The final cell is missing `n = count`.", call. = FALSE)
  }

  out <- dplyr::bind_rows(rows)
  out$n <- validate_cell_count_values(out$n, "n")
  tibble::as_tibble(out)
}


#' @keywords internal
#' @noRd
power_unbalanced <- function(cell_n,
                             within = NULL,
                             term,
                             target_pes,
                             n_sims = 10000,
                             alpha = 0.05,
                             ss_type = "III",
                             gpower = FALSE,
                             progress = interactive(),
                             parallel = FALSE,
                             cores = NULL,
                             seed = NULL) {
  sd <- 1
  r <- 0.5
  spec <- unbalanced_anova_design(cell_n = cell_n, within = within)
  term <- resolve_design_term(term, spec)
  ss_type <- validate_ss_type(ss_type)
  assert_unit_interval(target_pes, "target_pes")
  assert_unit_interval(alpha, "alpha")
  if (!is.numeric(n_sims) || length(n_sims) != 1L || n_sims < 1) {
    stop("`n_sims` must be a positive integer.", call. = FALSE)
  }
  if (!is.logical(gpower) || length(gpower) != 1L || is.na(gpower)) {
    stop("`gpower` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("`progress` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("`parallel` must be TRUE or FALSE.", call. = FALSE)
  }
  cores <- validate_parallel_cores(cores = cores, parallel = parallel)
  validate_calibration_n(spec$cell_n, spec, "cell_n$n")
  message_long_serial_run(as.integer(n_sims), parallel)
  if (!is.null(seed)) set.seed(seed)

  progress_bar <- make_progress_bar(
    enabled = progress,
    total = if (parallel) 1L else as.integer(n_sims),
    label = "Simulating power"
  )
  on.exit(close_progress_bar(progress_bar), add = TRUE)

  row <- run_unbalanced_power(
    spec = spec,
    term = term,
    target_pes = target_pes,
    n_sims = as.integer(n_sims),
    alpha = alpha,
    ss_type = ss_type,
    sd = sd,
    r = r,
    gpower = gpower,
    progress_bar = if (parallel) NULL else progress_bar,
    parallel = parallel,
    cores = cores
  )
  if (parallel) tick_progress_bar(progress_bar)
  warn_power_disagreement(row, as.integer(n_sims))

  structure(
    list(
      results = row,
      term = term,
      power = NA_real_,
      alpha = alpha,
      target_pes = target_pes,
      scale_factor = NA_real_,
      n_sims = as.integer(n_sims),
      n_needed = NA_integer_,
      total_n_needed = NA_integer_,
      gpower = gpower,
      ss_type = ss_type,
      design = spec,
      call = match.call()
    ),
    class = "anovapowersim_curve"
  )
}


#' Create a balanced factorial ANOVA design specification
#'
#' Builds the design object used by [power_curve()], [design_term_means()], and
#' [simulate_design_dataset()]. This object stores factor names, level counts,
#' generated factor levels, and the between/within cell grids.
#'
#' @param between Named integer vector of between-subject factor level counts,
#'   e.g. `c(group = 2)`. Use `NULL` for no between-subject factors.
#' @param within Named integer vector of within-subject factor level counts,
#'   e.g. `c(time = 3)`. Use `NULL` for no within-subject factors.
#'
#' @return An object of class `anovapowersim_design_spec`.
#'
#' @examples
#' d <- balanced_anova_design(between = c(group = 2), within = c(time = 3))
#' d$between_cells
#' d$within_cells
#'
#' @export
balanced_anova_design <- function(between = NULL, within = NULL) {
  spec <- validate_design_spec(between, within)
  class(spec) <- c("anovapowersim_design_spec", class(spec))
  spec
}


#' Build calibrated default means for a design term
#'
#' Creates the default contrast pattern for one ANOVA term and scales it so an
#' exact reference dataset has the requested partial eta squared under the
#' supplied balanced design assumptions.
#'
#' @param design An `anovapowersim_design_spec` from [balanced_anova_design()].
#' @param term Character scalar naming the ANOVA term to target. Interaction
#'   terms are order-insensitive.
#' @param target_pes Target partial eta squared.
#' @param n Sample size per between-subject cell. For pure within designs, this
#'   is the total sample size.
#' @param sd Common outcome standard deviation.
#' @param r Compound-symmetric correlation among within-subject cells.
#' @param gpower Logical; if `TRUE`, calibrate to the G*Power-style 
#'   noncentrality convention `lambda = total_n * f^2` (using the 'as in Cohen (1988) option for within-subjects designs).
#' @param ss_type Sums-of-squares type for the tested ANOVA term. `"III"` is
#'   the default for order-invariant tests in unbalanced designs. Use `"I"` to
#'   reproduce sequential `stats::aov()` tests.
#'
#' @return A numeric matrix of cell means, with rows indexing between cells and
#'   columns indexing within cells.
#'
#' @examples
#' d <- balanced_anova_design(between = c(group = 2), within = c(time = 2))
#' design_term_means(d, term = "group:time", target_pes = 0.2, n = 20)
#'
#' @export
design_term_means <- function(design, term, target_pes, n, sd = 1, r = 0.5,
                              gpower = FALSE, ss_type = "III") {
  assert_design_spec(design)
  term <- resolve_design_term(term, design)
  ss_type <- validate_ss_type(ss_type)
  assert_unit_interval(target_pes, "target_pes")
  if (!is.numeric(n) || length(n) != 1L || n < 1 || n != as.integer(n)) {
    stop("`n` must be a single positive integer.", call. = FALSE)
  }
  validate_calibration_n(as.integer(n), design, "n")
  if (!is.numeric(sd) || length(sd) != 1L || !is.finite(sd) || sd <= 0) {
    stop("`sd` must be a positive finite number.", call. = FALSE)
  }
  if (!is.numeric(r) || length(r) != 1L || !is.finite(r) || r <= -1 || r >= 1) {
    stop("`r` must be a single finite correlation in (-1, 1).",
         call. = FALSE)
  }
  if (!is.logical(gpower) || length(gpower) != 1L || is.na(gpower)) {
    stop("`gpower` must be TRUE or FALSE.", call. = FALSE)
  }
  calibrate_design_means(
    spec = design,
    term = term,
    target_pes = target_pes,
    n = as.integer(n),
    sd = sd,
    r = r,
    gpower = gpower,
    ss_type = ss_type
  )
}


#' Simulate data from a balanced ANOVA design
#'
#' Generates one long-format dataset from a balanced design. Supply means from
#' [design_term_means()] or any conformable matrix with one row per
#' between-subject cell and one column per within-subject cell.
#'
#' @param design An `anovapowersim_design_spec` from [balanced_anova_design()].
#' @param n Sample size per between-subject cell. For pure within designs, this
#'   is the total sample size.
#' @param means Numeric matrix of population cell means.
#' @param sd Common outcome standard deviation.
#' @param r Compound-symmetric correlation among within-subject cells.
#' @param empirical Logical; if `TRUE`, use `MASS::mvrnorm(empirical = TRUE)`
#'   so the generated sample closely matches the requested means/covariance.
#'
#' @return A tibble ready for `stats::aov()` with columns `id`, factor
#'   columns, and `value`.
#'
#' @examples
#' d <- balanced_anova_design(between = c(group = 2), within = c(time = 2))
#' m <- design_term_means(d, term = "group:time", target_pes = 0.2, n = 20)
#' sim <- simulate_design_dataset(d, n = 20, means = m)
#' head(sim)
#'
#' @export
simulate_design_dataset <- function(design, n, means, sd = 1, r = 0.5,
                                    empirical = FALSE) {
  assert_design_spec(design)
  if (!is.numeric(n) || length(n) != 1L || n < 1 || n != as.integer(n)) {
    stop("`n` must be a single positive integer.", call. = FALSE)
  }
  if (missing(means)) {
    stop("`means` is required. Use design_term_means() to generate defaults.",
         call. = FALSE)
  }
  if (!is.matrix(means)) means <- as.matrix(means)
  expected_dim <- c(design$n_between_cells, design$n_within_cells)
  if (!identical(dim(means), expected_dim)) {
    stop("`means` must be a ", expected_dim[[1]], " x ", expected_dim[[2]],
         " matrix for this design.", call. = FALSE)
  }
  if (any(!is.finite(means))) {
    stop("`means` must contain only finite values.", call. = FALSE)
  }
  if (!is.numeric(sd) || length(sd) != 1L || !is.finite(sd) || sd <= 0) {
    stop("`sd` must be a positive finite number.", call. = FALSE)
  }
  if (!is.numeric(r) || length(r) != 1L || !is.finite(r) || r <= -1 || r >= 1) {
    stop("`r` must be a single finite correlation in (-1, 1).",
         call. = FALSE)
  }
  if (!is.logical(empirical) || length(empirical) != 1L || is.na(empirical)) {
    stop("`empirical` must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(empirical)) {
    validate_calibration_n(as.integer(n), design, "n")
  }

  simulate_balanced_design_data(
    spec = design,
    n = as.integer(n),
    means = means,
    sd = sd,
    r = r,
    empirical = empirical
  )
}


#' @keywords internal
#' @noRd
validate_design_spec <- function(between, within) {
  parse_counts <- function(x, arg) {
    if (is.null(x)) return(stats::setNames(integer(0), character(0)))
    if (!is.numeric(x) || is.null(names(x)) || any(names(x) == "")) {
      stop("`", arg, "` must be a named integer vector of level counts.",
           call. = FALSE)
    }
    bad_names <- names(x)[make.names(names(x)) != names(x)]
    if (length(bad_names)) {
      stop(
        "Factor names in `", arg, "` must be syntactic R names. ",
        "Problem name", if (length(bad_names) == 1L) "" else "s", ": ",
        paste(shQuote(bad_names), collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    if (any(x < 2) || any(x != as.integer(x))) {
      stop("Every entry in `", arg, "` must be an integer >= 2.",
           call. = FALSE)
    }
    stats::setNames(as.integer(x), names(x))
  }

  between <- parse_counts(between, "between")
  within <- parse_counts(within, "within")
  all_names <- c(names(between), names(within))
  if (!length(all_names)) {
    stop("At least one between- or within-subject factor is required.",
         call. = FALSE)
  }
  if (anyDuplicated(all_names)) {
    stop("Factor names must be unique across `between` and `within`.",
         call. = FALSE)
  }

  levels <- stats::setNames(
    lapply(all_names, function(nm) {
      k <- c(between, within)[[nm]]
      factor(paste0(nm, seq_len(k)), levels = paste0(nm, seq_len(k)))
    }),
    all_names
  )

  between_cells <- if (length(between)) {
    tidyr::expand_grid(!!!levels[names(between)])
  } else {
    tibble::tibble(.dummy_between = factor("all"))
  }
  within_cells <- if (length(within)) {
    tidyr::expand_grid(!!!levels[names(within)])
  } else {
    tibble::tibble(.dummy_within = factor("dv"))
  }

  list(
    between = names(between),
    within = names(within),
    factor_names = all_names,
    level_counts = c(between, within),
    levels = levels,
    between_cells = between_cells,
    within_cells = within_cells,
    n_between_cells = if (length(between)) nrow(between_cells) else 1L,
    n_within_cells = if (length(within)) nrow(within_cells) else 1L
  )
}


#' @keywords internal
#' @noRd
unbalanced_anova_design <- function(cell_n, within = NULL) {
  if (!is.data.frame(cell_n)) {
    stop("`cell_n` must be a data frame created by cell_counts() or equivalent.",
         call. = FALSE)
  }
  if (!"n" %in% names(cell_n)) {
    stop("`cell_n` must include an `n` column.", call. = FALSE)
  }
  between <- setdiff(names(cell_n), "n")
  if (!length(between)) {
    stop("`cell_n` must include at least one between-subject factor column.",
         call. = FALSE)
  }
  bad_names <- between[make.names(between) != between]
  if (length(bad_names)) {
    stop(
      "Factor names in `cell_n` must be syntactic R names. Problem name",
      if (length(bad_names) == 1L) "" else "s", ": ",
      paste(shQuote(bad_names), collapse = ", "), ".",
      call. = FALSE
    )
  }

  counts <- validate_cell_count_values(cell_n$n, "cell_n$n")
  cell_n <- tibble::as_tibble(cell_n[, c(between, "n"), drop = FALSE])
  cell_n$n <- counts

  for (nm in between) {
    values <- cell_n[[nm]]
    if (any(is.na(values))) {
      stop("`cell_n` factor columns must not contain missing values.",
           call. = FALSE)
    }
    if (is.list(values)) {
      stop("`cell_n` factor columns must be atomic vectors.", call. = FALSE)
    }
  }

  key <- interaction_key(cell_n, between)
  duplicated_cells <- duplicated(key)
  if (any(duplicated_cells)) {
    stop("`cell_n` contains duplicate between-subject cells.", call. = FALSE)
  }

  between_levels <- stats::setNames(
    lapply(between, function(nm) unique(as.character(cell_n[[nm]]))),
    between
  )
  between_factors <- stats::setNames(
    lapply(between, function(nm) {
      factor(between_levels[[nm]], levels = between_levels[[nm]])
    }),
    between
  )
  between_cells <- tidyr::expand_grid(!!!between_factors)
  expected_key <- interaction_key(between_cells, between)
  missing_key <- setdiff(expected_key, key)
  if (length(missing_key)) {
    stop(
      "`cell_n` must include every combination of the between-subject ",
      "factor levels.",
      call. = FALSE
    )
  }
  extra_key <- setdiff(key, expected_key)
  if (length(extra_key)) {
    stop("`cell_n` contains cells outside the between-subject grid.",
         call. = FALSE)
  }
  cell_n <- cell_n[match(expected_key, key), , drop = FALSE]

  within <- parse_count_vector(within, "within")
  all_names <- c(between, names(within))
  if (anyDuplicated(all_names)) {
    stop("Factor names must be unique across `cell_n` and `within`.",
         call. = FALSE)
  }
  within_levels <- stats::setNames(
    lapply(names(within), function(nm) {
      k <- within[[nm]]
      paste0(nm, seq_len(k))
    }),
    names(within)
  )
  levels <- c(between_levels, within_levels)
  factor_levels <- stats::setNames(
    lapply(names(levels), function(nm) factor(levels[[nm]], levels = levels[[nm]])),
    names(levels)
  )
  within_cells <- if (length(within)) {
    tidyr::expand_grid(!!!factor_levels[names(within)])
  } else {
    tibble::tibble(.dummy_within = factor("dv"))
  }

  spec <- list(
    between = between,
    within = names(within),
    factor_names = all_names,
    level_counts = c(stats::setNames(lengths(between_levels), between), within),
    levels = factor_levels,
    between_cells = between_cells,
    within_cells = within_cells,
    n_between_cells = nrow(between_cells),
    n_within_cells = if (length(within)) nrow(within_cells) else 1L,
    cell_n = as.integer(cell_n$n)
  )
  class(spec) <- c("anovapowersim_unbalanced_design_spec",
                   "anovapowersim_design_spec", class(spec))
  spec
}


#' @keywords internal
#' @noRd
parse_count_vector <- function(x, arg) {
  if (is.null(x)) return(stats::setNames(integer(0), character(0)))
  if (!is.numeric(x) || is.null(names(x)) || any(names(x) == "")) {
    stop("`", arg, "` must be a named integer vector of level counts.",
         call. = FALSE)
  }
  bad_names <- names(x)[make.names(names(x)) != names(x)]
  if (length(bad_names)) {
    stop(
      "Factor names in `", arg, "` must be syntactic R names. Problem name",
      if (length(bad_names) == 1L) "" else "s", ": ",
      paste(shQuote(bad_names), collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (any(x < 2) || any(x != as.integer(x))) {
    stop("Every entry in `", arg, "` must be an integer >= 2.",
         call. = FALSE)
  }
  stats::setNames(as.integer(x), names(x))
}


#' @keywords internal
#' @noRd
validate_cell_count_values <- function(x, arg) {
  if (!is.numeric(x) || any(is.na(x)) || any(x < 1) ||
      any(x != as.integer(x))) {
    stop("`", arg, "` must contain positive integer counts.", call. = FALSE)
  }
  as.integer(x)
}


#' @keywords internal
#' @noRd
interaction_key <- function(data, columns) {
  if (!length(columns)) return(rep(".all", nrow(data)))
  do.call(paste, c(lapply(data[columns], as.character), sep = "\r"))
}


#' @keywords internal
#' @noRd
assert_design_spec <- function(design) {
  if (!inherits(design, "anovapowersim_design_spec")) {
    stop("`design` must be created by balanced_anova_design().",
         call. = FALSE)
  }
  invisible(design)
}


#' @keywords internal
#' @noRd
assert_term_name <- function(term) {
  if (!is.character(term) || length(term) != 1L || is.na(term) || !nzchar(term)) {
    stop("`term` must be a single non-empty character string.", call. = FALSE)
  }
  invisible(term)
}


#' @keywords internal
#' @noRd
resolve_design_term <- function(term, spec) {
  assert_term_name(term)

  requested <- strsplit(term, ":", fixed = TRUE)[[1L]]
  if (any(!nzchar(requested))) {
    stop("`term` must not contain empty ':' components.", call. = FALSE)
  }
  if (anyDuplicated(requested)) {
    stop("`term` must not repeat factor names.", call. = FALSE)
  }

  unknown <- setdiff(requested, spec$factor_names)
  if (length(unknown)) {
    stop("Unknown factor", if (length(unknown) > 1L) "s" else "", " in `term`: ",
         paste(shQuote(unknown), collapse = ", "), ". Available factors: ",
         paste(shQuote(spec$factor_names), collapse = ", "), call. = FALSE)
  }

  paste(spec$factor_names[spec$factor_names %in% requested], collapse = ":")
}


#' @keywords internal
#' @noRd
prepare_power_curve_inputs <- function(between, within, term, target_pes,
                                       n_sims, alpha, ss_type, sd, r, gpower,
                                       progress, parallel, cores) {
  spec <- balanced_anova_design(between = between, within = within)
  term <- resolve_design_term(term, spec)
  ss_type <- validate_ss_type(ss_type)
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
  assert_unit_interval(alpha, "alpha")
  if (!is.numeric(n_sims) || length(n_sims) != 1L || n_sims < 1) {
    stop("`n_sims` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(sd) || length(sd) != 1L || !is.finite(sd) || sd <= 0) {
    stop("`sd` must be a positive finite number.", call. = FALSE)
  }
  if (!is.numeric(r) || length(r) != 1L || !is.finite(r) || r <= -1 || r >= 1) {
    stop("`r` must be a single finite correlation in (-1, 1).",
         call. = FALSE)
  }
  if (!is.logical(gpower) || length(gpower) != 1L || is.na(gpower)) {
    stop("`gpower` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("`progress` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("`parallel` must be TRUE or FALSE.", call. = FALSE)
  }
  cores <- validate_parallel_cores(cores = cores, parallel = parallel)

  list(
    spec = spec,
    term = term,
    n_sims = as.integer(n_sims),
    ss_type = ss_type,
    gpower = gpower,
    progress = progress,
    parallel = parallel,
    cores = cores
  )
}


#' @keywords internal
#' @noRd
normalize_n_range <- function(n_range) {
  if (!is.numeric(n_range) || any(n_range < 1) ||
      any(n_range != as.integer(n_range))) {
    stop("`n_range` must contain positive integers.", call. = FALSE)
  }
  sort(unique(as.integer(n_range)))
}


#' @keywords internal
#' @noRd
validate_ss_type <- function(ss_type) {
  if (!is.character(ss_type) || length(ss_type) != 1L || is.na(ss_type)) {
    stop("`ss_type` must be one of 'III', 'II', or 'I'.", call. = FALSE)
  }
  ss_type <- toupper(ss_type)
  if (!ss_type %in% c("III", "II", "I")) {
    stop("`ss_type` must be one of 'III', 'II', or 'I'.", call. = FALSE)
  }
  ss_type
}


#' @keywords internal
#' @noRd
parallel_worker_helpers <- function(names) {
  env <- new.env(parent = globalenv())
  for (nm in names) {
    value <- get(nm, mode = "function")
    if (is.function(value)) environment(value) <- env
    assign(nm, value, envir = env)
  }
  stats::setNames(lapply(names, get, envir = env), names)
}


#' @keywords internal
#' @noRd
validate_parallel_cores <- function(cores, parallel) {
  available <- as.integer(future::availableCores()[[1L]])
  if (!is.finite(available) || available < 1L) available <- 1L

  if (!is.null(cores) &&
      (!is.numeric(cores) || length(cores) != 1L ||
       !is.finite(cores) || cores < 1 || cores != as.integer(cores))) {
    stop("`cores` must be a single positive integer.", call. = FALSE)
  }

  if (!is.null(cores)) {
    cores <- as.integer(cores)
    if (cores > available) {
      stop(
        "`cores` must not exceed the number of available cores (",
        available, ").",
        call. = FALSE
      )
    }
    return(cores)
  }

  if (!isTRUE(parallel)) return(NULL)

  cores <- max(1L, available - 1L)
  message(
    "`parallel = TRUE` and `cores` was not set; using ",
    cores,
    " of ",
    available,
    " available core",
    if (available == 1L) "" else "s",
    "."
  )
  cores
}


#' @keywords internal
#' @noRd
message_long_serial_run <- function(n_sims, parallel) {
  if (isTRUE(parallel) || n_sims < 5000L) return(invisible(NULL))
  message(
    "`n_sims` is 5000 or more and `parallel = FALSE`; this may take a ",
    "while. Consider setting `parallel = TRUE` to run simulations in parallel."
  )
  invisible(NULL)
}


#' @keywords internal
#' @noRd
run_design_power_at_n <- function(spec, term, target_pes, n, n_sims,
                                  alpha, ss_type, sd, r, gpower,
                                  progress_bar = NULL,
                                  parallel = FALSE, cores = NULL) {
  means <- design_term_means(
    design = spec,
    term = term,
    target_pes = target_pes,
    n = n,
    sd = sd,
    r = r,
    gpower = gpower,
    ss_type = ss_type
  )
  sanity <- sanity_check_term_effect(
    spec = spec,
    term = term,
    target_pes = target_pes,
    n = n,
    means = means,
    sd = sd,
    r = r,
    alpha = alpha,
    gpower = gpower,
    ss_type = ss_type
  )

  helpers <- parallel_worker_helpers(c(
    "simulate_balanced_design_data",
    "fit_design_term_stats",
    "validate_ss_type",
    "fit_design_model",
    "fit_car_term_stats",
    "extract_term_stats",
    "set_sum_contrasts",
    "design_aov_formula",
    "extract_aov_rows",
    "car_between_rows",
    "car_repeated_rows",
    "repeated_measures_wide_data",
    "make_cell_labels",
    "interaction_key",
    "validate_calibration_n",
    "minimum_calibration_n",
    "compound_symmetric_sigma",
    "tick_progress_bar"
  ))
  fit_one_stats <- helpers$fit_design_term_stats
  simulate_one_dataset <- helpers$simulate_balanced_design_data
  tick_progress <- helpers$tick_progress_bar
  simulate_one <- function(i) {
    sim <- simulate_one_dataset(
      spec = spec,
      n = n,
      means = means,
      sd = sd,
      r = r,
      empirical = FALSE
    )
    tick_progress(progress_bar)
    stats <- tryCatch(
      fit_one_stats(sim, spec, term, ss_type = ss_type),
      error = function(e) NULL
    )
    if (is.null(stats) || !is.finite(stats$p_value)) return(NA)
    isTRUE(stats$p_value < alpha)
  }

  successes <- if (parallel) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = cores)
    unlist(
      future.apply::future_lapply(
        seq_len(n_sims),
        simulate_one,
        future.seed = TRUE
      ),
      use.names = FALSE
    )
  } else {
    purrr::map_lgl(seq_len(n_sims), simulate_one)
  }

  failed_count <- sum(is.na(successes))
  valid_count <- n_sims - failed_count
  if (valid_count == 0L) {
    stop(
      "All simulated ANOVA fits failed. This usually indicates an internal ",
      "model-fitting error rather than zero power.",
      call. = FALSE
    )
  }
  success_count <- sum(successes, na.rm = TRUE)
  power <- if (valid_count > 0L) success_count / valid_count else NA_real_

  tibble::tibble(
    n_per_cell = as.integer(n),
    total_n = as.integer(n * max(1L, spec$n_between_cells)),
    n_sims = n_sims,
    num_df = sanity$num_df,
    den_df = sanity$den_df,
    ncp = round(sanity$ncp, 3),
    power_calc = round(sanity$power_calc, 3),
    power_sim = round(power, 3)
  )
}


#' @keywords internal
#' @noRd
run_unbalanced_power <- function(spec, term, target_pes, n_sims, alpha,
                                 ss_type, sd, r, gpower, progress_bar = NULL,
                                 parallel = FALSE, cores = NULL) {
  means <- calibrate_unbalanced_design_means(
    spec = spec,
    term = term,
    target_pes = target_pes,
    sd = sd,
    r = r,
    gpower = gpower,
    ss_type = ss_type
  )
  sanity <- sanity_check_unbalanced_term_effect(
    spec = spec,
    term = term,
    target_pes = target_pes,
    means = means,
    sd = sd,
    r = r,
    alpha = alpha,
    gpower = gpower,
    ss_type = ss_type
  )

  helpers <- parallel_worker_helpers(c(
    "simulate_unbalanced_design_data",
    "fit_design_term_stats",
    "validate_ss_type",
    "fit_design_model",
    "fit_car_term_stats",
    "extract_term_stats",
    "set_sum_contrasts",
    "design_aov_formula",
    "extract_aov_rows",
    "car_between_rows",
    "car_repeated_rows",
    "repeated_measures_wide_data",
    "make_cell_labels",
    "interaction_key",
    "validate_calibration_n",
    "minimum_calibration_n",
    "compound_symmetric_sigma",
    "tick_progress_bar"
  ))
  fit_one_stats <- helpers$fit_design_term_stats
  simulate_one_dataset <- helpers$simulate_unbalanced_design_data
  tick_progress <- helpers$tick_progress_bar
  simulate_one <- function(i) {
    sim <- simulate_one_dataset(
      spec = spec,
      means = means,
      sd = sd,
      r = r,
      empirical = FALSE
    )
    tick_progress(progress_bar)
    stats <- tryCatch(
      fit_one_stats(sim, spec, term, ss_type = ss_type),
      error = function(e) NULL
    )
    if (is.null(stats) || !is.finite(stats$p_value)) return(NA)
    isTRUE(stats$p_value < alpha)
  }

  successes <- if (parallel) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = cores)
    unlist(
      future.apply::future_lapply(
        seq_len(n_sims),
        simulate_one,
        future.seed = TRUE
      ),
      use.names = FALSE
    )
  } else {
    purrr::map_lgl(seq_len(n_sims), simulate_one)
  }

  failed_count <- sum(is.na(successes))
  valid_count <- n_sims - failed_count
  if (valid_count == 0L) {
    stop(
      "All simulated ANOVA fits failed. This usually indicates an internal ",
      "model-fitting error rather than zero power.",
      call. = FALSE
    )
  }
  success_count <- sum(successes, na.rm = TRUE)
  power <- if (valid_count > 0L) success_count / valid_count else NA_real_

  tibble::tibble(
    total_n = as.integer(sum(spec$cell_n)),
    min_cell_n = as.integer(min(spec$cell_n)),
    max_cell_n = as.integer(max(spec$cell_n)),
    n_sims = n_sims,
    num_df = sanity$num_df,
    den_df = sanity$den_df,
    ncp = round(sanity$ncp, 3),
    power_calc = round(sanity$power_calc, 3),
    power_sim = round(power, 3)
  )
}


#' @keywords internal
#' @noRd
adaptive_design_search <- function(run_one, target, n_start, n_max, tol,
                                   max_iter = 25L, progress_bar = NULL) {
  visited <- list()
  n <- max(1L, as.integer(n_start))
  lo <- NULL
  lo_p <- NA_real_
  hi <- NULL
  hi_p <- NA_real_

  repeat {
    row <- run_one(n)
    tick_progress_bar(progress_bar)
    visited[[length(visited) + 1L]] <- row
    p <- row$power_sim
    if (!is.na(p) && p >= target) {
      hi <- n
      hi_p <- p
      break
    }
    if (n >= n_max) {
      break
    }
    if (!is.na(p)) {
      lo <- n
      lo_p <- p
    }
    n <- min(n_max, max(n + 1L, n * 2L))
  }

  if (!is.null(lo) && !is.null(hi)) {
    lo_n <- lo
    hi_n <- hi
    iter <- 0L
    while (hi_n > lo_n + 1L && iter < max_iter) {
      visited_n <- vapply(
        visited,
        function(x) as.integer(x$n_per_cell[[1L]]),
        integer(1L)
      )
      next_n <- next_adaptive_n(
        lo_n = lo_n,
        lo_p = lo_p,
        hi_n = hi_n,
        hi_p = hi_p,
        target = target,
        visited_n = visited_n
      )
      if (is.na(next_n)) break

      row <- run_one(next_n)
      tick_progress_bar(progress_bar)
      visited[[length(visited) + 1L]] <- row
      p <- row$power_sim
      if (is.na(p)) break
      if (p >= target) {
        hi_n <- next_n
        hi_p <- p
      } else {
        lo_n <- next_n
        lo_p <- p
      }
      iter <- iter + 1L
    }
  }

  dplyr::bind_rows(visited) |>
    dplyr::arrange(.data$n_per_cell) |>
    dplyr::distinct(.data$n_per_cell, .keep_all = TRUE)
}


#' @keywords internal
#' @noRd
estimate_design_n_needed <- function(curve, target) {
  curve <- dplyr::arrange(curve, .data$n_per_cell)
  above <- which(curve$power_sim >= target)
  if (length(above) == 0L) return(NA_integer_)
  as.integer(curve$n_per_cell[[above[1L]]])
}


#' @keywords internal
#' @noRd
power_is_in_precision_band <- function(power_sim, target, tol) {
  is.finite(power_sim) && power_sim >= target && power_sim <= target + tol
}


#' @keywords internal
#' @noRd
next_adaptive_n <- function(lo_n, lo_p, hi_n, hi_p, target, visited_n) {
  candidates <- unsimulated_bracket_n(lo_n, hi_n, visited_n)
  if (!length(candidates)) return(NA_integer_)

  interp <- NA_integer_
  if (is.finite(lo_p) && is.finite(hi_p) && hi_p > lo_p) {
    frac <- (target - lo_p) / (hi_p - lo_p)
    raw_interp <- lo_n + frac * (hi_n - lo_n)
    interp <- as.integer(ceiling(raw_interp - sqrt(.Machine$double.eps)))
    interp <- min(max(interp, lo_n + 1L), hi_n - 1L)
  }
  if (!is.na(interp)) {
    return(as.integer(candidates[which.min(abs(candidates - interp))]))
  }

  mid <- as.integer(floor((lo_n + hi_n) / 2L))
  as.integer(candidates[which.min(abs(candidates - mid))])
}


#' @keywords internal
#' @noRd
unsimulated_bracket_n <- function(lo_n, hi_n, visited_n) {
  if (hi_n <= lo_n + 1L) return(integer())
  candidates <- seq.int(lo_n + 1L, hi_n - 1L)
  candidates[!(candidates %in% visited_n)]
}


#' @keywords internal
#' @noRd
warn_precision_band_not_reached <- function(curve, target, tol, n_needed) {
  if (is.na(n_needed)) return(invisible(NULL))
  if (any(vapply(
    curve$power_sim,
    power_is_in_precision_band,
    logical(1L),
    target = target,
    tol = tol
  ))) {
    return(invisible(NULL))
  }

  reported <- curve[curve$n_per_cell == n_needed, , drop = FALSE]
  if (!nrow(reported)) return(invisible(NULL))
  warning(
    sprintf(
      paste(
        "Requested precision band was not reached for target power %.3f",
        "with tolerance %.3f; reporting the closest explicitly simulated",
        "n_per_cell at or above target power: %d (power_sim = %.3f)."
      ),
      target,
      tol,
      as.integer(reported$n_per_cell[[1L]]),
      reported$power_sim[[1L]]
    ),
    call. = FALSE
  )
  invisible(NULL)
}


#' @keywords internal
#' @noRd
warn_power_disagreement <- function(results, n_sims, threshold = 0.05) {
  if (!all(c("power_sim", "power_calc") %in% names(results))) return(invisible(NULL))
  diff <- abs(results$power_sim - results$power_calc)
  bad <- which(is.finite(diff) & diff > threshold)
  if (!length(bad)) return(invisible(NULL))

  n_label <- if ("n_per_cell" %in% names(results)) "n_per_cell" else "total_n"
  ns <- paste(results[[n_label]][bad], collapse = ", ")
  max_diff <- max(diff[bad], na.rm = TRUE)
  msg <- paste0(
    "`power_sim` and `power_calc` differ by more than ",
    threshold * 100,
    " percentage points for ", n_label, " = ",
    ns,
    " (largest difference = ",
    sprintf("%.3f", max_diff),
    "). "
  )
  if (n_sims < 5000) {
    msg <- paste0(msg, "Try increasing `n_sims` for a more stable simulation estimate.")
  } else {
    msg <- paste0(
      msg,
      "You are already using 5000+ simulations; consider raising a GitHub issue ",
      "and treat this power estimate as potentially unstable."
    )
  }
  warning(msg, call. = FALSE)
  invisible(NULL)
}


#' @keywords internal
#' @noRd
default_term_pattern <- function(spec, term) {
  grid <- tidyr::expand_grid(!!!spec$levels[spec$factor_names])
  for (nm in spec$factor_names) {
    grid[[nm]] <- factor(grid[[nm]], levels = levels(spec$levels[[nm]]))
    stats::contrasts(grid[[nm]]) <- stats::contr.sum(nlevels(grid[[nm]]))
  }
  grid$.mu <- 0
  form <- stats::reformulate(paste(spec$factor_names, collapse = "*"),
                             response = ".mu")
  mm <- stats::model.matrix(form, data = grid)
  term_labels <- attr(stats::terms(form), "term.labels")
  idx <- match(term, term_labels)
  if (is.na(idx)) {
    stop("Term '", term, "' is not available. Available terms: ",
         paste(shQuote(term_labels), collapse = ", "), call. = FALSE)
  }
  cols <- which(attr(mm, "assign") == idx)
  pattern <- as.numeric(mm[, cols, drop = FALSE] %*% rep(1, length(cols)))
  if (all(abs(pattern) < sqrt(.Machine$double.eps))) {
    stop("The default contrast pattern for term '", term, "' is zero.",
         call. = FALSE)
  }

  between_key <- if (length(spec$between)) {
    do.call(paste, c(grid[spec$between], sep = "\r"))
  } else {
    rep(".all_between", nrow(grid))
  }
  within_key <- if (length(spec$within)) {
    do.call(paste, c(grid[spec$within], sep = "\r"))
  } else {
    rep(".all_within", nrow(grid))
  }
  between_levels <- if (length(spec$between)) {
    do.call(paste, c(spec$between_cells[spec$between], sep = "\r"))
  } else {
    ".all_between"
  }
  within_levels <- if (length(spec$within)) {
    do.call(paste, c(spec$within_cells[spec$within], sep = "\r"))
  } else {
    ".all_within"
  }

  out <- matrix(
    NA_real_,
    nrow = spec$n_between_cells,
    ncol = spec$n_within_cells
  )
  for (i in seq_along(pattern)) {
    row <- match(between_key[[i]], between_levels)
    col <- match(within_key[[i]], within_levels)
    out[row, col] <- pattern[[i]]
  }
  if (anyNA(out)) {
    stop("Could not align the default contrast pattern to the design cells.",
         call. = FALSE)
  }
  out
}


#' @keywords internal
#' @noRd
calibrate_design_means <- function(spec, term, target_pes, n, sd, r,
                                   gpower = FALSE, ss_type = "III") {
  base <- default_term_pattern(spec, term)
  exact <- simulate_balanced_design_data(
    spec = spec,
    n = n,
    means = base,
    sd = sd,
    r = r,
    empirical = TRUE
  )
  stats <- fit_design_term_stats(exact, spec, term, ss_type = ss_type)
  old_pes <- stats$pes
  if (is.na(old_pes) || old_pes <= 0 || old_pes >= 1) {
    stop("Could not calibrate the default means for term '", term, "'.",
         call. = FALSE)
  }
  calibration_pes <- calibration_pes_for_ncp(
    target_pes = target_pes,
    total_n = n * max(1L, spec$n_between_cells),
    num_df = stats$num_df,
    den_df = stats$den_df,
    gpower = gpower
  )
  k <- compute_scale_factor(old_pes, calibration_pes)
  base * k
}


#' @keywords internal
#' @noRd
calibrate_unbalanced_design_means <- function(spec, term, target_pes, sd, r,
                                              gpower = FALSE,
                                              ss_type = "III") {
  base <- default_term_pattern(spec, term)
  exact <- simulate_unbalanced_design_data(
    spec = spec,
    means = base,
    sd = sd,
    r = r,
    empirical = TRUE
  )
  stats <- fit_design_term_stats(exact, spec, term, ss_type = ss_type)
  old_pes <- stats$pes
  if (is.na(old_pes) || old_pes <= 0 || old_pes >= 1) {
    stop("Could not calibrate the default means for term '", term, "'.",
         call. = FALSE)
  }
  calibration_pes <- calibration_pes_for_ncp(
    target_pes = target_pes,
    total_n = sum(spec$cell_n),
    num_df = stats$num_df,
    den_df = stats$den_df,
    gpower = gpower
  )
  k <- compute_scale_factor(old_pes, calibration_pes)
  base * k
}


#' @keywords internal
#' @noRd
calibration_pes_for_ncp <- function(target_pes, total_n, num_df, den_df,
                                    gpower) {
  if (!isTRUE(gpower)) return(target_pes)
  f2 <- target_pes / (1 - target_pes)
  target_ncp <- total_n * f2
  calibration_f2 <- target_ncp / den_df
  calibration_pes <- calibration_f2 / (1 + calibration_f2)
  if (!is.finite(calibration_pes) || calibration_pes <= 0 ||
      calibration_pes >= 1) {
    stop("Could not calibrate G*Power-style means for this design.",
         call. = FALSE)
  }
  calibration_pes
}


#' @keywords internal
#' @noRd
sanity_check_term_effect <- function(spec, term, target_pes, n, means, sd, r,
                                     alpha, gpower, ss_type = "III") {
  exact <- simulate_balanced_design_data(
    spec = spec,
    n = n,
    means = means,
    sd = sd,
    r = r,
    empirical = TRUE
  )
  stats <- fit_design_term_stats(exact, spec, term, ss_type = ss_type)
  total_n <- n * max(1L, spec$n_between_cells)
  expected_pes <- calibration_pes_for_ncp(
    target_pes = target_pes,
    total_n = total_n,
    num_df = stats$num_df,
    den_df = stats$den_df,
    gpower = gpower
  )
  tolerance <- max(1e-6, sqrt(.Machine$double.eps) * 10)
  if (!is.finite(stats$pes) || abs(stats$pes - expected_pes) > tolerance) {
    stop(
      "Sanity check failed: empirical reference data produced partial eta ",
      "squared ", signif(stats$pes, 6), " for term '", term,
      "', expected ", signif(expected_pes, 6), ".",
      call. = FALSE
    )
  }
  ncp <- ncp_from_pes(
    pes = target_pes,
    total_n = total_n,
    den_df = stats$den_df,
    gpower = gpower
  )
  stats$ncp <- ncp
  stats$power_calc <- stats::pf(
    stats::qf(1 - alpha, stats$num_df, stats$den_df),
    stats$num_df,
    stats$den_df,
    ncp = ncp,
    lower.tail = FALSE
  )
  stats
}


#' @keywords internal
#' @noRd
sanity_check_unbalanced_term_effect <- function(spec, term, target_pes, means,
                                                sd, r, alpha, gpower,
                                                ss_type = "III") {
  exact <- simulate_unbalanced_design_data(
    spec = spec,
    means = means,
    sd = sd,
    r = r,
    empirical = TRUE
  )
  stats <- fit_design_term_stats(exact, spec, term, ss_type = ss_type)
  total_n <- sum(spec$cell_n)
  expected_pes <- calibration_pes_for_ncp(
    target_pes = target_pes,
    total_n = total_n,
    num_df = stats$num_df,
    den_df = stats$den_df,
    gpower = gpower
  )
  tolerance <- max(1e-6, sqrt(.Machine$double.eps) * 10)
  if (!is.finite(stats$pes) || abs(stats$pes - expected_pes) > tolerance) {
    stop(
      "Sanity check failed: empirical reference data produced partial eta ",
      "squared ", signif(stats$pes, 6), " for term '", term,
      "', expected ", signif(expected_pes, 6), ".",
      call. = FALSE
    )
  }
  ncp <- ncp_from_pes(
    pes = target_pes,
    total_n = total_n,
    den_df = stats$den_df,
    gpower = gpower
  )
  stats$ncp <- ncp
  stats$power_calc <- stats::pf(
    stats::qf(1 - alpha, stats$num_df, stats$den_df),
    stats$num_df,
    stats$den_df,
    ncp = ncp,
    lower.tail = FALSE
  )
  stats
}


#' @keywords internal
#' @noRd
ncp_from_pes <- function(pes, total_n, den_df, gpower) {
  f2 <- pes / (1 - pes)
  if (isTRUE(gpower)) total_n * f2 else den_df * f2
}


#' @keywords internal
#' @noRd
estimate_ncp_n_needed <- function(spec, term, target_pes, target_power, alpha,
                                  ss_type, sd, r, gpower, n_min, n_max) {
  power_at <- function(n) {
    power_calc_at_n(
      spec = spec,
      term = term,
      target_pes = target_pes,
      n = n,
      alpha = alpha,
      ss_type = ss_type,
      sd = sd,
      r = r,
      gpower = gpower
    )
  }

  p_min <- power_at(n_min)
  if (!is.finite(p_min)) return(NA_integer_)
  if (p_min >= target_power) return(as.integer(n_min))

  lo <- n_min
  hi <- min(n_max, max(n_min + 1L, n_min * 2L))
  repeat {
    p_hi <- power_at(hi)
    if (!is.finite(p_hi)) return(NA_integer_)
    if (p_hi >= target_power || hi >= n_max) break
    lo <- hi
    hi <- min(n_max, max(hi + 1L, hi * 2L))
  }
  if (hi >= n_max && power_at(hi) < target_power) return(NA_integer_)

  while (hi > lo + 1L) {
    mid <- as.integer(floor((lo + hi) / 2L))
    p_mid <- power_at(mid)
    if (!is.finite(p_mid)) return(NA_integer_)
    if (p_mid < target_power) lo <- mid else hi <- mid
  }
  as.integer(hi)
}


#' @keywords internal
#' @noRd
power_calc_at_n <- function(spec, term, target_pes, n, alpha, ss_type, sd, r,
                            gpower) {
  base <- default_term_pattern(spec, term)
  exact <- simulate_balanced_design_data(
    spec = spec,
    n = n,
    means = base,
    sd = sd,
    r = r,
    empirical = TRUE
  )
  term_stats <- fit_design_term_stats(exact, spec, term, ss_type = ss_type)
  ncp <- ncp_from_pes(
    pes = target_pes,
    total_n = n * max(1L, spec$n_between_cells),
    den_df = term_stats$den_df,
    gpower = gpower
  )
  stats::pf(
    stats::qf(1 - alpha, term_stats$num_df, term_stats$den_df),
    term_stats$num_df,
    term_stats$den_df,
    ncp = ncp,
    lower.tail = FALSE
  )
}


#' @keywords internal
#' @noRd
simulate_balanced_design_data <- function(spec, n, means, sd, r,
                                          empirical = FALSE) {
  if (isTRUE(empirical)) {
    validate_calibration_n(n, spec, "n")
  }

  sigma <- compound_symmetric_sigma(spec$n_within_cells, sd = sd, r = r)
  between_labels <- make_cell_labels("b", spec$n_between_cells)
  within_labels <- make_cell_labels("w", spec$n_within_cells)
  rownames(means) <- between_labels
  colnames(means) <- within_labels

  subject_rows <- vector("list", spec$n_between_cells)
  y_rows <- vector("list", spec$n_between_cells)

  for (i in seq_len(spec$n_between_cells)) {
    y <- MASS::mvrnorm(
      n = n,
      mu = means[i, ],
      Sigma = sigma,
      empirical = empirical
    )
    if (spec$n_within_cells == 1L) {
      y <- matrix(y, nrow = n, ncol = 1L)
    } else if (n == 1L) {
      y <- matrix(y, nrow = 1L)
    }
    colnames(y) <- within_labels

    b <- spec$between_cells[i, spec$between, drop = FALSE]
    subject_rows[[i]] <- dplyr::bind_cols(
      tibble::tibble(id = seq_len(n) + (i - 1L) * n),
      b[rep(1L, n), , drop = FALSE]
    )
    y_rows[[i]] <- y
  }

  subjects <- dplyr::bind_rows(subject_rows)
  y_mat <- do.call(rbind, y_rows)
  wide <- dplyr::bind_cols(subjects, tibble::as_tibble(y_mat))

  if (!length(spec$within)) {
    out <- wide |>
      dplyr::rename(value = dplyr::all_of(within_labels[1L]))
  } else {
    within_map <- dplyr::bind_cols(
      tibble::tibble(.within_cell = within_labels),
      spec$within_cells[, spec$within, drop = FALSE]
    )
    out <- wide |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(within_labels),
        names_to = ".within_cell",
        values_to = "value"
      ) |>
      dplyr::left_join(within_map, by = ".within_cell") |>
      dplyr::select(-".within_cell")
  }

  out$id <- factor(out$id)
  for (nm in spec$factor_names) {
    out[[nm]] <- factor(out[[nm]], levels = levels(spec$levels[[nm]]))
  }
  tibble::as_tibble(out)
}


#' @keywords internal
#' @noRd
simulate_unbalanced_design_data <- function(spec, means, sd, r,
                                            empirical = FALSE) {
  if (isTRUE(empirical)) {
    validate_calibration_n(spec$cell_n, spec, "cell_n$n")
  }

  sigma <- compound_symmetric_sigma(spec$n_within_cells, sd = sd, r = r)
  between_labels <- make_cell_labels("b", spec$n_between_cells)
  within_labels <- make_cell_labels("w", spec$n_within_cells)
  rownames(means) <- between_labels
  colnames(means) <- within_labels

  subject_rows <- vector("list", spec$n_between_cells)
  y_rows <- vector("list", spec$n_between_cells)
  id_offset <- 0L

  for (i in seq_len(spec$n_between_cells)) {
    n_i <- spec$cell_n[[i]]
    y <- MASS::mvrnorm(
      n = n_i,
      mu = means[i, ],
      Sigma = sigma,
      empirical = empirical
    )
    if (spec$n_within_cells == 1L) {
      y <- matrix(y, nrow = n_i, ncol = 1L)
    } else if (n_i == 1L) {
      y <- matrix(y, nrow = 1L)
    }
    colnames(y) <- within_labels

    b <- spec$between_cells[i, spec$between, drop = FALSE]
    subject_rows[[i]] <- dplyr::bind_cols(
      tibble::tibble(id = seq_len(n_i) + id_offset),
      b[rep(1L, n_i), , drop = FALSE]
    )
    y_rows[[i]] <- y
    id_offset <- id_offset + n_i
  }

  subjects <- dplyr::bind_rows(subject_rows)
  y_mat <- do.call(rbind, y_rows)
  wide <- dplyr::bind_cols(subjects, tibble::as_tibble(y_mat))

  if (!length(spec$within)) {
    out <- wide |>
      dplyr::rename(value = dplyr::all_of(within_labels[1L]))
  } else {
    within_map <- dplyr::bind_cols(
      tibble::tibble(.within_cell = within_labels),
      spec$within_cells[, spec$within, drop = FALSE]
    )
    out <- wide |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(within_labels),
        names_to = ".within_cell",
        values_to = "value"
      ) |>
      dplyr::left_join(within_map, by = ".within_cell") |>
      dplyr::select(-".within_cell")
  }

  out$id <- factor(out$id)
  for (nm in spec$factor_names) {
    out[[nm]] <- factor(out[[nm]], levels = levels(spec$levels[[nm]]))
  }
  tibble::as_tibble(out)
}


#' @keywords internal
#' @noRd
minimum_calibration_n <- function(spec) {
  if (spec$n_within_cells > 1L) spec$n_within_cells + 1L else 2L
}


#' @keywords internal
#' @noRd
validate_calibration_n <- function(n, spec, arg) {
  min_n <- minimum_calibration_n(spec)
  too_small <- n < min_n
  if (!any(too_small)) return(invisible(n))

  bad <- paste(sort(unique(n[too_small])), collapse = ", ")
  stop(
    "`", arg, "` is too small for this design. This design has ",
    spec$n_within_cells, " within-subject cell",
    if (spec$n_within_cells == 1L) "" else "s",
    ", so empirical calibration needs at least ", min_n,
    " participant", if (min_n == 1L) "" else "s",
    " per between-subject cell. ",
    "Increase `", arg, "`",
    if (length(n) > 1L) " values" else "",
    " (too small: ", bad, ").",
    call. = FALSE
  )
}


#' @keywords internal
#' @noRd
compound_symmetric_sigma <- function(n, sd, r) {
  sigma <- matrix(r * sd^2, nrow = n, ncol = n)
  diag(sigma) <- sd^2
  sigma
}


#' @keywords internal
#' @noRd
fit_design_model <- function(data, spec) {
  data <- set_sum_contrasts(data, spec)
  formula <- design_aov_formula(spec)
  suppressWarnings(stats::aov(formula, data = data))
}


#' @keywords internal
#' @noRd
fit_design_term <- function(data, spec, term, alpha, ss_type = "III") {
  stats <- tryCatch(
    fit_design_term_stats(data, spec, term, ss_type = ss_type),
    error = function(e) NULL
  )
  if (is.null(stats) || !is.finite(stats$p_value)) return(NA)
  isTRUE(stats$p_value < alpha)
}


#' @keywords internal
#' @noRd
fit_design_term_stats <- function(data, spec, term, ss_type = "III") {
  ss_type <- validate_ss_type(ss_type)
  if (identical(ss_type, "I")) {
    fit <- fit_design_model(data, spec)
    return(extract_term_stats(fit, term))
  }
  fit_car_term_stats(data, spec, term, ss_type = ss_type)
}


#' @keywords internal
#' @noRd
fit_car_term_stats <- function(data, spec, term, ss_type) {
  data <- set_sum_contrasts(data, spec)
  type <- switch(ss_type, II = 2, III = 3)
  if (!length(spec$within)) {
    fixed <- paste(spec$factor_names, collapse = " * ")
    formula <- stats::as.formula(paste("value ~", fixed))
    fit <- stats::lm(formula, data = data)
    tab <- as.data.frame(car::Anova(fit, type = type))
    rows <- car_between_rows(tab)
  } else {
    wide_setup <- repeated_measures_wide_data(data, spec)
    fit <- stats::lm(wide_setup$formula, data = wide_setup$wide)
    idesign <- stats::reformulate(paste(spec$within, collapse = " * "))
    av <- car::Anova(
      fit,
      idata = wide_setup$idata,
      idesign = idesign,
      type = type,
      icontrasts = c("contr.sum", "contr.poly")
    )
    tab <- as.data.frame.matrix(
      suppressWarnings(summary(av, multivariate = FALSE))$univariate.tests
    )
    rows <- car_repeated_rows(tab)
  }

  idx <- which(rows$term == term)
  if (!length(idx)) {
    stop("Term '", term, "' was not found in the fitted ANOVA table.",
         call. = FALSE)
  }
  row <- rows[idx[[1L]], , drop = FALSE]
  if (!is.finite(row$f_value) || !is.finite(row$den_df)) {
    stop("Term '", term, "' does not have a finite F test in the fitted ANOVA.",
         call. = FALSE)
  }
  pes <- row$f_value * row$num_df / (row$f_value * row$num_df + row$den_df)
  list(
    num_df = as.numeric(row$num_df),
    den_df = as.numeric(row$den_df),
    pes = as.numeric(pes),
    p_value = as.numeric(row$p_value)
  )
}


#' @keywords internal
#' @noRd
car_between_rows <- function(tab) {
  row_terms <- trimws(rownames(tab))
  residual_idx <- which(row_terms == "Residuals")
  den_df <- if (length(residual_idx)) {
    as.numeric(tab[residual_idx[[1L]], "Df"])
  } else {
    NA_real_
  }
  keep <- row_terms != "Residuals" & row_terms != "(Intercept)" & nzchar(row_terms)
  tibble::tibble(
    term = row_terms[keep],
    num_df = as.numeric(tab[keep, "Df"]),
    den_df = den_df,
    f_value = as.numeric(tab[keep, "F value"]),
    p_value = as.numeric(tab[keep, "Pr(>F)"])
  )
}


#' @keywords internal
#' @noRd
car_repeated_rows <- function(tab) {
  row_terms <- trimws(rownames(tab))
  keep <- row_terms != "(Intercept)" & nzchar(row_terms)
  tibble::tibble(
    term = row_terms[keep],
    num_df = as.numeric(tab[keep, "num Df"]),
    den_df = as.numeric(tab[keep, "den Df"]),
    f_value = as.numeric(tab[keep, "F value"]),
    p_value = as.numeric(tab[keep, "Pr(>F)"])
  )
}


#' @keywords internal
#' @noRd
repeated_measures_wide_data <- function(data, spec) {
  within_labels <- make_cell_labels("w", spec$n_within_cells)
  within_keys <- interaction_key(spec$within_cells, spec$within)
  data$.within_key <- interaction_key(data, spec$within)
  data$.within_cell <- within_labels[match(data$.within_key, within_keys)]
  data$.within_key <- NULL

  id_cols <- c("id", spec$between)
  wide <- data |>
    dplyr::select(dplyr::all_of(c(id_cols, ".within_cell", "value"))) |>
    tidyr::pivot_wider(
      names_from = ".within_cell",
      values_from = "value"
    )
  for (nm in spec$between) {
    wide[[nm]] <- factor(wide[[nm]], levels = levels(spec$levels[[nm]]))
    stats::contrasts(wide[[nm]]) <- stats::contr.sum(nlevels(wide[[nm]]))
  }

  response <- paste0("cbind(", paste(within_labels, collapse = ", "), ")")
  rhs <- if (length(spec$between)) {
    paste(spec$between, collapse = " * ")
  } else {
    "1"
  }
  idata <- spec$within_cells[, spec$within, drop = FALSE]
  for (nm in spec$within) {
    idata[[nm]] <- factor(idata[[nm]], levels = levels(spec$levels[[nm]]))
    stats::contrasts(idata[[nm]]) <- stats::contr.sum(nlevels(idata[[nm]]))
  }

  list(
    wide = wide,
    idata = as.data.frame(idata),
    formula = stats::as.formula(paste(response, "~", rhs))
  )
}


#' @keywords internal
#' @noRd
extract_term_stats <- function(fit, term) {
  rows <- extract_aov_rows(fit)
  idx <- which(rows$term == term)
  if (!length(idx)) {
    stop("Term '", term, "' was not found in the fitted ANOVA table.",
         call. = FALSE)
  }
  row <- rows[idx[[1L]], , drop = FALSE]
  if (!is.finite(row$f_value) || !is.finite(row$den_df)) {
    stop("Term '", term, "' does not have a finite F test in the fitted ANOVA.",
         call. = FALSE)
  }
  pes <- row$f_value * row$num_df / (row$f_value * row$num_df + row$den_df)
  list(
    num_df = as.numeric(row$num_df),
    den_df = as.numeric(row$den_df),
    pes = as.numeric(pes),
    p_value = as.numeric(row$p_value)
  )
}


#' @keywords internal
#' @noRd
set_sum_contrasts <- function(data, spec) {
  for (nm in spec$factor_names) {
    if (nm %in% names(data)) {
      data[[nm]] <- factor(data[[nm]], levels = levels(spec$levels[[nm]]))
      stats::contrasts(data[[nm]]) <- stats::contr.sum(nlevels(data[[nm]]))
    }
  }
  data$id <- factor(data$id)
  data
}


#' @keywords internal
#' @noRd
design_aov_formula <- function(spec) {
  fixed <- paste(spec$factor_names, collapse = " * ")
  if (!length(spec$within)) {
    return(stats::as.formula(paste("value ~", fixed)))
  }
  within <- paste(spec$within, collapse = " * ")
  stats::as.formula(
    paste("value ~", fixed, "+ Error(id / (", within, "))")
  )
}


#' @keywords internal
#' @noRd
extract_aov_rows <- function(fit) {
  s <- summary(fit)
  if (inherits(fit, "aovlist")) {
    tables <- lapply(s, `[[`, 1L)
  } else {
    tables <- list(s[[1L]])
  }
  rows <- vector("list", length(tables))
  for (i in seq_along(tables)) {
    tab <- as.data.frame(tables[[i]])
    row_terms <- trimws(rownames(tab))
    residual_idx <- which(row_terms == "Residuals")
    den_df <- if (length(residual_idx)) {
      as.numeric(tab[residual_idx[[1L]], "Df"])
    } else {
      NA_real_
    }
    keep <- row_terms != "Residuals" & nzchar(row_terms)
    rows[[i]] <- tibble::tibble(
      term = row_terms[keep],
      num_df = as.numeric(tab[keep, "Df"]),
      den_df = den_df,
      f_value = as.numeric(tab[keep, "F value"]),
      p_value = as.numeric(tab[keep, "Pr(>F)"])
    )
  }
  dplyr::bind_rows(rows)
}


#' @keywords internal
#' @noRd
make_progress_bar <- function(enabled, total, label) {
  if (!isTRUE(enabled)) return(NULL)
  total <- max(1L, as.integer(total))
  bar <- utils::txtProgressBar(min = 0, max = total, style = 3)
  env <- new.env(parent = emptyenv())
  env$bar <- bar
  env$total <- total
  env$value <- 0L
  env$label <- label
  env
}


#' @keywords internal
#' @noRd
tick_progress_bar <- function(progress_bar, by = 1L) {
  if (is.null(progress_bar)) return(invisible(NULL))
  progress_bar$value <- min(progress_bar$total, progress_bar$value + by)
  utils::setTxtProgressBar(progress_bar$bar, progress_bar$value)
  invisible(NULL)
}


#' @keywords internal
#' @noRd
close_progress_bar <- function(progress_bar) {
  if (is.null(progress_bar)) return(invisible(NULL))
  utils::setTxtProgressBar(progress_bar$bar, progress_bar$total)
  close(progress_bar$bar)
  invisible(NULL)
}


#' @keywords internal
#' @noRd
estimate_adaptive_progress_total <- function(n_start, n_max, max_iter = 25L) {
  n <- max(1L, as.integer(n_start))
  n_max <- as.integer(n_max)
  steps <- 1L
  while (n < n_max) {
    n_next <- min(n_max, max(n + 1L, n * 2L))
    steps <- steps + 1L
    n <- n_next
  }
  steps + max_iter
}
