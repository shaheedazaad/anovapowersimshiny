#' @keywords internal
#' @noRd
balanced_anova_design <- function(between = NULL, within = NULL) {
  spec <- validate_design_spec(between, within)
  class(spec) <- c("anovapowersim_design_spec", class(spec))
  spec
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
    expand_factor_grid(levels[names(between)])
  } else {
    data.frame(.dummy_between = factor("all"), check.names = FALSE)
  }
  within_cells <- if (length(within)) {
    expand_factor_grid(levels[names(within)])
  } else {
    data.frame(.dummy_within = factor("dv"), check.names = FALSE)
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
expand_factor_grid <- function(factors) {
  out <- do.call(
    expand.grid,
    c(factors, list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
  )
  for (nm in names(factors)) {
    out[[nm]] <- factor(out[[nm]], levels = levels(factors[[nm]]))
  }
  out
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
ncp_from_pes <- function(pes, total_n, den_df, gpower) {
  f2 <- pes / (1 - pes)
  if (isTRUE(gpower)) total_n * f2 else den_df * f2
}
