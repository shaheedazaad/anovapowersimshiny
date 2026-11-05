#' Convert Cohen's f to partial eta squared
#'
#' Converts Cohen's f for an effect to its corresponding partial eta squared:
#' \deqn{\eta_p^2 = \frac{f^2}{1 + f^2}.}
#' The same formula applies to between-subject, within-subject, and interaction
#' effects when f uses the error variance corresponding to that effect's
#' partial eta squared. This helper does not translate between G*Power's
#' repeated-measures effect-size conventions.
#'
#' @param f A single finite, nonnegative numeric value representing Cohen's f,
#'   not an ANOVA F statistic or f-squared.
#'
#' @return A single numeric partial eta squared. Zero maps to zero; power
#'   functions retain their own restrictions on `target_pes`. Very large f
#'   values may yield exactly one because of floating-point rounding.
#'
#' @references Cohen, J. (1988). Statistical power analysis for the behavioral
#'   sciences (2nd ed.). Lawrence Erlbaum Associates.
#'
#' @examples
#' f_to_pes(0.25) # approximately 0.05882
#' power_n_calc(
#'   between = c(group = 2),
#'   term = "group",
#'   target_pes = f_to_pes(0.25)
#' )
#'
#' @seealso [compute_scale_factor()], [power_n_calc()]
#' @export
f_to_pes <- function(f) {
  if (!is.numeric(f) || length(f) != 1L || !is.finite(f) || f < 0) {
    stop("`f` must be a single finite, nonnegative numeric value.",
         call. = FALSE)
  }

  f <- as.numeric(f)
  if (f > 1) {
    return(1 / (1 + (1 / f)^2))
  }
  f^2 / (1 + f^2)
}
