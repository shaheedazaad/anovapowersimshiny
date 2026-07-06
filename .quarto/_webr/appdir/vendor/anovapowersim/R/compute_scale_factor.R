#' Compute the mean-deviation scaling factor from a change in partial eta squared
#'
#' Given an existing partial eta squared for a term and a target partial eta
#' squared, returns the multiplier `k` that must be applied to that term's
#' additive contribution to the cell means in order to obtain the target effect
#' size under the same residual structure.
#'
#' The derivation is straightforward: partial eta squared can be written as
#' `pes = SS_effect / (SS_effect + C)`, where `C` is the part of the
#' denominator held fixed by this package's rescaling. Thus
#' `pes / (1 - pes)` scales as the target effect's sum of squares. Scaling the
#' term's deviations by `k` scales the target effect's sum of squares by
#' `k^2`, so the required multiplier is
#' \deqn{k = \sqrt{\frac{p_{\mathrm{new}} / (1 - p_{\mathrm{new}})}
#'                       {p_{\mathrm{old}} / (1 - p_{\mathrm{old}})}}.}
#'
#' @param old_pes Numeric scalar in (0, 1), or a numeric-looking character
#'   scalar such as `".310"`. The current partial eta squared for the term of
#'   interest.
#' @param new_pes Numeric scalar in (0, 1), or a numeric-looking character
#'   scalar such as `".200"`. The target partial eta squared.
#'
#' @return A single positive numeric value `k`. `k > 1` amplifies the
#'   effect, `k < 1` shrinks it, and `k == 1` leaves it unchanged.
#'
#' @examples
#' compute_scale_factor(0.10, 0.05)   # shrink
#' compute_scale_factor(0.05, 0.10)   # amplify
#'
#' @seealso [design_term_means()], [power_curve()]
#' @export
compute_scale_factor <- function(old_pes, new_pes) {
  old_pes <- as_unit_interval(old_pes, "old_pes")
  new_pes <- as_unit_interval(new_pes, "new_pes")

  old_f2 <- old_pes / (1 - old_pes)
  new_f2 <- new_pes / (1 - new_pes)

  sqrt(new_f2 / old_f2)
}

