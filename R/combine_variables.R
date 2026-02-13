#' @title Compute Weighted Combined Score from Multiple ariables
#' @description Creates a single combined variable from multiple ariables,
#' weighted by relative importance scores. Variables are scaled before combination.
#' @param formula A formula specifying the ariables, e.g., `1 ~ cov1 + cov2 + cov3`.
#'        The left side should be 1 (placeholder).
#' @param data A data frame containing the variables specified in the formula.
#' @param weights A numeric vector of relative importance weights for each variable.
#'        Must have the same length as the number of variables in the formula.
#'        Higher weights indicate greater contribution (correlation) to the combined score.
#'        If NULL, equal weights (all 1s) are used. Default: NULL.
#' @return A numeric vector of the same length as nrow(data), representing the combined
#'         weighted score for each observation. The score correlates with each input
#'         variable proportionally to its importance weight.
#' @details Variables are first scaled to mean 0 and standard deviation 1, then
#'          multiplied by their importance weights and summed. The final score is
#'          a weighted linear combination: score = w1*z1 + w2*z2 + ... + wk*zk,
#'          where zi are the scaled variables and wi are the weights.
#' @examples
#' dat <- data.frame(
#'   age    = rnorm(80, 5, 2),
#'   weight = rnorm(80, 11, 2),
#'   class  = rbinom(80, 3, 0.5)
#' )
#'
#' # Equal weights
#' score <- combine_variables(1 ~ age + weight + class, data = dat)
#'
#' # Custom weights: age contributes 2x, weight 1.5x, class 1x
#' score <- combine_variables(1 ~ age + weight + class, data = dat,
#'                            weights = c(2, 1.5, 1))
#'
#' @rdname combine_variables
#' @export



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Function to compute combined distance score with variable importance weights
combine_variables <- function(formula, data, weights = NULL) {

  # Parse formula - with "1 ~" syntax, all.vars returns only the RHS variables
  fnames <- all.vars(as.formula(formula))
  vars   <- fnames

  # Check for unsupported formula syntax
  if (any(grepl('|', as.character(formula), fixed = TRUE))) {
    stop('The combined_distance function does not support stratum analysis (|).')
  }
  if (any(grepl('*', as.character(formula), fixed = TRUE))) {
    stop('Only + for covariates is allowed.')
  }
  if (any(grepl(':', as.character(formula), fixed = TRUE))) {
    stop('Only + for covariates is allowed.')
  }

  # Validate that variables exist in data
  missing_vars <- vars[!vars %in% colnames(data)]
  if (length(missing_vars) > 0) {
    stop(paste0("Variables not found in data: ", paste(missing_vars, collapse = ", ")))
  }

  # Subset data to relevant columns
  data_subset <- data[, vars, drop = FALSE]

  # Handle missing values - track which rows have NA
  na_rows <- apply(data_subset, 1, function(x) any(is.na(x)))
  if (any(na_rows)) {
    warning('NA values present in the data. Rows with NA will have NA in the result.')
  }

  # Scale variables (mean=0, sd=1)
  data_scaled <- scale(data_subset)

  # Set default weights if not provided
  if (is.null(weights)) {
    weights <- rep(1, length(vars))
  }

  # Validate weights
  if (length(weights) != length(vars)) {
    stop(paste0("Length of weights (", length(weights),
                ") must match number of variables (", length(vars), ")"))
  }
  if (any(weights <= 0)) {
    stop('Weights must be positive. Variables with zero importance should be removed from the formula.')
  }
  if (any(!is.finite(weights))) {
    stop('Weights must be finite numeric values.')
  }

  # Apply weights to scaled data (same approach as in VCG_sampler)
  data_weighted <- sweep(data_scaled, 2, weights, `*`)

  # Compute weighted sum (linear combination) for each row
  # This preserves the directional relationship with input variables
  # Score correlates with each variable proportionally to its weight
  result <- rowSums(data_weighted, na.rm=TRUE)

  return(result)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

