#' @title Find Outlier Groups Based on Energy Distance
#' @description Identifies groups (e.g., studies) that are most distant from the average
#' group based on energy distance across multiple variables.
#' @param formula A formula specifying the group variable and variables.
#' e.g., `study ~ var1 + var2 +...`. The group variable should be a factor or will be converted to one.
#' @param data A data frame containing the variables specified in the formula.
#' @param cutoff Numeric. Percentile threshold for the permutation-based cutoff (default 0.99).
#' The cutoff is determined by permuting group labels and calculating the percentile of
#' permuted median distances.
#' @param R Integer. Number of permutations for determining the cutoff (default 500).
#' @param plot Logical. If TRUE (default), returns a visualization of the outlier analysis.
#' @return If `plot = TRUE`, returns a list with:
#' \itemize{
#'   \item `cutoff_value`: The permutation-based cutoff value used for outlier detection.
#'   \item `summary`: Data frame with group, median_distance, outlier_score, and is_outlier columns.
#'   \item `heatmap`: A ggplot2 heatmap of pairwise energy distances.
#'   \item `barplot`: A ggplot2 bar plot showing median distance to other groups.
#' }
#' If `plot = FALSE`, returns only the elements without plots.
#' @details
#' Groups with high median distance to other groups are identified as potential outliers.
#' The outlier_score is a z-score that indicates how many standard deviations a group's
#' median distance is from the overall median distance.
#'
#' Before distance calculation, all covariates are scaled to mean 0 and standard deviation 1.
#' @examples
#'
#' # Example 1: 10 studies with real outliers (Study-8, Study-9, Study-10)
#' set.seed(123)
#' dat <- data.frame(
#'   study = factor(rep(paste0("Study-", 1:10), each = 20)),
#'   var1 = c(rnorm(20, 10, 1), rnorm(20, 10, 1), rnorm(20, 10, 1), rnorm(20, 10, 1),
#'            rnorm(20, 10, 1), rnorm(20, 10, 1), rnorm(20, 10, 1), rnorm(20, 15, 1),
#'            rnorm(20, 10, 1), rnorm(20, 16, 1)),
#'   var2 = c(rnorm(20, 5, 1), rnorm(20, 5, 1), rnorm(20, 5, 1), rnorm(20, 5, 1),
#'            rnorm(20, 5, 1), rnorm(20, 5, 1), rnorm(20, 5, 1), rnorm(20, 5, 1),
#'            rnorm(20, 10, 1), rnorm(20, 5, 1))
#' )
#' out <- find_outliers(study ~ var1 + var2, data = dat, R = 200)
#' out$summary      # Study-8, Study-9, Study-10 should be flagged
#' out$cutoff_value # Permutation-based threshold
#'
#' # Example 2: 20 studies with NO real outliers (all from same distribution)
#' set.seed(456)
#' dat_no_outliers <- data.frame(
#'   study = factor(rep(paste0("Study-", 1:20), each = 15)),
#'   var1 = rnorm(300, 10, 2),
#'   var2 = rnorm(300, 5, 1)
#' )
#' out2 <- find_outliers(study ~ var1 + var2, data = dat_no_outliers, R = 200)
#' out2$summary     # Should have few or no outliers flagged
#' sum(out2$is_outlier)  # Count of flagged outliers (expected: 0 or very few)
#'
#' @rdname find_outliers
#' @export
#' @importFrom stats dist sd na.omit quantile
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient2 theme_minimal labs coord_fixed geom_bar geom_hline geom_vline element_text




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
find_outliers <- function(formula, data, cutoff = 0.99, R = 500, plot = TRUE) {
  Group1 <- Group2 <- Distance <- Group <- Median_Distance <- NULL

  # Parse formula
  if (any(grepl('|', as.character(formula), fixed = TRUE))) {
    stop("The find_outliers function does not support stratum notation. Use only: group ~ var1 + var2 + ...")
  }

  fnames <- all.vars(as.formula(formula))
  group_var <- fnames[1]
  vars <- fnames[-1]

  # Validate inputs
  if (length(vars) < 1) stop("At least one covariate must be specified in the formula.")

  # Prepare data
  data <- data[, c(group_var, vars), drop = FALSE]
  if (any(is.na(data))) warning("NA values were removed from the data frame.")
  data <- na.omit(data)

  # Ensure group is a factor
  data[[group_var]] <- as.factor(data[[group_var]])
  groups <- levels(data[[group_var]])
  n_groups <- length(groups)

  if (n_groups < 2) stop("At least two groups are required for comparison.")
  if (n_groups < 3) warning("With only 2 groups, outlier detection is not meaningful.")

  # Scale covariates (mean=0, sd=1)
  data[, vars] <- scale(data[, vars])

  # Compute full Euclidean distance matrix ONCE (efficient approach)
  all_covs <- as.matrix(data[, vars, drop = FALSE])
  full_dist <- as.matrix(dist(all_covs, method = "euclidean"))

  # Get indices for each group
  group_indices <- lapply(groups, function(g) which(data[[group_var]] == g))
  names(group_indices) <- groups

  # Internal function: compute energy distance from precomputed distance matrix
  energy_dist_from_matrix <- function(dist_mat, idx1, idx2) {
    n <- length(idx1)
    m <- length(idx2)

    D_XY <- dist_mat[idx1, idx2, drop = FALSE]
    D_XX <- dist_mat[idx1, idx1, drop = FALSE]
    D_YY <- dist_mat[idx2, idx2, drop = FALSE]

    A <- 2 * mean(D_XY)
    B <- if (n > 1) mean(D_XX[lower.tri(D_XX)]) / (n / (n - 1)) else 0
    C <- if (m > 1) mean(D_YY[lower.tri(D_YY)]) / (m / (m - 1)) else 0

    ED <- A - B - C
    if (A > 0) ED <- ED / A
    return(ED)
  }

  # Compute pairwise energy distances efficiently (only upper triangle, then mirror)
  dist_matrix <- matrix(0, nrow = n_groups, ncol = n_groups)
  rownames(dist_matrix) <- colnames(dist_matrix) <- groups

  for (i in 1:(n_groups - 1)) {
    for (j in (i + 1):n_groups) {
      ed <- energy_dist_from_matrix(full_dist, group_indices[[i]], group_indices[[j]])
      dist_matrix[i, j] <- ed
      dist_matrix[j, i] <- ed
    }
  }


  # Calculate median distance for each group to all others (robust to outlier comparisons)
  median_distances <- apply(dist_matrix, 1, function(x) median(x[x > 0]))

  # Validate parameters
  if (cutoff < 0 || cutoff > 1) stop("cutoff must be between 0 and 1.")
  if (R < 100) warning("R should be at least 100 for reliable permutation-based cutoff.")

  # Ranking (most outlier = highest median distance first)
  ranking <- names(sort(median_distances, decreasing = TRUE))

  # Permutation-based cutoff: permute group labels and calculate null distribution
  # Collect ALL median distances from all permutations to build null distribution
  perm_all_distances <- numeric(R * n_groups)
  original_groups <- data[[group_var]]

  for (r in 1:R) {
    # Permute group labels
    perm_groups <- sample(original_groups)

    # Get indices for each permuted group
    perm_indices <- lapply(groups, function(g) which(perm_groups == g))

    # Calculate pairwise distances for permuted groups
    perm_dist_matrix <- matrix(0, nrow = n_groups, ncol = n_groups)
    for (i in 1:(n_groups - 1)) {
      for (j in (i + 1):n_groups) {
        ed <- energy_dist_from_matrix(full_dist, perm_indices[[i]], perm_indices[[j]])
        perm_dist_matrix[i, j] <- ed
        perm_dist_matrix[j, i] <- ed
      }
    }
    # Store ALL median distances from this permutation
    perm_all_distances[((r - 1) * n_groups + 1):(r * n_groups)] <- apply(perm_dist_matrix, 1, function(x) median(x[x > 0]))
  }

  # Calculate cutoff and median from permutation distribution
  cutoff_value <- quantile(perm_all_distances, probs = cutoff)
  null_median <- median(perm_all_distances)
  is_outlier <- median_distances > cutoff_value

  # Create summary data frame
  summary_df <- data.frame(
    group = names(median_distances),
    median_distance = median_distances,
    is_outlier = is_outlier,
    row.names = NULL
  )
  summary_df <- summary_df[order(summary_df$median_distance, decreasing = TRUE), ]

  result <- list(
    cutoff_value = as.numeric(cutoff_value),
    null_median = as.numeric(null_median),
    summary = summary_df
  )

  if (plot) {
    # Create heatmap of pairwise distances
    dist_df <- expand.grid(Group1 = factor(groups, levels = groups),
                           Group2 = factor(groups, levels = rev(groups)))
    dist_df$Distance <- sapply(1:nrow(dist_df), function(i) {
      dist_matrix[as.character(dist_df$Group1[i]), as.character(dist_df$Group2[i])]
    })
    dist_df$Distance <- round(dist_df$Distance, 2)

    p_heatmap <- ggplot(dist_df, aes(x = Group1, y = Group2, fill = Distance)) +
      geom_tile(color = "white") +
      geom_text(aes(label = round(Distance, 3)), size = 3) +
      scale_fill_gradient2(low = "#00A091", mid = "#00A091", high = "#FE7500",
                           midpoint = null_median) +
      theme_minimal() +
      coord_fixed() +
      labs(title = "Pairwise Energy Distances", x = "", y = "", fill = "E-Distance") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    # Create bar plot of median distances (sorted by outlier rank, vertical layout)
    bar_df <- data.frame(
      Group = factor(ranking, levels = rev(ranking)),
      Median_Distance = median_distances[ranking]
    )

    p_barplot <- ggplot(bar_df, aes(x = Median_Distance, y = Group, fill = Median_Distance)) +
      geom_bar(stat = "identity", color = "white") +
      geom_vline(xintercept = cutoff_value, linetype = "dashed",
                 color = "#FE7500", linewidth = 1) +
      scale_fill_gradient2(low = "#00A091", mid = "#00A091", high = "#FE7500", midpoint = cutoff_value) +
      theme_minimal() +
      labs(title = "Median Distance to Other Groups",
           subtitle = paste0("Dashed line = ", cutoff * 100, "th percentile"),
           x = "Median Energy Distance", y = "", fill = "E-Distance")

    result$heatmap <- p_heatmap
    result$barplot <- p_barplot
  }

  return(result)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

