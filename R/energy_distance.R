#' @title Compute Energy Distance Between Two Groups
#' @description Calculates the energy distance between two groups.
#' @param formula A formula specifying the treated and covariates, e.g., `treated ~ cov1 + cov2`. The treated variable must be binary (0=pool, 1=treated)
#' @param data A data frame containing the variables specified in the formula.
#' @param standardized If TRUE, the standardized energy distance that lies in the range 0 to 1 is returned, the so-called E-coefficient.
#' If FALSE, non-scaled energy distance is returned that can be >1.
#' @return A numeric value representing the energy distance between the two groups.
#' @details Energy distance is a non-parametric measure of distributional difference. It is sensitive to differences
#' in location, scale, and shape between groups. Before calculation, the covariates are scaled to a mean value of 0 and a standard deviation of 1.
#' @examples
#' dat <- data.frame(
#'  treated = rep(0:1, c(50, 30)),
#'  age    = c(rnorm(50, 5, 2),   rnorm(30, 5, 1)),
#'  weight = c(rnorm(50, 11, 2),  rnorm(30, 10, 1)),
#'  class  = c(rbinom(50, 3, 0.6),   rbinom(30, 3, 0.4))
#'  )
#'
#'  energy_distance(treated ~ age + weight + class, data=dat)
#'
#' @rdname energy_distance
#' @export




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Function to compute energy distance
energy_distance <- function(formula, data, standardized=TRUE){

  if(any(grepl('|', as.character(formula), fixed = TRUE))) stop('The energy_distance function does not support stratum analysis, but your formula contains | ')
  fnames   <- all.vars(as.formula(formula))
  leftpart <- fnames[1]   # The response variable
  vars     <- fnames[-1]  # Exclude the response variable

      data <- data[,c(leftpart,vars)]
      if(any(is.na(data))) warning('NA values were removed from the data frame.')
      data <- na.omit(data[,c(leftpart,vars)])

      data[,vars] <- scale(data[,vars])
      lvl         <- levels(factor(as.factor(data[,leftpart])))
      X           <- data[which(data[,leftpart]==lvl[1]),vars]
      Y           <- data[which(data[,leftpart]==lvl[2]),vars]
      X <- as.matrix(X)
      Y <- as.matrix(Y)
      n <- nrow(X)
      m <- nrow(Y)
      dist_XY <- as.matrix(stats::dist(rbind(X, Y)))
      D_XY <- dist_XY[1:n, (n+1):(n+m)]
      D_XX <- dist_XY[1:n, 1:n]
      D_YY <- dist_XY[(n+1):(n+m), (n+1):(n+m)]
      A <- 2 * mean(D_XY)
      B <- mean(D_XX[lower.tri(D_XX)]) / (n / (n - 1))
      C <- mean(D_YY[lower.tri(D_YY)]) / (m / (m - 1))
      ED <- A-B-C
      if(standardized){ED <- if(A > 0) ED/A else 0 }
      return(ED)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#




