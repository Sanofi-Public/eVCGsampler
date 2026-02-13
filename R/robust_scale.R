#' @title Robust Scaling of Numeric and Categorical Variables
#' @description Applies robust scaling to numeric and categorical variables. For numeric variables,
#' the function centers by the median and scales by the MAD. For categorical variables
#' with 2–4 unique levels, it applies a custom transformation to map them to numeric values.
#' @param x A numeric vector, factor, matrix, or data frame. If a matrix or data frame is provided, scaling is applied
#'          column-wise.
#'
#' @param group vector indicating which group is the TG to scale to
#' @return A scaled numeric vector or a data frame with scaled columns.
#' @details This function is designed to make numeric and categorical variables comparable.
#' This is an internal function that should not be used by package users.
#' @examples
#'
#' dat<-data.frame(x=rnorm(100, 10, 3), sex=factor(rbinom(100, 1, 0.5), labels=c("M","F")))
#'
#' x<- robust_scale(dat$x, dat$sex)
#' round(median(x), 2)
#' round(mad(x), 2)
#'
#' @rdname robust_scale
#' @importFrom stats median mad sd
#' @export


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
robust_scale <- function(x, group) {
  #if (!is.numeric(x)) stop("Input must be numeric.")
  if (is.matrix(x) || is.data.frame(x)) {
    x <- as.data.frame(x)
    return(as.data.frame(lapply(x, robust_scale, group)))
  }
  nn <- length(unique(x))
  group <- as.factor(group)
  lvl2  <- levels(group)[2]

  if(nn>4){
    x <- as.numeric(x)
    x <- x - median(x[which(group==lvl2)], na.rm = TRUE)
    mad_val <- mad(x[which(group==lvl2)], na.rm = TRUE)
    if(mad_val == 0) mad_val <- sd(x[which(group==lvl2)], na.rm = TRUE)
    x <- x / mad_val
  }
  if(nn==2){
    x <- as.numeric(as.factor(x))
    #x <- (x-1.5)*4
    x <- ((x-1.5)*2)*0.675
  }
  if(nn==3){
    x <- as.numeric(as.factor(x))
    #x <- (x-2)*2
    x <- (x-2)*0.675
  }
  if(nn==4){
    x <- as.numeric(as.factor(x))
    #x <- round((x-2.5)*1.33, 2)
    x <- (x-2.5)/1.5
  }
  return(x)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#library(sinew)
#makeOxygen(robust_scale)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# robust_scale <- function(x) {
#   #if (!is.numeric(x)) stop("Input must be numeric.")
#   if (is.matrix(x) || is.data.frame(x)) {
#     x <- as.data.frame(x)
#     return(as.data.frame(lapply(x, robust_scale)))
#   }
#   nn <- length(unique(x))
#
#   if(nn>4){
#     x <- as.numeric(x)
#     x <- x - median(x, na.rm = TRUE)
#     mad_val <- mad(x, na.rm = TRUE)
#     if(mad_val == 0) mad_val <- sd(x, na.rm = TRUE)
#     x <- x / mad_val
#   }
#   if(nn==2){
#     x <- as.numeric(as.factor(x))
#     #x <- (x-1.5)*4
#     x <- ((x-1.5)*2)*0.675
#   }
#   if(nn==3){
#     x <- as.numeric(as.factor(x))
#     #x <- (x-2)*2
#     x <- (x-2)*0.675
#   }
#   if(nn==4){
#     x <- as.numeric(as.factor(x))
#     #x <- round((x-2.5)*1.33, 2)
#     x <- (x-2.5)/1.5
#   }
#   return(x)
# }
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
