#' @title Permutation Energy Test for Covariate Imbalance
#' @description Performs a permutation-based energy distance test to assess whether two groups (defined by a binary treated variable)
#' are balanced across a set of covariates. Optionally, it visualizes the distribution of permuted energy distances
#' and highlights the observed test statistic and critical value.
#' @param formula A formula specifying the treated and covariates, e.g., `treated ~ cov1 + cov2 | stratum`.
#' @param data A data frame containing the variables specified in the formula.
#' @param alpha Significance level for the test (default is 0.05).
#' @param R  Number of permutations to perform (default is 2000).
#' @param plot Logical. If `TRUE`, returns a ggplot2 visualization of the permutation distribution.
#' @return If `plot = TRUE`, returns a list with:
#' \itemize{
#'   \item A list of class `"htest"` containing:
#'     \itemize{
#'       \item `p.value`: The permutation p-value.
#'       \item `estimate`: The observed energy distance.
#'       \item `critical.value`: The critical value at the specified alpha level.
#'       \item `alternative`: The alternative hypothesis ("one.sided").
#'       \item `method`: Description of the test.
#'       \item `n.permutations`: Number of permutations performed.
#'       \item `permutations`: Vector of permuted energy distances.
#'     }
#'   \item A ggplot2 object showing the histogram of permuted distances, with vertical lines for the observed
#'         statistic and critical value.
#' }
#' If `plot = FALSE`, returns only the `"htest"` result list.
#' @details
#' The energy distance is a non-parametric measure of distributional difference. This test evaluates whether
#' the covariate distributions between two groups are statistically distinguishable. A small p-value indicates
#' imbalance between groups. A one-sided test is used because the energy distance is strictly positive; only values greater than the observed statistic in the permutation distribution are relevant.
#' @examples
#'
#' dat <- data.frame(
#'  treated = rep(0:1, c(50, 30)),
#'  age    = c(rnorm(50, 5, 2),   rnorm(30, 5, 1)),
#'  weight = c(rnorm(50, 11, 2),  rnorm(30, 10, 1)),
#'  class  = c(rbinom(50, 3, 0.6),   rbinom(30, 3, 0.4))
#'  )
#'
#'  energy_test(treated ~ age + weight + class, data=dat, R = 500)
#'
#' @seealso
#'  \code{\link[ggplot2]{element}}
#' @rdname energy_test
#' @export
#' @importFrom ggplot2 margin ggplot geom_histogram geom_vline annotate scale_x_continuous theme_minimal margin labs element_text coord_cartesian after_stat aes theme
#' @importFrom stats quantile



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Permutation energy test (test if samples are balanced)
energy_test <- function(formula, data, alpha=0.05, R=2000, plot=TRUE){
  if(any(grepl('*', as.character(formula), fixed = TRUE))) stop('Only + for covariates and | for stratum variables are allowed.')
  if(any(grepl(':', as.character(formula), fixed = TRUE))) stop('Only + for covariates and | for stratum variables are allowed.')
  fnames   <- all.vars(as.formula(formula))
  fparts   <- all.names(as.formula(formula))
  leftpart <- fnames[1]   # The response variable
  vars     <- fnames[-1]  # Exclude the response variable

  if(any(fparts=='|')){
    if(plot) warning('Plot is not supported with stratum.')
    stratum <- strsplit(as.character(as.formula(formula))[3], "\\|")[[1]]
    stratum <- trimws(strsplit(stratum[2], "\\+")[[1]])
    if(length(stratum)>3) stop('max 3 stratum variables are allowed!')
    vars    <- vars[which(!vars%in%stratum)]
    data$in_stratum <- interaction(data[stratum])
    if(nlevels(data$in_stratum)>8) stop('Too many strata, a maximum of 8 are allowed!')
    new_formula <- as.formula(paste(leftpart, '~', paste(vars, collapse='+')))
    tab <-  NULL

    for(k in 1:nlevels(data$in_stratum)){
      data_k <- data[which(data$in_stratum==levels(data$in_stratum)[k]), ]
      out    <- energy_test(new_formula, data=data_k, alpha=alpha, R=R, plot=FALSE)
      tab  <- rbind(tab, c(levels(data$in_stratum)[k], out$estimate, out$critical.value[2], out$p.value))
    }

    colnames(tab) <- c('Stratum', 'estimate', 'critical.value', 'p-value')
    return(tab)
  }else{
    data <- na.omit(data[, c(leftpart, vars)])
    obs_dists  <- energy_distance(formula, data)
    dists <- rep(NA, R)
    for(i in 1:R){
      data$per_group  <- sample(data[, leftpart], replace = FALSE)
      dists[i]  <- energy_distance(as.formula(paste('per_group ~', paste(vars, collapse='+'))), data)
    }
    range   <- quantile(dists, probs = 1-alpha, na.rm = TRUE)
    p_value <- mean(dists >= obs_dists)

    result <- list(
      p.value = p_value,
      estimate = c(e_distance = round(obs_dists, 4)),
      critical.value = c(paste0((1-alpha)*100, '%: '), round(range, 4)),
      alternative = "one.sided",
      method = "Permutation-test for energy distance",
      n.permutations = R,
      permutations = dists
    )
    class(result) <- "htest"

    if(plot){
      s_colors <- c("#00A091", "#FE7500", "#7A00E6", "#F5C142","#694A3E","#1F70A6","#EBA0AB","#6CADD1")
      df <- data.frame(dists = dists)
      p <- ggplot(df, aes(x = dists)) +
        geom_histogram(aes(y = after_stat(density)),
                       bins = 30,
                       fill = "#7A00E6",
                       color = "white") +
        geom_vline(xintercept = range,      color = "#FE7500", linewidth = 1.2) +
        geom_vline(xintercept = obs_dists,  color = "#00A091", linewidth = 1.2) +
        scale_x_continuous(name = "Energy Distance", limits = c(min(c(obs_dists, dists)), max(c(obs_dists, dists)))) +
        theme_minimal() +
        theme(plot.margin = ggplot2::margin(1, 0.5, 0.5, 0.5, "cm"), plot.subtitle = element_text(margin = ggplot2::margin(b = 15))) +
        labs(title="Permutation test", subtitle = paste('p-value: ', signif(p_value, digits=3))) +
        coord_cartesian(clip = "off") +
        annotate("text", x = range, y = Inf, label = paste0((1-alpha) * 100, "%: ", round(range, 2)),
                 vjust = -1.2, color = "#FE7500", size = 4) +
        annotate("text", x = obs_dists, y = Inf, label = round(obs_dists, 2),
                 vjust = -0.1, color = "#00A091",  size = 4)
      return(list(result, p))
    }else{
      return(result)
    }
  }
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
