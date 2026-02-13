#' @title The function attempts to find the optimal size for VCG.
#' @description The function tries out different sizes of VCG and searches for the smallest distance.
#' @param formula A formula specifying the treated and covariates, e.g., `treated ~ cov1 + cov2 | stratum`. The treated variable must be binary (0=pool, 1=treated)
#' @param data A data frame containing the variables specified in the formula.
#' @param plot Logical. If `TRUE`, returns a ggplot2 plot. Default: TRUE
#' @return If `plot = TRUE`, returns a list with:
#' \item{optimal_n}{The estimated optimal VCG size (integer).}
#' \item{plot}{A ggplot2 object visualizing the energy distance curve and plateau.}
#' @details It is only intended for exploratory purposes, as the VCG size is normally given.
#' But it can be used to see how well the given size fits.
#' The recommendation for VCG size is based solely on distance and does not take into account other aspects such as power or validity.
#' @examples
#'
#' set.seed(2342)
#' dat <- data.frame(
#'   treat = rep(0:1, c(50, 30)),
#'   cov1 = c(rnorm(50, 11, 2),  rnorm(30, 10, 1)),
#'   cov2 = c(rnorm(50, 12, 2),  rnorm(30, 10, 1)),
#'   cov3 = c(rnorm(50, 9,  2),  rnorm(30, 10, 1))
#' )
#'  BestVCGsize(treat ~ cov1 + cov2 + cov3, data=dat)
#'
#' @rdname BestVCGsize
#' @importFrom stats lm summary.lm
#' @importFrom ggplot2 geom_line scale_color_manual
#' @export



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
BestVCGsize <- function(formula, data=data, plot=TRUE){
  n <- NULL
  distance <- NULL

  # Function to find plateau using t-tests
  find_increase <- function(data, alpha = 0.05, maxN) {
    pval <- slope <- rep(NA, nrow(data)-2)
    for (i in 3:(nrow(data)-2)) {
      fit     <- summary(lm(distance ~ n, data=data[i:min(i+round(maxN/10), nrow(data)), ]))
      pval[i] <- fit$coef[2, 4]
      slope[i] <- fit$coef[2, 1]
    }
    at <- which(pval<0.05 & slope>0)[1]
    return(min(at+round(maxN/20), maxN))
  }
  # Function to find plateau using t-tests
  find_plateau <- function(data, alpha = 0.05) {
    pval <- slope <- rep(NA, nrow(data)-2)
    for (i in 1:(nrow(data)-2)) {
      fit     <- summary(lm(distance ~ n, data=data[i:nrow(data), ]))
      pval[i] <- fit$coef[2, 4]
      slope[i] <- fit$coef[2, 1]
    }
    return(which(pval>0.05)[1])
  }

  if(any(grepl('*', as.character(formula), fixed = TRUE))) stop('Only + for covariates and | for stratum variables are allowed.')
  if(any(grepl(':', as.character(formula), fixed = TRUE))) stop('Only + for covariates and | for stratum variables are allowed.')
  fnames   <- all.vars(as.formula(formula))
  fparts   <- all.names(as.formula(formula))
  leftpart <- fnames[1]   # The response variable
  vars     <- fnames[-1]  # Exclude the response variable
  leftpart_b <- paste0(leftpart, '_balanced')

  if(any(fparts=='|')){
    stratum <- strsplit(as.character(as.formula(formula))[3], "\\|")[[1]]
    stratum <- trimws(strsplit(stratum[2], "\\+")[[1]])
    if(length(stratum)>3) stop('max 3 stratum variables are allowed!')
    vars    <- vars[which(!vars%in%stratum)]
    data$in_stratum <- interaction(data[stratum])
    if(nlevels(data$in_stratum)>8) stop('Too many strata, a maximum of 8 are allowed!')
    new_formula <- as.formula(paste(leftpart, '~', paste(vars, collapse='+')))
    f2  <- as.formula(paste0(leftpart, '_balanced ~', paste(vars, collapse='+')))
    tab <- dall <- NULL
    cuts <- rep(NA, nlevels(data$in_stratum))
    mnn <- 0

    for(k in 1:nlevels(data$in_stratum)){
      data_k <- data[which(data$in_stratum==levels(data$in_stratum)[k]), ]
      maxN     <- table(data_k[,leftpart])[1]
      Dist     <- rep(NA, maxN)
      Ns <- 1:maxN
      for(i in 1:maxN){
        out     <- VCG_sampler(new_formula, data=data_k, n=Ns[i], plot=FALSE)
        Dist[i] <- energy_distance(f2, data=na.omit(out[, c(paste0(leftpart, '_balanced'), vars)]))
      }
      df <- data.frame(distance=Dist, n=Ns, stratum=levels(data$in_stratum)[k])
      mn <- find_increase(df, maxN=maxN)
      if(!is.na(mn)) df <- df[1:mn, ]
      plateau  <- find_plateau(df)
      if(plateau<3) plateau <- NA
      best_n <- NA
      if(!is.na(plateau)) best_n <- df$n[plateau]
      cuts[k] <- best_n
      tab <- rbind(tab, c(levels(data$in_stratum)[k], maxN, best_n))
      dall<- rbind(dall, df)
      mnn <- max(mnn, df$n, na.rm=T)
    }
    colnames(tab) <- c('Stratum', 'available N', 'optimal VCG N')

    if(plot){
      s_colors <- c("#00A091", "#FE7500", "#7A00E6", "#F5C142","#694A3E","#1F70A6","#EBA0AB","#6CADD1")
      p <- ggplot(dall, aes(x = n, y = distance, color=stratum)) +
        geom_line(linewidth = 1) +
        geom_point(size=3) +
        geom_vline(xintercept = cuts, linetype = "dashed", linewidth = rep(1, length(cuts)), color = s_colors[1:length(cuts)]) +
        labs(title = "Best VCG size", x = "VCG size (n)", y = "Energy distance") +
        theme_minimal() +
        scale_color_manual(values = s_colors[1:length(cuts)])+
        scale_x_continuous(breaks = seq(1, max(mnn), by = 2))
      return(list(tab, p))
    }else{
      return(tab)
    }


  }else{

    f2       <- as.formula(paste0(leftpart_b, '~', paste(vars, collapse='+')))
    maxN     <- table(data[,leftpart])[1]
    Dist     <- rep(NA, maxN)
    Ns       <- 1:maxN
    data     <- na.omit(data[, c(leftpart, vars)])
    for(i in 1:maxN){
      out     <- VCG_sampler(formula, data=data, n=Ns[i], plot=FALSE)
      Dist[i] <- energy_distance(f2, data=na.omit(out[, c(leftpart_b, vars)]))
    }

    df <- data.frame(distance=Dist, n=Ns)
    mn  <- find_increase(df, maxN=maxN)
    if(!is.na(mn)) df       <- df[1:mn, ]
    plateau  <- find_plateau(df)
    if(plateau<3) plateau <- NA

    if(plot){
      p <- ggplot(df, aes(x = n, y = distance)) +
        geom_line(color = "#7A00E6", linewidth = 1) +
        geom_point(color = "#7A00E6", size=3) +
        {if(!is.na(plateau)) geom_vline(xintercept = df$n[plateau], linetype = "dashed", linewidth = 1, color = "#FE7500")} +
        {if(!is.na(plateau)) annotate("text", x = df$n[plateau]+0.5, y = Inf,
                                      label = paste("No significant\nimprovement\nafter n =", df$n[plateau]),
                                      color = "#FE7500", hjust = 0, vjust = 1.0, size=5)} +
        labs(title = "Best VCG size",
             x = "VCG size (n)", y = "Energy distance") +
        theme_minimal() +
        scale_x_continuous(breaks = seq(1, max(df$n), by = 2))
      return(list(optimal_n=df$n[plateau], p))
    }
    return(optimal_n=df$n[plateau])
  }
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

