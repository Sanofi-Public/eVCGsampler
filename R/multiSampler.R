#' @title Multi-Sample VCG Generator and Overlap Visualization
#' @description
#' Repeatedly samples VCGs (via `VCG_sampler` with `random=TRUE`) from the pool,
#' optionally plots the overlap of VCGs.
#' @param formula A formula specifying the treated and covariates, e.g., `treated ~ cov1 + cov2`. The treated variable must be binary (0=pool, 1=treated)
#' @param data A data frame containing the variables specified in the formula.
#' @param n Integer. Number of observations to sample from the pool. Or a vector of n for each stratum.
#' @param c_w Optional: Vector of positive weights for covariates, reflecting the relative importance of the covariates for the balancing.
#' @param Nsamples Number of VCGs to generate (default is 20).
#' @param plot Logical; if `TRUE`, returns a ggplot2 plot showing the overlap of VCGs (default is `TRUE`).
#'
#' @return
#' If `plot = TRUE`, returns a list with:
#' \describe{
#'   \item{data}{The original data frame with additional VCG columns (`VCG_1`, ..., `VCG_Nsamples`).}
#'   \item{p}{A `ggplot2` object showing the number of times each observation was selected across VCG samples.}
#' }
#' If `plot = FALSE`, returns the modified data frame only.
#' @details
#' The function repeatedly calls `VCG_sampler` with `random` set to TRUE
#' to generate multiple VCGs. It calculates the frequency of selection for each observation
#' and computes the average percentage of overlapping observations.
#' This function should only be used if you really need multiple VCGs, e.g. for PoC studies.
#' It is not intended for selecting one VCG from them afterwards!
#' In this case, the VCG_sampler function should be used directly and only one VCG should be generated.
#'
#' @examples
#'
#'   dat <- data.frame(
#'   treat = rep(0:1, c(50, 30)),
#'   cov_1 = c(rnorm(50, 5, 2),   rnorm(30, 5, 1)),
#'   cov_2 = c(rnorm(50, 11, 2),  rnorm(30, 10, 1))
#'   )
#'
#'   result <- multiSampler(treat ~ cov_1 + cov_2, data = dat, n = 10, Nsamples = 10)
#'
#' @rdname multiSampler
#' @export
#' @importFrom ggplot2 margin scale_y_discrete



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
multiSampler <- function(formula, data, n, c_w=NULL, Nsamples=20, plot=TRUE){
  Value <- ID <- NULL

  if(any(grepl('*', as.character(formula), fixed = TRUE))) stop('Only + for covariates and | for stratum variables are allowed.')
  if(any(grepl(':', as.character(formula), fixed = TRUE))) stop('Only + for covariates and | for stratum variables are allowed.')
  fnames   <- all.vars(as.formula(formula))
  fparts   <- all.names(as.formula(formula))
  leftpart <- fnames[1]   # The response variable
  vars     <- fnames[-1]  # Exclude the response variable
  data[, leftpart] <- as.factor(data[, leftpart])
  data <- data[which(!is.na(data[,leftpart])), ]
  lvls       <- levels(data[, leftpart])
  leftpart_b <- paste0(leftpart, '_balanced')

  if(any(fparts=='|')){
    stratum <- strsplit(as.character(as.formula(formula))[3], "\\|")[[1]]
    stratum <- trimws(strsplit(stratum[2], "\\+")[[1]])
    if(length(stratum)>3) stop('max 3 stratum variables are allowed!')
    vars    <- vars[which(!vars%in%stratum)]
    data$in_stratum <- interaction(data[stratum])
    if(nlevels(data$in_stratum)>8) stop('Too many strata, a maximum of 8 are allowed!')
    if(length(n)!=nlevels(data$in_stratum)) n <- n[1]
    if(length(n)==1) n <- rep(n, nlevels(data$in_stratum))
    new_formula <- as.formula(paste(leftpart, '~', paste(vars, collapse='+')))
    data_out <-  NULL
    p_all <- list()

    for(k in 1:nlevels(data$in_stratum)){
      data_k   <- data[which(data$in_stratum==levels(data$in_stratum)[k]), ]
      out      <- multiSampler(new_formula, data=data_k, n=n[k], c_w=c_w, Nsamples=Nsamples, plot=plot)
      if(plot){
      data_out <- rbind(data_out, out[[1]])
      p_all[[k]] <- out[[2]]
      }else{
        data_out <- rbind(data_out, out)
      }
    }

    if(!plot) return(data_out)
    if(plot)  return(list(data_out, p_all))

    #data  <- na.omit(data[, c(leftpart, vars, stratum)])
  }else{

  for(i in 1:Nsamples){
    out   <- VCG_sampler(formula, data=data, n=n, c_w=c_w, random=TRUE, plot=FALSE)
    if(i==1) data_out <- out
    eval(parse(text=paste0('data_out$VCG_', i, '<- out$VCG')))
  }
  data_out <- data_out[, which(!colnames(data_out)%in%c('VCG', 'e_weights', leftpart_b))]
  dat2 <- data_out[which(data_out[, leftpart]==lvls[1]), ]
  vars <- paste0('VCG_', 1:Nsamples)

  freq <- apply(dat2[, vars], 1, sum)
  p_overlap <- round((mean(freq[which(freq>0)])/Nsamples) * 100,1)


  if(plot){
    df <- data.frame(Value=as.vector(freq), ID=1:nrow(dat2))
    df$ID <- as.factor(df$ID)
    df$ID <- factor(df$ID, levels = rev(levels(df$ID)))

    p <- ggplot(df, aes(x = Value, y = ID)) +
      geom_vline(xintercept = 1,  color = "#00A091", linewidth = 1) +
      geom_bar(stat = "identity", fill = "#7A00E6", colour ='white') +
      geom_vline(xintercept = 0,  color = "#7A00E6", linewidth = 1) +
      geom_vline(xintercept = Nsamples,  color = "#FE7500", linewidth = 1) +
      scale_x_continuous(breaks = seq(0, Nsamples, 1), minor_breaks = NULL) +
      scale_y_discrete(breaks = df$ID,  labels = ifelse(seq_along(df$ID)%%2==0, '', rev(df$ID))) +
      labs(title = "VCGs overlapping", subtitle=paste0('percentage of overlapping: ', p_overlap, '%'), x = "Number of times selected", y = "Index (POOL)") +
      theme_minimal()
    return(list(data_out, p))
  }else{
    return(data_out)
  }
  }

}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#




