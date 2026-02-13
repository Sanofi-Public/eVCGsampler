#' @title Visualize Covariate Distribution Across TG, VCG, and POOL
#' @description Creates a plot to compare the distribution of a selected variable across three groups:
#' TG (treated groups), VCG (virtual control group), and POOL (data pool).
#' @param data A balanced data frame (output of the VCG_sampler function)
#' @param what A string specifying the name of the variable to be visualized.
#' @param stratum A string specifying the name of the stratum variable (default is `"in_stratum"`)
#' @param group A string specifying the column name used to define group membership (default is `"VCG"`).
#' @param title Optional title for the plot.
#' @return A ggplot2 object showing either:
#' \itemize{
#'   \item A boxplot for continuous variables (more than 4 unique values).
#'   \item A proportional bar chart for categorical variables (2–4 unique values).
#' }
#' @details The function uses energy distance to quantify distributional differences between groups.
#' For continuous variables, it overlays dashed lines for TG group statistics (mean, min, max) and displays sample sizes.
#' For categorical variables, it uses color-coded bars and cumulative proportion lines to highlight imbalance.
#'
#' @examples
#'
#' dat   <- data.frame(
#'   cov1  = rnorm(50, 10, 1),
#'   cov2  = rnorm(50, 7,  1),
#'   cov3  = rnorm(50, 5,  1),
#'   treated = rep(c(0, 1), c(35, 15))
#' )
#'   out <- VCG_sampler(treated ~ cov1 + cov2 + cov3, data=dat, n=5, plot=FALSE)
#'   plot_var(out, what='cov1', group='VCG')
#'   plot_var(out, what='cov2', group='VCG')
#'
#' @rdname plot_var
#' @export
#' @importFrom ggplot2 ggplot geom_boxplot geom_bar geom_hline geom_text labs theme_minimal scale_fill_manual aes stat_summary theme position_dodge
#' @importFrom stats na.omit
#' @importFrom patchwork plot_layout
#'



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
plot_var <- function(data, what=NULL, stratum='in_stratum', group='VCG', title=''){
  p1<-p2<-p3<-p4<-p5<-p6<-p7<-p8<- plot_g <- plot_x <- count <- group_c <- Var1 <- Freq <- cum_prop <- NULL
  hh <- NULL
  #grid.arrange(p1, p2, ncol = 2)
  if(!is.null(stratum)){ if(!stratum%in%names(data)) stratum <- NULL }


  if(!is.null(stratum)){
    data$in_stratum <- data[,stratum]
    nl <- nlevels(data$in_stratum)
    for(k in 1:nl){
      data_k <- data[which(data$in_stratum==levels(data$in_stratum)[k]), ]
      p0     <- plot_var(data_k, what=what, stratum='no_stratum', group=group, title=levels(data$in_stratum)[k])
      eval(parse(text=paste0('p', k, '<- p0')))
    }

    #require(patchwork)
    if(nl==2) (p <- p1 / p2)
    if(nl==3) (p <- p1 / p2 / p3)
    if(nl==4) (p <- (p1 | p2) / (p3 | p4))
    if(nl==5) (p <- (p1 | p2) / (p3 | p4) / (p5))
    if(nl==6) (p <- (p1 | p2) / (p3 | p4) / (p5 | p6))
    if(nl==7) (p <- (p1 | p2) / (p3 | p4) / (p5 | p6)/ (p7))
    if(nl==8) (p <- (p1 | p2) / (p3 | p4) / (p5 | p6)/ (p7 | p8))
    return((p))

  }
  if(is.null(stratum)){
    {
      if(nchar(title)>0)  title_text <- paste0(title, ': ', what)
      if(nchar(title)==0) title_text <- what


      data$left_p <- data[, grep("_balanced", names(data), value = TRUE)]
      data$plot_x <- as.numeric(data[, what])
      data_tg   <- data[which(is.na(data$VCG) & !is.na(data$plot_x) & !is.na(data$left_p)), ]
      data_vcg  <- data[which(data$VCG==1), ]
      data_pool <- data[which(!is.na(data$VCG)), ]
      data_tg$plot_g   <-1
      data_vcg$plot_g  <-2
      data_pool$plot_g <-3
      data <- rbind(data_tg, data_vcg, data_pool)
      data$plot_g <- factor(data$plot_g, levels = c(1, 2, 3), labels = c('TG', 'VCG', 'POOL'))
      data <- na.omit(data[,c('plot_x', 'plot_g')])

      n_df <- data.frame(table(data$plot_g))
      n_df$group_c <-"POOL_1"
      stat_n <- function(y){return(data.frame(y = -Inf, label = paste0('n=', sum(!is.na(y)))))}



      dist1 <- round(energy_distance(plot_g ~ plot_x, data=data[which(data$plot_g!='VCG'), ]),  3)
      dist2 <- round(energy_distance(plot_g ~ plot_x, data=data[which(data$plot_g!='POOL'), ]), 3)
      dist_text <- paste('Distance: ', dist1, ' (TG-POOL) >> ', dist2, ' (TG-VCG)')

      if(length(unique(data$plot_x))>4){

        ggplot(data, aes(x = plot_g, y = plot_x, fill = plot_g)) +
          geom_boxplot(outlier.size = 3) +
          geom_hline(yintercept = mean(data$plot_x[which(data$plot_g=='TG')], na.rm=T), linetype = "dashed", color = "#7A00E6") +
          geom_hline(yintercept = min(data$plot_x[which(data$plot_g=='TG')],  na.rm=T), linetype = "dashed", color = "#7A00E6") +
          geom_hline(yintercept = max(data$plot_x[which(data$plot_g=='TG')],  na.rm=T), linetype = "dashed", color = "#7A00E6") +
          stat_summary(fun = mean, geom = "point", shape = 21, size = 5, color =  1, bg='white') +
          scale_fill_manual(values = c("TG" = "#7A00E6", "VCG" = "#FE7500", "POOL" = "#00A091")) +
          theme_minimal() +  theme(legend.position = "none") + labs(title = title_text, subtitle = dist_text) +
          xlab('') + ylab("value") +
          stat_summary(aes(group=plot_g), fun.data = stat_n, geom="text", position=position_dodge(width = 0.75), hjust =  0.5, vjust = -0.5, col=rgb(0.4, 0.4, 0.4))
      }else{

        data$plot_x <- as.factor(as.numeric(as.factor(data$plot_x)))
        nn <-nlevels(data$plot_x)

        if(nn==2){
          custom_colors <- c(
            "TG_1" =  "#A366FF",  "TG_2" = "#5C00B0",
            "VCG_1" = "#FFA64D", "VCG_2" = "#C75E00",
            "POOL_1"= "#33C6B8","POOL_2" = "#00796B")
        }else{
          custom_colors <- c(
            "TG_1" = "#C080FF",    "TG_2" = "#A366FF",  "TG_3" = "#5C00B0",  "TG_4" = "#4B0080",
            "VCG_1" = "#FFD199",   "VCG_2" = "#FFA64D", "VCG_3" = "#C75E00", "VCG_4" = "#994C00",
            "POOL_1" = "#66E0D4",  "POOL_2"= "#33C6B8", "POOL_3"= "#00796B", "POOL_4" = "#00564F")
        }

        df <- as.data.frame(prop.table(table(data$plot_g, data$plot_x), margin = 1))
        colnames(df) <- c("plot_g", "plot_x", "count")
        df$group_c <- paste0(df$plot_g, '_', df$plot_x)

        group1_props <- df[which(df$plot_g=='TG'), ]
        group1_props$cum_prop <- cumsum(rev(group1_props$count))
        group1_props <- group1_props[which(group1_props$cum_prop<1),]
        ggplot(df, aes(x = plot_g, y = count, fill = group_c)) +
          geom_bar(stat = "identity") +
          geom_hline(data = group1_props, aes(yintercept = cum_prop), linetype = "dashed", color = "#7A00E6") +
          scale_fill_manual(values = custom_colors) +
          labs(title = title_text, subtitle = dist_text) +
          theme_minimal() + xlab('') + ylab("proportion") + theme(legend.position = "none") +
          geom_text(data = n_df, aes(x = Var1, y = 0, label = paste0("n=", Freq)),
                    vjust = 1.5, fontface = "italic", size = 3.5, col=rgb(0.4, 0.4, 0.4))

      }
    }
  }
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

