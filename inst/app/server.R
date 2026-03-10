


server <- function(input, output, session) {


  output$fileUploaded <- reactive({
    return(!is.null(input$file))
  })
  outputOptions(output, "fileUploaded", suspendWhenHidden = FALSE)


  pool <- reactive({
    req(input$file)
    ext <- tools::file_ext(input$file$name)
    if (ext == "csv") {
      first_lines <- readLines(input$file$datapath, n = 10)
      comma_count <- sum(sapply(first_lines, function(x) length(strsplit(x, ",")[[1]])))
      semicolon_count <- sum(sapply(first_lines, function(x) length(strsplit(x, ";")[[1]])))
      sep <- if(semicolon_count > comma_count) ";" else ","
      read.csv(input$file$datapath, sep = sep, dec = ".")
    } else if (ext == "xlsx") {
      as.data.frame(readxl::read_xlsx(input$file$datapath))
    } else if (ext == "xls") {
      as.data.frame(readxl::read_xls(input$file$datapath))
    } else {
      stop("Unsupported file type")
    }
  })


  tg <- reactive({
    req(input$file2)
    ext <- tools::file_ext(input$file2$name)
    if (ext == "csv") {
      first_lines <- readLines(input$file2$datapath, n = 10)
      comma_count <- sum(sapply(first_lines, function(x) length(strsplit(x, ",")[[1]])))
      semicolon_count <- sum(sapply(first_lines, function(x) length(strsplit(x, ";")[[1]])))
      sep <- if(semicolon_count > comma_count) ";" else ","
      read.csv(input$file2$datapath, sep = sep, dec = ".")
    } else if (ext == "xlsx") {
      as.data.frame(readxl::read_xlsx(input$file2$datapath))
    } else if (ext == "xls") {
      as.data.frame(readxl::read_xls(input$file2$datapath))
    } else {
      stop("Unsupported file type")
    }
  })



  data <- reactive({
    req(pool())
    dat1 <- pool()
    if (!is.null(input$file2)) {
      dat2 <- tg()
      combined <- eVCGsampler::combine_data(dat1, dat2, indicator_name = 'POOL_vs_TG')
    } else {
      combined <- dat1
    }
    return(combined)
  })


  output$outcome_ui <- renderUI({
    req(data())
    datx  <- data()
    uniqs <- apply(datx, 2, FUN=function(x){length(unique(x[which(!is.na(x))]))})
    bin_names <- names(datx)[which(uniqs==2)]
    if ('POOL_vs_TG' %in% bin_names) bin_names <- c('POOL_vs_TG', bin_names[bin_names!='POOL_vs_TG'])
    selectInput("outcome", "Select Binary group variable (POOL vs TG)", choices = bin_names)
  })

  output$stratum_ui <- renderUI({
    req(data())
    datx  <- data()
    uniqs <- apply(datx, 2, FUN=function(x){length(unique(x[which(!is.na(x))]))})
    cat_names <- names(datx)[which(uniqs<=8)]
    cat_names <- cat_names[which(cat_names!='POOL_vs_TG')]
    selectInput("stratum", "Select stratum variable", choices = cat_names, multiple = TRUE)
  })

  output$covariates_ui <- renderUI({
    req(data())
    cnames <- names(data())
    cnames <- cnames[which(cnames!='POOL_vs_TG')]
    selectInput("covariates", "Select covariates to be balanced", choices = cnames, multiple = TRUE)
  })



  observe({
    req(input$outcome,  input$covariates)
    n1 <- strsplit(input$n1, ",")[[1]]
    n1 <- as.integer(trimws(n1))
    n1 <- n1[which(!is.na(n1))]

    c_w <- strsplit(input$c_w, ",")[[1]]
    c_w <- abs(as.numeric(trimws(c_w)))
    c_w <- c_w[which(!is.na(c_w))]
    c_w <- c_w[which(is.finite(c_w))]
    c_w <- c_w[which(c_w>0)]
    if(length(c_w)<2) c_w <- NULL
    if(!is.null(c_w)) c_w <- c_w[1:length(input$covariates)]


    dat <- data()
    dat <- na.omit(dat[,c(input$outcome, input$stratum, input$covariates)])
    if(is.null(input$stratum))  formule <- as.formula(paste0(input$outcome, '~', paste(input$covariates, collapse='+')))
    if(!is.null(input$stratum)) formule <- as.formula(paste0(input$outcome, '~', paste(input$covariates, collapse='+'), ' | ', paste(input$stratum, collapse='+')))

    out  <- eVCGsampler::VCG_sampler(formule, data=dat, n=n1, c_w=c_w, random=input$rand)
    dat2 <- out[[1]]
    vars <- input$covariates

    output$Plot1 <- renderPlot({
      out[[2]]
    })


    output$plots_n <- renderPlot({
      eVCGsampler::BestVCGsize(formule, data=dat)
    })

    output$Plot2 <- renderPlot({
      formula_x <- as.formula(paste0(input$outcome, '~', paste(input$covariates, collapse='+')))
      eVCGsampler::energy_test(formula_x, data=dat, R = 1000)
    })

    output$Plot3 <- renderPlot({
      formule <- as.formula(paste0(input$outcome, '_balanced ~', paste(input$covariates, collapse='+')))
      eVCGsampler::energy_test(formule, data=dat2, R = 1000)
    })

    # Dynamically generate plotOutput elements
    output$plots_ui <- renderUI({
      plot_output_list <- lapply(1:length(input$covariates), function(number) {
        plotOutput(outputId = paste0("splot_", number))
      })
      do.call(tagList, plot_output_list)
    })

    # Dynamically generate renderPlot expressions
    lapply(1:length(input$covariates), function(i) {
      output[[paste0("splot_", i)]] <- renderPlot({
        eVCGsampler::plot_var(dat2, what=input$covariates[i], group='VCG')
      })
    })
    tab <- NULL


    dat3_tg   <- dat2[which(is.na(dat2$VCG)), ]
    dat3_vcg  <- dat2[which(dat2$VCG == 1), ]
    dat3_pool <- dat2[which(!is.na(dat2$VCG)), ]
    dat3_tg$plot_g <- 1
    dat3_vcg$plot_g <- 2
    dat3_pool$plot_g <- 3
    dat3 <- rbind(dat3_tg, dat3_vcg, dat3_pool)
    dat3$plot_g <- factor(dat3$plot_g, levels = c(1, 2, 3), labels = c("TG", "VCG", "POOL"))

    if(length(vars)>1){
      tab_m <- apply(dat3[, vars], 2, FUN=function(x){tapply(x, dat3$plot_g, FUN=median, na.rm=T)})
      tab_mad <- apply(dat3[, vars], 2, FUN=function(x){tapply(x, dat3$plot_g, FUN=mad, na.rm=T)})
      tab <- paste0(signif(tab_m, 3), ' (', signif(tab_mad, 3), ')')
      tab <- matrix(tab, nrow = nrow(tab_m), ncol = ncol(tab_m))
      rownames(tab) <- rownames(tab_m)
      colnames(tab) <- colnames(tab_m)
      tab <- cbind(tab, table(dat3$plot_g))
      colnames(tab)[ncol(tab)]<- 'N'
    }else{
      tab_m<-tapply(dat3[, vars], dat3$plot_g, FUN=median, na.rm=T)
      tab_mad<-tapply(dat3[, vars], dat3$plot_g, FUN=mad, na.rm=T)
      tab <- paste0(signif(tab_m, 3), ' (', signif(tab_mad, 3), ')')
      tab <- cbind(tab, table(dat3$plot_g))
      colnames(tab)<- c(vars, 'N')
      rownames(tab) <- levels(dat3$plot_g)
    }

    output$s_table <- renderTable({
      tab
    }, rownames = TRUE)


    output$downloadData <- downloadHandler(
      filename = function() {
        paste0("eVCGsampler_results-", Sys.Date(), ".csv")
      },
      content = function(file) {
        write.table(dat2, file, row.names = FALSE, sep = ";", dec='.')
      }
    )


  })



}
