setup_server <- function(
        id,
        data
) { moduleServer(id, function(input, output, session) {

    ns <- session$ns

    # ------ VARIABLE PICKERS ----------------------------------------------
    # Rebuilt whenever a new file is uploaded. The choices are restricted to
    # what eVCGsampler can actually use, so an unusable column cannot be
    # selected in the first place.

    # The animal ID. Pre-selected with the first complete, per-group unique
    # column so the common case needs no interaction. It takes no part in the
    # balancing and only travels through to the results, so a sampled animal
    # can be traced back to the source data.
    output$id_ui <- renderUI({
        d <- data()
        req(d)
        ch <- id_vars(d, input$group %||% group_indicator)
        selectInput(
            inputId  = ns("id"),
            label    = "Animal ID (identifies the sampled animals)",
            choices  = c("(none)" = "", ch),
            selected = if (length(ch)) ch[1] else ""
        )
    })

    output$group_ui <- renderUI({
        d <- data()
        req(d)
        ch <- binary_vars(d)
        # combine_data() puts the indicator in front - keep it as the default.
        if (group_indicator %in% ch) ch <- c(group_indicator, setdiff(ch, group_indicator))
        selectInput(
            inputId  = ns("group"),
            label    = "Binary group variable (POOL vs TG)",
            choices  = ch,
            selected = ch[1]
        )
    })

    output$covariates_ui <- renderUI({
        d <- data()
        req(d)
        ch <- setdiff(covariate_vars(d), c(group_indicator, input$group, input$id))
        selectizeInput(
            inputId  = ns("covariates"),
            label    = "Covariates to balance",
            choices  = ch,
            multiple = TRUE,
            options  = multi_select_options
        )
    })

    output$stratum_ui <- renderUI({
        d <- data()
        req(d)
        ch <- setdiff(stratum_vars(d, max_strata),
                      c(group_indicator, input$group, input$id, input$covariates))
        selectizeInput(
            inputId  = ns("stratum"),
            label    = paste0("Stratum variables (optional, max ", max_strata, " strata in total)"),
            choices  = ch,
            multiple = TRUE,
            options  = multi_select_options
        )
    })

    # ------ DERIVED SELECTION ---------------------------------------------
    # "" (the "(none)" choice) is normalised to NULL so downstream code has a
    # single "no ID" representation to test for.
    id_var     <- reactive(if (is.null(input$id) || !nzchar(input$id)) NULL else input$id)
    covariates <- reactive(input$covariates)
    stratum    <- reactive(input$stratum)
    sizes      <- reactive(parse_sizes(input$sizes))
    weights    <- reactive(parse_weights(input$weights, length(covariates())))

    # The complete-case data set for the selected variables. NULL rather than
    # req() so the Continue gate below can read it without being interrupted
    # by a silent error when nothing is selected yet.
    adata <- reactive({
        d <- tryCatch(data(), error = function(e) NULL)
        if (is.null(d) || is.null(input$group) || length(covariates()) == 0) return(NULL)
        analysis_data(d, input$group, covariates(), stratum(), id = id_var())
    })

    overview <- reactive({
        ad <- adata()
        if (is.null(ad) || nrow(ad) == 0) return(NULL)
        size_overview(ad, input$group, stratum(), sizes())
    })

    ftext <- reactive(formula_text(input$group, covariates(), stratum()))

    # ------ FORMULA -------------------------------------------------------
    output$formula_box <- renderUI({
        txt <- ftext()
        if (!nzchar(txt)) {
            return(div(class = "vcg-note note-info",
                       "Select at least one covariate to build the formula."))
        }
        div(class = "vcg-formula", txt)
    })

    output$weights_msg <- renderUI({
        if (!nzchar(trimws(input$weights %||% ""))) return(NULL)
        w <- weights()
        n <- length(covariates())
        if (is.null(w)) {
            div(class = "mt-3",
                HTML(paste0(
                    pill("weights ignored", "warn"),
                    " Give one positive weight per covariate (", n,
                    " needed). At least two weights are required - a single value would",
                    " scale every covariate equally and change nothing."
                )))
        } else {
            div(class = "mt-3",
                HTML(paste0(
                    pill("weights applied", "ok"), " ",
                    paste(paste0(covariates(), " x", w), collapse = ", "),
                    " - a higher weight pulls that covariate closer to the TG at the",
                    " expense of the others."
                )))
        }
    })

    # ------ SIZE OVERVIEW -------------------------------------------------
    output$tiles <- renderUI({
        d  <- tryCatch(data(), error = function(e) NULL)
        ad <- adata()
        ov <- overview()
        req(d, ad, ov)
        n_drop <- nrow(d) - nrow(ad)
        tiles(
            tile("POOL", format(sum(ov$N_POOL), big.mark = " "), "animals to sample from"),
            tile("TG", format(sum(ov$N_TG), big.mark = " "), "animals to balance towards"),
            tile("VCG", if (all(is.na(ov$VCG_size))) "-" else format(sum(ov$VCG_size), big.mark = " "),
                 "animals to be selected", accent = TRUE),
            tile("Strata", nrow(ov),
                 if (length(stratum()) > 0) paste(stratum(), collapse = " x ") else "no stratification"),
            tile("Dropped", format(n_drop, big.mark = " "), "rows with a missing value")
        )
    })

    output$dt_sizes <- DT::renderDataTable({
        ov <- overview()
        req(ov)
        show <- ov
        show$Share_of_POOL <- ifelse(is.na(show$Share_of_POOL), "-",
                                     paste0(show$Share_of_POOL, "%"))
        show$VCG_size <- ifelse(is.na(show$VCG_size), "-", as.character(show$VCG_size))
        show$Status <- vapply(ov$Status, function(s) pill(s, status_level(s)), character(1))
        names(show) <- c("Stratum", "N POOL", "N TG", "VCG size", "Share of POOL", "Status")
        vcg_table(show, dom = "t", escape = FALSE)
    })

    output$size_msg <- renderUI({
        ov <- overview()
        req(ov)
        msg <- size_messages(ov, length(covariates()))
        if (length(msg) == 0) {
            return(div(class = "mt-3",
                       HTML(paste0(pill("ready to test", "ok"),
                                   " The selection is usable - continue to the pre-tests."))))
        }
        div(class = "mt-3", lapply(msg, function(m) div(class = "vcg-note", m)))
    })

    # ------ GATE ----------------------------------------------------------
    ready <- reactive(isTRUE(size_is_runnable(overview())))

    observe({
        if (ready()) shinyjs::enable("continue") else shinyjs::disable("continue")
    })

    # ------ RETURNS -------------------------------------------------------
    list(
        continue     = reactive(input$continue),
        ready        = ready,
        group        = reactive(input$group),
        id           = id_var,
        covariates   = covariates,
        stratum      = stratum,
        sizes        = sizes,
        weights      = weights,
        adata        = adata,
        overview     = overview,
        formula_text = ftext
    )

})}
