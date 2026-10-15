import_server <- function(
        id
) { moduleServer(id, function(input, output, session) {

    ns <- session$ns

    # ------ EXAMPLE DATA --------------------------------------------------
    # The "Load the example data set" link makes the two files that ship in
    # data/ stand in for an upload, so the whole workflow can be walked
    # through without a real data set. Everything downstream reads the two
    # reactives below rather than input$file_* directly, so it cannot tell
    # the difference.
    demo <- reactiveVal(FALSE)

    observeEvent(input$load_example, {
        missing <- example_data[!file.exists(unlist(example_data))]
        if (length(missing) > 0) {
            showNotification(
                paste0("The example files are missing (", paste(unlist(missing), collapse = ", "),
                       "). Regenerate them with source(\"data/make_example_data.R\")."),
                type = "error", duration = 12
            )
            return()
        }
        demo(TRUE)
    })

    # A real upload always wins over the example data. ignoreInit is implicit:
    # observeEvent skips the NULL these inputs start out with.
    observeEvent(input$file_pool, demo(FALSE))
    observeEvent(input$file_tg,   demo(FALSE))

    file_pool <- reactive({
        if (isTRUE(demo()))
            list(name = basename(example_data$pool),
                 datapath = normalizePath(example_data$pool))
        else input$file_pool
    })

    file_tg <- reactive({
        if (isTRUE(demo()))
            list(name = basename(example_data$tg),
                 datapath = normalizePath(example_data$tg))
        else input$file_tg
    })

    # ------ RAW UPLOADS ---------------------------------------------------
    # Read errors (wrong extension, unreadable xlsx, ...) are caught and
    # reported in the UI instead of taking the session down.
    read_err <- reactiveValues(pool = NULL, tg = NULL)

    pool <- reactive({
        req(file_pool())
        out <- tryCatch(
            read_upload(file_pool()$datapath, file_pool()$name),
            error = function(e) e
        )
        if (inherits(out, "error")) {
            read_err$pool <- conditionMessage(out)
            return(NULL)
        }
        read_err$pool <- NULL
        out
    })

    tg <- reactive({
        req(file_tg())
        out <- tryCatch(
            read_upload(file_tg()$datapath, file_tg()$name),
            error = function(e) e
        )
        if (inherits(out, "error")) {
            read_err$tg <- conditionMessage(out)
            return(NULL)
        }
        read_err$tg <- NULL
        out
    })

    # ------ COMBINING -----------------------------------------------------
    # With a single upload the file is used as it is; with two uploads
    # combine_data() stacks them on the intersection of the column names and
    # adds the POOL / TG indicator.
    combine_error <- reactiveVal(NULL)

    data <- reactive({
        p <- pool()
        req(p)
        if (is.null(file_tg())) {
            combine_error(NULL)
            return(p)
        }
        t <- tg()
        req(t)
        out <- tryCatch(
            withCallingHandlers(
                eVCGsampler::combine_data(p, t, indicator_name = group_indicator),
                warning = function(w) invokeRestart("muffleWarning")
            ),
            error = function(e) e
        )
        if (inherits(out, "error")) {
            combine_error(conditionMessage(out))
            return(NULL)
        }
        combine_error(NULL)
        out
    })

    # Which variable came from where - shown in the overview table.
    origin <- reactive({
        p <- pool()
        req(p)
        if (is.null(file_tg())) {
            return(stats::setNames(rep("single file", ncol(p)), names(p)))
        }
        t <- tg()
        req(t)
        nms <- union(names(p), names(t))
        o <- ifelse(nms %in% names(p) & nms %in% names(t), "both",
                    ifelse(nms %in% names(p), "POOL only", "TG only"))
        o <- stats::setNames(o, nms)
        c(o, stats::setNames("created", group_indicator))
    })

    # The TG data being bigger than the POOL is the classic "files swapped"
    # mistake, so it is called out explicitly.
    swapped <- reactive({
        if (is.null(file_tg())) return(FALSE)
        p <- pool(); t <- tg()
        !is.null(p) && !is.null(t) && nrow(t) > nrow(p)
    })

    # ------ MESSAGES ------------------------------------------------------
    output$upload_msg <- renderUI({
        msgs <- list()
        # Flagged up front so a demo run can never be mistaken for a real one.
        # The file names carry through to the report as well.
        if (isTRUE(demo()))
            msgs <- c(msgs, list(div(
                class = "vcg-note note-info",
                HTML(paste0(pill("example data", "neutral"),
                            " Simulated data for trying the app out - not a real study. ",
                            "Upload a POOL file to switch to your own data."))
            )))
        if (!is.null(read_err$pool))
            msgs <- c(msgs, list(div(class = "vcg-note note-bad",
                                     paste0("POOL file could not be read: ", read_err$pool))))
        if (!is.null(read_err$tg))
            msgs <- c(msgs, list(div(class = "vcg-note note-bad",
                                     paste0("TG file could not be read: ", read_err$tg))))
        if (!is.null(combine_error()))
            msgs <- c(msgs, list(div(class = "vcg-note note-bad",
                                     paste0("The two files could not be combined: ", combine_error()))))
        if (isTRUE(swapped()))
            msgs <- c(msgs, list(div(
                class = "vcg-note note-info",
                HTML(paste0(pill("check the order", "warn"),
                            " The TG file has more rows than the POOL file. Normally the ",
                            "historical control pool is the larger one - make sure the POOL ",
                            "went into the first field."))
            )))
        if (length(msgs) == 0) return(NULL)
        div(class = "mt-3", msgs)
    })

    output$combine_msg <- renderUI({
        req(file_tg())
        chk <- combine_check(pool(), tg())
        req(chk)
        n_drop <- sum(chk$Status == "dropped")
        n_bad  <- sum(chk$Status == "type mismatch")
        n_keep <- sum(chk$Status == "combined")
        div(
            class = "mb-3",
            HTML(paste(
                pill(paste0(n_keep, " combined"), if (n_keep > 0) "ok" else "bad"),
                if (n_drop > 0) pill(paste0(n_drop, " dropped"), "warn"),
                if (n_bad  > 0) pill(paste0(n_bad, " type mismatch"), "bad")
            ))
        )
    })

    # ------ TILES ---------------------------------------------------------
    output$tiles <- renderUI({
        d <- data()
        req(d)
        p <- pool()
        n_tg <- if (is.null(file_tg())) NA_integer_ else nrow(tg())
        div(
            class = "mt-3",
            tiles(
                tile("POOL file", format(nrow(p), big.mark = " "), "rows uploaded"),
                tile("TG file",
                     if (is.na(n_tg)) "-" else format(n_tg, big.mark = " "),
                     if (is.na(n_tg)) "included in the POOL file" else "rows uploaded"),
                tile("Combined", format(nrow(d), big.mark = " "), "rows for balancing", accent = TRUE),
                tile("Variables", ncol(d), "columns available")
            )
        )
    })

    # ------ TABLES --------------------------------------------------------
    output$dt_vars <- DT::renderDataTable({
        d <- data()
        req(d)
        vcg_table(variable_overview(d, origin()), dom = "tp", page = 15, scrollx = TRUE)
    })

    output$dt_combine <- DT::renderDataTable({
        req(file_tg())
        chk <- combine_check(pool(), tg())
        req(chk)
        vcg_table(chk, dom = "tp", page = 15, scrollx = TRUE)
    })

    output$dt_data <- DT::renderDataTable({
        d <- data()
        req(d)
        vcg_table(d, dom = "tp", page = 10, scrollx = TRUE)
    })

    # ------ GATE ----------------------------------------------------------
    # Continue needs a data set that at least one variable pair could be
    # balanced on: a two-valued group variable plus one usable covariate.
    ready <- reactive({
        d <- data()
        if (is.null(d) || nrow(d) == 0) return(FALSE)
        length(binary_vars(d)) > 0 && length(covariate_vars(d)) > 0
    })

    observe({
        d <- data()
        if (is.null(d)) {
            shinyjs::hide("overview_card")
            shinyjs::disable("continue")
            return()
        }
        shinyjs::show("overview_card")
        if (isTRUE(ready())) shinyjs::enable("continue") else shinyjs::disable("continue")
    })

    observeEvent(ready(), {
        if (!isTRUE(ready()) && !is.null(data())) {
            showNotification(
                paste("This data set has no variable with exactly two distinct values.",
                      "The sampler needs a POOL / TG group variable - upload the treated",
                      "animals as a second file, or add a two-level group column."),
                type = "warning", duration = 12
            )
        }
    })

    # ------ RETURNS -------------------------------------------------------
    list(
        continue  = reactive(input$continue),
        ready     = ready,
        data      = data,
        pool_name = reactive(file_pool()$name),
        tg_name   = reactive(file_tg()$name),
        n_pool    = reactive(if (is.null(pool())) NA_integer_ else nrow(pool())),
        n_tg      = reactive(if (is.null(file_tg()) || is.null(tg())) NA_integer_ else nrow(tg())),
        combined  = reactive(!is.null(file_tg())),
        check     = reactive(if (is.null(file_tg())) NULL else combine_check(pool(), tg()))
    )

})}
