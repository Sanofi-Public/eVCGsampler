import_ui <- function(id) {

    ns <- NS(id)

    fluidRow(
        class = "full-screen-home border-top border-2 border-dark",
        column(
            width = 12,
            fixedPage(
                # ------ UPLOADS ----------------------------------------------
                div(
                    class = "card mt-4",
                    div(
                        class = "card-container",
                        h4("1. Data import", class = "vcg-section-title"),
                        p(
                            class = "vcg-section-hint",
                            "Upload the historical control POOL. The treated animals (TG) can be in",
                            "the same file - then it already holds a group variable with two values -",
                            "or in a second file, in which case the group variable
                            ",
                            tags$code(group_indicator), " is created for you.",
        
                            "Data must be provided in the standard wide (rectangular) format, with one row per animal. SENT-format is not supported."
                          
                        ),
                        fluidRow(
                            column(
                                width = 6,
                                fileInput(
                                    inputId = ns("file_pool"),
                                    label   = "POOL data, or POOL and TG combined (.csv, .txt, .xlsx, .xls)",
                                    accept  = c(".csv", ".txt", ".tsv", ".xlsx", ".xls")
                                )
                            ),
                            column(
                                width = 6,
                                fileInput(
                                    inputId = ns("file_tg"),
                                    label   = "Optional: TG data - covariates must match POOL file names.",
                                    accept  = c(".csv", ".txt", ".tsv", ".xlsx", ".xls")
                                )
                            )
                        ),
                        # ------ * Example data -------------------------
                        # One click into the workflow, for trying the app out
                        # and for demos. Handled by input$load_example in
                        # mod_import_server.R.
                        p(
                            class = "vcg-section-hint mb-0",
                            "No data at hand? ",
                            actionLink(
                                inputId = ns("load_example"),
                                label   = "Load the example data set",
                                class   = "link-primary fw-bold"
                            ),
                            " (simulated, for trying the app out)."
                        ),
                        uiOutput(ns("upload_msg")),
                        uiOutput(ns("tiles"))
                    )
                ),
                # ------ OVERVIEW ---------------------------------------------
                shinyjs::hidden(
                    div(
                        id = ns("overview_card"),
                        class = "card mt-4",
                        div(
                            class = "card-container",
                            tabsetPanel(
                                id = ns("overview_tabs"),
                                tabPanel(
                                    title = "Variable overview",
                                    div(
                                        class = "mt-3",
                                        p(
                                            class = "vcg-section-hint",
                                            "What each column looks like, and which role it can play in the",
                                            "balancing. A group variable needs exactly two distinct values, a",
                                            "stratum variable at most ", max_strata, ", and a covariate has to be",
                                            "numeric or a low-cardinality category."
                                        ),
                                        shinycssloaders::withSpinner(
                                            DT::dataTableOutput(ns("dt_vars"))
                                        )
                                    )
                                ),
                                tabPanel(
                                    title = "Combining check",
                                    div(
                                        class = "mt-3",
                                        p(
                                            class = "vcg-section-hint",
                                            "Only variables present in both files survive the merge. Anything",
                                            "listed as dropped is not available for balancing; a type mismatch",
                                            "usually means one file has text in a numeric column."
                                        ),
                                        uiOutput(ns("combine_msg")),
                                        DT::dataTableOutput(ns("dt_combine"))
                                    )
                                ),
                                tabPanel(
                                    title = "Data preview",
                                    div(
                                        class = "mt-3",
                                        p(
                                            class = "vcg-section-hint",
                                            "The data as the sampler will see it, group variable included."
                                        ),
                                        DT::dataTableOutput(ns("dt_data"))
                                    )
                                )
                            )
                        )
                    )
                ),
                p(
                    class = "text-end mt-4 mb-5",
                    shinyjs::disabled(
                        actionButton(
                            inputId = ns("continue"),
                            label = "Continue",
                            class = "btn-primary text-uppercase",
                            onclick = "$('#tabcontent').data('shiny-input-binding').setValue($('#tabcontent'), 'setup');"
                        )
                    )
                )
            )
        )
    )
}
