setup_ui <- function(id) {

    ns <- NS(id)

    fluidRow(
        class = "full-screen-home border-top border-2 border-dark",
        column(
            width = 12,
            fixedPage(
                fluidRow(
                    class = "mt-4",
                    # ------ SETTINGS -------------------------------------------
                    column(
                        width = 4,
                        div(
                            class = "card h-100",
                            div(
                                class = "card-container vcg-side",
                                h4("2. Variables and VCG size", class = "vcg-section-title"),
                                p(
                                    class = "vcg-section-hint",
                                    "The selection below is turned into the formula that eVCGsampler works with."
                                ),
                                uiOutput(ns("id_ui")),
                                uiOutput(ns("group_ui")),
                                uiOutput(ns("covariates_ui")),
                                uiOutput(ns("stratum_ui")),
                                textInput(
                                    inputId = ns("sizes"),
                                    label   = "VCG size (n) - one value, or one per stratum, comma-separated e.g. 7, 5",
                                    value   = "5"
                                ),
                                textInput(
                                    inputId = ns("weights"),
                                    label   = "Covariate importance weights - one per covariate, comma-separated e.g. 2, 1, 1 (optional)",
                                    value   = ""
                                )
                            )
                        )
                    ),
                    # ------ FORMULA + SIZES ------------------------------------
                    column(
                        width = 8,
                        div(
                            class = "card",
                            div(
                                class = "card-container",
                                h4("Formula", class = "vcg-section-title"),
                                p(
                                    class = "vcg-section-hint",
                                    "Covariates are added with ", tags$code("+"), " and strata are separated",
                                    "with ", tags$code("|"), ". Balancing runs within each stratum."
                                ),
                                uiOutput(ns("formula_box")),
                                uiOutput(ns("weights_msg"))
                            )
                        ),
                        div(
                            class = "card mt-4",
                            div(
                                class = "card-container",
                                h4("Sample size overview", class = "vcg-section-title"),
                                p(
                                    class = "vcg-section-hint",
                                    "Counts are after dropping rows with a missing value in any selected",
                                    "variable - exactly the data the sampler and the tests will use."
                                ),
                                uiOutput(ns("tiles")),
                                div(
                                    class = "mt-3",
                                    DT::dataTableOutput(ns("dt_sizes"))
                                ),
                                uiOutput(ns("size_msg"))
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
                            onclick = "$('#tabcontent').data('shiny-input-binding').setValue($('#tabcontent'), 'pretest');"
                        )
                    )
                )
            )
        )
    )
}
