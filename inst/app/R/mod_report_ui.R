report_ui <- function(id) {

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
                        width = 5,
                        div(
                            class = "card h-100",
                            div(
                                class = "card-container vcg-side",
                                h4("5. Report", class = "vcg-section-title"),
                                p(
                                    class = "vcg-section-hint",
                                    "One self-contained HTML file with the data checks, the settings,",
                                    "the pre-tests and the sampling result - everything needed to",
                                    "document how this VCG was drawn."
                                ),
                                textInput(
                                    inputId = ns("study"),
                                    label   = "Study title (optional)",
                                    value   = ""
                                ),
                                textInput(
                                    inputId = ns("author"),
                                    label   = "Prepared by (optional)",
                                    value   = ""
                                ),
                                textAreaInput(
                                    inputId = ns("comment"),
                                    label   = "Comment (optional)",
                                    value   = "",
                                    rows    = 3
                                ),
                                checkboxInput(
                                    inputId = ns("inc_covplots"),
                                    label   = "Include the per-covariate distribution plots",
                                    value   = TRUE
                                ),
                                checkboxInput(
                                    inputId = ns("inc_outlier"),
                                    label   = "Include the outlier study test",
                                    value   = TRUE
                                ),
                                checkboxInput(
                                    inputId = ns("inc_data"),
                                    label   = "Include the sampled VCG animals as a table",
                                    value   = TRUE
                                ),
                                tags$hr(class = "my-4"),
                                shinyjs::disabled(
                                    downloadButton(
                                        outputId = ns("download_html"),
                                        label    = "Download HTML report",
                                        class    = "btn-primary"
                                    )
                                ),
                                div(
                                    class = "mt-2",
                                    shinyjs::disabled(
                                        downloadButton(
                                            outputId = ns("download_xlsx"),
                                            label    = "Download all data (Excel)",
                                            class    = "btn-primary"
                                        )
                                    )
                                ),
                                p(
                                    class = "vcg-section-hint mt-2 mb-0",
                                    "The complete data set - the selected VCG, the POOL animals that",
                                    "were not selected and the TG animals, with their energy weights.",
                                    "The Sampling step exports the selected animals on their own."
                                )
                            )
                        )
                    ),
                    # ------ WHAT GOES IN --------------------------------------
                    column(
                        width = 7,
                        div(
                            class = "card",
                            div(
                                class = "card-container",
                                h4("Report content", class = "vcg-section-title"),
                                p(
                                    class = "vcg-section-hint",
                                    "Steps you skipped are simply left out of the report."
                                ),
                                uiOutput(ns("content_table"))
                            )
                        )
                    )
                ),
                p(
                    class = "text-end mt-4 mb-5",
                    actionButton(
                        inputId = ns("restart"),
                        label = "Start over",
                        class = "btn-primary text-uppercase",
                        onclick = "$('#tabcontent').data('shiny-input-binding').setValue($('#tabcontent'), 'home');"
                    )
                )
            )
        )
    )
}
