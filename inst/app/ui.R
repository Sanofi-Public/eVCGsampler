

ui <- fluidPage(
  fluidRow(
    column(1,
           tags$img(src = "logo.png", height = "100px"),
    ),
    column(8,
           titlePanel("Energy distance based VCG sampler"),
    )
  ),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload POOL or combined POOL+TG data here, in excel or csv format.", accept = c(".csv", ".xlsx")),
      conditionalPanel(condition = "output.fileUploaded",
                       fileInput("file2", "Optional: Upload the treated (TG) data here. The app will create POOL_vs_TG variable. Covariates must have the same names in both data files.", accept = c(".csv", ".xlsx")),
                       uiOutput("outcome_ui"),
                       uiOutput("covariates_ui"),
                       uiOutput("stratum_ui"),
                       textInput("n1",  "VCG size (n), for different size in stratum comma-separated e.g. (5, 10, 6),",  value = 5),
                       textInput("c_w", "Covariate importance weights: e.g. (1, 2, 1),",  value = NULL),
                       checkboxInput("rand", 'Use random subset sampling? Otherwise, the best possible will be used.', value = FALSE),
                       downloadButton("downloadData", "Download Results (data set)")
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Results",
                 plotOutput("Plot1", width = "600px", height = "600px"),
                 h4("Median (MAD):"),
                 tableOutput("s_table")
        ),
        tabPanel("Distance tests",
                 h4("TG vs POOL permutation test"),
                 plotOutput("Plot2", width = "600px", height = "400px"),
                 h4("TG vs VCG permutation test"),
                 plotOutput("Plot3", width = "600px", height = "400px")
        ),
        tabPanel("Plot covariates",
                 uiOutput("plots_ui")

        ),
        tabPanel("Optimal VCG size",
                 h4("Recommendation for optimal VCG size (purely exploratory approach)"),
                 plotOutput("plots_n")

        )
      )
    )
  )
)

