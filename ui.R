library(shinyMatrix)
library(magrittr)

# options set to make datatables work correctly 
options(shiny.legacy.datatable = TRUE)

nrw <- 11
ncl <- 3

kD  <-  c(35 * 30, 80 * 30, (55 / 2) * 0.075) %>% rep(nrw) %>% matrix(ncol = nrw) %>% t() # Transmissivity [m2/day]
c_  <- c(50, 400, 85 / 0.075) %>% rep(nrw) %>% matrix(ncol = nrw) %>% t() # Vertical Hydraulic resistance [day]
h   <- matrix(c(-1.10, -3.85, -1.20, -1.00, -0.80, -0.40, 0.00, 0.40, 0.80, 1.20, 1.60), nrw, 1) # Head on top of each section [m]
x   <- matrix(c(-1000, 1000, 3250, 4500, 5500, 6500, 7250, 8750, 9750, 10500), (nrw-1), 1) # Coordinates of intersection points [m]
Q   <- matrix(0, nrow = (nrw-1), ncol = ncl) # Matrix of nodal injections [m2/day]
f_ <- 1 # Number of sub-sections per intersection interval to create (f=1: no sub-sections. [-] integer)

grid_min <- -2500 # Minimum coordinate where values will be computed (m)
grid_max <- 11000 # Maximum coordinate where values will be computed (m)
nr_grid_points <- 1 + (grid_max-grid_min)/100  # Number of points where values will be computed (-)

shiny::tagList(
    fluidPage(
        titlePanel("mLay1D"),
        tabsetPanel(
            id="all_tabs",
            tabPanel(
                  title = "Control",
                  sidebarPanel(
                        width = 4,
                        actionButton("go", "Refresh results & plots", class = "btn-warning"),
                        numericInput("ncl", "Nr. of aquifers", ncl, min = 1, max = 5),
                        numericInput("nrw", "Nr. of sections", nrw, min = 2, max = 26),
                        matrixInput(
                              inputId = "x",
                              label = "Intersection points [m]",
                              value = x,
                              class = "numeric",
                              rows = list(names = FALSE),
                              cols = list(names = FALSE)
                        ),
                        numericInput(
                              "f_",
                              "Divide sections in f-parts [-]",
                              f_,
                              step = 1,
                              min = 1,
                              max = 100
                        ),
                        numericInput("grid_min", "Presentation grid min [m]", grid_min),
                        numericInput("grid_max", "Presentation grid max [m]", grid_max),
                        numericInput(
                              "nr_grid_points",
                              "Nr. of presentation points (-)",
                              nr_grid_points,
                              min = 1,
                              max = 10000
                        ),
                        textInput("download_filename", "Name for the file with input data", value = "mlay1D.rds"),
                        downloadButton("download", "Download input data."),
                        br(),
                        br(),
                        fileInput("upload", "Upload input data", accept =
                                        ".rds")
                  )
            ),
            tabPanel(title = "Transmissivity [m2/day]",
                     sidebarPanel(
                         width="100%",
                         matrixInput(
                             inputId = "kD",
                             value = kD,
                             class = "numeric",
                             cols = list(
                                 names = TRUE
                             ),
                             rows = list(
                                 names = TRUE
                             )
                         ))),
            tabPanel(title = "Vertical Hydraulic resistance [day]",
                     sidebarPanel(
                         width="100%",
                         matrixInput(
                             inputId = "c_",
                             value = c_,
                             class = "numeric",
                             cols = list(
                                 names = TRUE
                             ),
                             rows = list(
                                 names = TRUE
                             )
                         ))),
            tabPanel(title = "Head on top of each section [m]",
                     sidebarPanel(
                         width=4,
                         matrixInput(
                             inputId = "h",
                             value = h,
                             class = "numeric",
                             rows = list(
                                 names = TRUE
                             )
                         ))),
            tabPanel(title = "Nodal injections [m2/day]",
                     sidebarPanel(
                         width="100%",
                         matrixInput(
                             inputId = "Q",
                             value = Q,
                             class = "numeric",
                             cols = list(
                                 names = TRUE
                             ),
                             rows = list(
                                 names = TRUE
                             )
                         ))),
            
            tabPanel(title = "Measurements",
                     # Sidebar layout with options to upload measurements
                     sidebarLayout(
                           # Sidebar panel for inputs ----
                           sidebarPanel(
                                 # Input: Select a file ----
                                 fileInput(
                                       "fname_measurements",
                                       "Choose CSV File with measurements",
                                       accept = c("text/csv",
                                                  "text/comma-separated-values,text/plain",
                                                  ".csv")
                                 ),
                                 # Horizontal line ----
                                 tags$hr(),
                                 # Input: Select separator ----
                                 radioButtons(
                                       "sep",
                                       "Separator",
                                       choices = c(
                                             Comma = ",",
                                             Semicolon = ";",
                                             Tab = "\t"
                                       ),
                                       selected = ";"
                                 ),
                                 
                                 # Input: Select quotes ----
                                 radioButtons(
                                       "quote",
                                       "Quote",
                                       choices = c(
                                             None = "",
                                             "Double Quote" = '"',
                                             "Single Quote" = "'"
                                       ),
                                       selected = ''
                                 )
                           ),
                           # Main panel for displaying outputs ----
                           mainPanel(# Output: Data file ----
                                     tableOutput("contents"))
                     )),                           
                        
            tabPanel(title="Results",
                     sidebarLayout(
                     sidebarPanel(
                         width="3",
                         textInput("export_filename", "Name for the file with results", value = "mlay1D.csv"),
                         downloadButton("export", "Download results")
                         ),
                     mainPanel(# Output: Data file ----
                               dataTableOutput("matrix"))
                     )),
            tabPanel(title = "Plots",
                     sidebarLayout(
                           sidebarPanel(
                                 width = "4",
                                 checkboxInput("Labels", "Label measurements", value = TRUE),
                                 br(), br(),
                  textInput("phi_plot_filename", "Name for the Head plot file", value = "Head.png"),
                  downloadButton("dwnld_phi_plot", "Download Head plot"),
                                 br(),
                                 br(),
                  textInput("latflx_filename", "Name for the Lateral Flux plot file", value = "Lateral Flux.png"),
                  downloadButton("dwnld_latflx_plot", "Download Lateral Flux plot"),
                                 br(),
                                 br(),
                  textInput("seepage_filename", "Name for the Seepage plot file", value = "Seepage.png"),
                  downloadButton("dwnld_seepage_plot", "Download Seepage plot")
                           ),
                           mainPanel(
                                 verbatimTextOutput("phi_info"),
                                 plotOutput(outputId = "phi_plot", click = "plot_click_phi"),
                                 br(),
                                 br(),
                                 verbatimTextOutput("q_info"),
                                 plotOutput(outputId = "q_plot", click = "plot_click_q"),
                                 br(),
                                 br(),
                                 verbatimTextOutput("s_info"),
                                 plotOutput(outputId = "s_plot", click = "plot_click_s")
                           )
                     )),
            tabPanel(title = "Documentation", 
                     shiny::includeHTML("readme.html")
                     )
        )  
    ) #shiny::fluidPage
)

