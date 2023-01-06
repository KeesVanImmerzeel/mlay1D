library(shinyMatrix)
library(magrittr)

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
            tabPanel(title="Control",
                     sidebarPanel(
                         width=4,
                         actionButton("go", "Refresh results & plots", class = "btn-warning"),
                         numericInput("ncl", "Nr. of aquifers", ncl, min = 1, max = 5),
                         numericInput("nrw", "Nr. of sections", nrw, min = 2, max = 26),
                         matrixInput(
                             inputId = "x",
                             label = "Intersection points [m]",
                             value = x,
                             class = "numeric",
                             rows = list(
                                 names = FALSE
                            ), 
                            cols=list(
                                names= FALSE
                            )
                         ),
                         numericInput("f_", "Divide sections in f-parts [-]", f_, step=1, min = 1, max = 100),
                         numericInput("grid_min", "Presentation grid min [m]", grid_min),
                         numericInput("grid_max", "Presentation grid max [m]", grid_max),
                         numericInput("nr_grid_points", "Nr. of presentation points (-)", nr_grid_points, min=1, max=10000),
                         shinyFiles::shinySaveButton('download', 'Download', 'Save input data.', 
                                                     multiple=FALSE, filename="mlay1D.rds", filetype=list(picture=c('rds')), icon = icon("download")),
                         shinyFiles::shinyFilesButton('upload', 'Upload', 'Readinput data.', 
                                                      multiple=FALSE, filename="mlay1D.rds", filetype=list(picture=c('rds')), icon = icon("upload"))
                         )),
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
            tabPanel(title="Results",
                     sidebarPanel(
                         width="100%",
                         dataTableOutput("matrix"))
                     ),
            tabPanel(title = "Plots",
                     sidebarLayout(
                           sidebarPanel(
                                 width = "4",
                                 shinyFiles::shinySaveButton(
                                       'dwnld_phi_plot',
                                       'Download Head plot',
                                       'Save Head plot.',
                                       multiple =
                                             FALSE,
                                       filename = "Head.png",
                                       filetype = list(picture = c('png')),
                                       icon = icon("download")
                                 )
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
            tabPanel(title="Export",
                     sidebarPanel(
                           width="3",
                           shinyFiles::shinySaveButton('export', 'Export results', 'Save results to file.', 
                                                       multiple=FALSE, filename="mlay1D", filetype=list(picture=c('csv')), icon = icon("download")))
                     
            ),            
            tabPanel(title = "Documentation", 
                     shiny::includeHTML("readme.html")
                     )
        )  
    ) #shiny::fluidPage
)

