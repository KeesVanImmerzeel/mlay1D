library(shinyMatrix)
library(magrittr)

nrw <- 11
ncl <- 3

kD  <-  c(35 * 30, 80 * 30, (55 / 2) * 0.075) %>% rep(nrw) %>% matrix(ncol = nrw) %>% t() # Transmissivity [m2/day]
c_  <- c(50, 400, 85 / 0.075) %>% rep(nrw) %>% matrix(ncol = nrw) %>% t() # Vertical Hydraulic resistance [day]
h   <- matrix(c(-1.10, -3.85, -1.20, -1.00, -0.80, -0.40, 0.00, 0.40, 0.80, 1.20, 1.60), nrw, 1) # Head on top of each section [m]
x   <- matrix(c(-1000, 1000, 3250, 4500, 5500, 6500, 7250, 8750, 9750, 10500), (nrw-1), 1) # Coordinates of intersection points [m]
Q   <- matrix(0, nrow = (nrw-1), ncol = ncl) # Matrix of nodal injections [m2/day]

grid_min <- -2500 # Minimum coordinate where values will be computed (m)
grid_max <- 11000 # Maximum coordinate where values will be computed (m)
nr_grid_points <- 1 + (grid_max-grid_min)/100  # Number of points where values will be computed (-)

shiny::tagList(
    fluidPage(
        titlePanel("mLay1D"),
        tabsetPanel(
            tabPanel(title="Control",
                     sidebarPanel(
                         width=3,
                         actionButton("go", "Refresh results & plots", class = "btn-warning"),
                         numericInput("ncl", "Nr. of aquifers", ncl, min = 1, max = 5),
                         numericInput("nrw", "Nr. of sections", nrw, min = 2, max = 5),
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
                         numericInput("grid_min", "Calc. grid min [m]", grid_min),
                         numericInput("grid_max", "Calc. grid max [m]", grid_max),
                         numericInput("nr_grid_points", "Nr. of calc. grid points (-)", nr_grid_points, min=1, max=10000)
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
            tabPanel(title="Plots",
                     sidebarPanel(
                         width="100%",
                         plotOutput(outputId = "phi_plot"),
                         plotOutput(outputId = "q_plot"),
                         plotOutput(outputId = "s_plot")
                    )
            ),
            tabPanel(title = "Documentation", 
                     shiny::includeHTML("readme.html")
                     )
        )  
    ) #shiny::fluidPage
)

