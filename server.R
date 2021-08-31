library(shinyMatrix)
library(expm)
library(pracma)
library(ggplot2)
library(reshape2)
library(dplyr)
library(magrittr)

# Constants
min_nr_sections <- 2; max_nr_sections <- 25
min_nr_aquifers <- 1; max_nr_aquifers <- 5
min_kD <- 1;         max_kD <- 10000
min_Q <- -10^6;       max_Q <- 10^6
min_c <- 1;          max_c <- 10000
min_h <- -1000;       max_h <- 1000


# Create a new matrix based on the values in matrices x and y.
# The resulting matrix has the dimensions of the matrix y.
# Values of the matrix x are copied into the upper left corner of the resulting matrix.
# If the dimensions of matrix y is smaller then the dimensions of matrix x, 
# the values not fitting are lost.
resize_matrix <- function(x, y) {
  if (nrow(y) > nrow(x)) {
    if (ncol(y) > ncol(x)) {
      y[1:nrow(x), 1:ncol(x)] <- x
    } else {
      y[1:nrow(x), 1:ncol(y)] <- x[1:nrow(x), 1:ncol(y)]
    }
  } else {
    if (ncol(y) > ncol(x)) {
      y[1:nrow(y), 1:ncol(x)] <- x[1:nrow(y), 1:ncol(x)]
    } else {
      y[1:nrow(y), 1:ncol(y)] <- x[1:nrow(y), 1:ncol(y)]
    }
  }
  return(y)
}

# Limit values in a vector x.
limit_vector <- function(x, min_value, max_value) {
  sapply(x, function(y)
    min(max(y, min_value), max_value))
}

# Limit values in matrix x.
limit_matrix <- function(x, min_value, max_value) {
  apply(x, c(1,2), function(y) min(max(y,min_value),max_value))
}

#' Allokate parameters (kD, c_, h) on refined intersection point vector.
#' @param par_org Original parameter (kD, c_, h) 
#' @param x_org  Original vector of coordinates of intersection points
#' @param x Vector of coordinates of intersection points (except +/-inf, [L], nNod by 1)
#' @return Refined model parameter (matrix).
.refine_parm <- function( par_org, x_org, x ){
      cbind(par_org[,1], 
            par_org %>% apply(1, function(y)
                  yorg <-
                        y[2:length(y)] %>% as.array() %>% stats::approx(
                              x = x_org,
                              y = .,
                              xout = x,
                              method = "constant",
                              f = 0
                        )) %>% lapply(function(a) a$y) %>% do.call(rbind, .))
}

#' Allokate all model parameters on refined intersection point vector specified by 'f' and 'x'.
#' @inheritParams .refine_parm      
#' @param f  Number of sub-sections per intersection interval to create (f=1: no sub-sections. [-] integer)
#' @param kD Matrix of transmissivity values  ([L2/T], nLay by nSec)
#' @param c_ Matrix of vertical resistance values of overlaying aquitards ([T], nLay by nSec) 
#' @param Q  Matrix of nodal injections ([L2/T], nLay by nSec-1) 
#' @param h  Vector of given heads on top of each sections  ([L], 1 by nSec)
#' @return All model parameters (list).
.refine_model <- function(f, kD, c_, Q, h, x) {
      res <- list()
      x_org <- x
      dx <- diff(x_org) / f
      midpoints <-
            1:(f - 1) %>% as.array() %>% apply(1, function(a)
                  x_org[1:(length(x_org) - 1)] + a * dx) %>% as.vector()
      res$x <- c(x_org, midpoints) %>% sort()
      res$kD <- kD %>% .refine_parm(x_org, res$x)
      res$c_ <- c_ %>% .refine_parm(x_org, res$x)
      res$h <-
            h %>% as.matrix() %>% t() %>% .refine_parm(x_org, res$x) %>% as.vector()
      res$Q <- matrix(0, nrow = nrow(Q), ncol = length(res$x))
      res$Q[, match(x_org, res$x)] <- Q
      return(res)
}

#' Calculate the head and fluxes of an analytic steady-state, one-dimensional, 
#' leaky multi-aquifer model using n connected sections. 
#' Left and right most sections run to -inf and +inf resp.
#' @inheritParams .refine_model     
#' @param X  Vector of points where values will be computed ([L] =npoints)
#' @return   Matrix with calculated heads ([L] rows (1-nLay)); lateral fluxes ([L2/T] rows (nLay+1 - 2xnLay)); seepage ([L/T] rows (2xnLay+1 - 3xnLay))
solve_mlay1d <- function(kD, c_, Q, h, x, X, f=1) {
      if (f > 1) {
            # Create sub-sections per intersection interval to increase numerical stability
            #print(paste("Refine modelgrid f=",as.character(f)))
            res <- .refine_model(f, kD, c_, Q, h, x)
            kD <- res$kD
            c_ <- res$c_
            Q <- res$Q
            h <- res$h
            x <- res$x
      }
      
      Q %<>% cbind(0, ., 0)
      x %<>% c(-Inf, ., Inf)
      nSec <- length(h)
      nNod <- nSec - 1
      nLay <- nrow(c_)
      
      Nx <- x %>% length()
      H <- h %>% rep(nLay) %>% matrix(nrow = nLay, byrow = TRUE)
      
      #  are used to compute relative coordinates within sections
      xMidSec <- .5 * (x[-length(x)] + x[-1])
      xMidSec[1] <- x[2]
      xMidSec[length(xMidSec)] <- x[length(x) - 1]
      
      withProgress(message="System matrices for all sections.", min=1, max=nSec, value=0, {
      
      # System matrices for all sections.
      A <- array(dim = c(nLay, nLay, nSec))
      for (iSec in 1:nSec) {
            a <- 1 / (kD[, iSec] * c_[, iSec])
            if (nLay > 1) {
                  b <- 1 / (kD[, iSec] * c(c_[2:nLay, iSec], Inf))
                  A[, , iSec] <-
                        diag(a + b) - pracma::Diag(c(a[2:nLay]), -1) - pracma::Diag(c(b[1:(nLay - 1)]), 1)
            } else {
                  b <- 0 %>% as.matrix()
                  A[, , iSec] <- diag(a + b)
            }
            incProgress(1)
      }
      
      })
      
      C <- matrix(0, nrow = nLay * (2 * (Nx - 2)),
                  ncol = nLay * (2 * (Nx - 2) + 2))
      R <-  matrix(0, nrow = nLay * (2 * (Nx - 2)), ncol = 1)
      
      withProgress(message="Creating matrices C & R", min=1, max=nNod, value=0, {
      
      for (i in 1:nNod) {
            j <- i + 1
            ii <- 2 * nLay * (i - 1) + 1
            jj <- ii + nLay - 1
            C[ii:jj, ii:jj] <-
                  +expm::expm(-(x[j] - xMidSec[i]) * expm::sqrtm(as.matrix(A[, , i])))
            C[ii:jj, (ii + nLay):(jj + nLay)] <-
                  +expm::expm(+(x[j] - xMidSec[i]) * expm::sqrtm(as.matrix(A[, , i])))
            C[ii:jj, (ii + 2 * nLay):(jj + 2 * nLay)] <-
                  -expm::expm(-(x[j] - xMidSec[j]) * expm::sqrtm(as.matrix(A[, , j])))
            
            C[ii:jj, (ii + 3 * nLay):(jj + 3 * nLay)] <-
                  -expm::expm(+(x[j] - xMidSec[j]) * expm::sqrtm(as.matrix(A[, , j])))
            R[ii:jj] <- -H[, i] + H[, j]
            
            di <- kD[, i]
            dj <- kD[, j]
            if (nLay == 1) {
                  di %<>% as.matrix()
                  dj %<>% as.matrix()
            }
            di %<>% diag()
            dj %<>% diag()
            
            exp1 <- (x[j] - xMidSec[i]) * expm::sqrtm(as.matrix(A[, , i]))
            exp2 <- (x[j] - xMidSec[j]) * expm::sqrtm(as.matrix(A[, , j]))
            #crit <- range(c(range(exp1),c(range(exp2)))) %>% abs() %>% max()
            #if (crit>20) {
            #      m <- solve_mlay1d(kD, c_, Q[,2:(ncol(Q)-1)], h, x[2:(length(x)-1)], X, f=3) 
            #      return( m )
            #}            
            
            C[(ii + nLay):(jj + nLay), (ii:jj)] <-
                  -di %*% expm::sqrtm(as.matrix(A[, , i])) %*% expm::expm(-exp1)
            C[(ii + nLay):(jj + nLay), (ii + nLay):(jj + nLay)] <-
                  +di %*% expm::sqrtm(as.matrix(A[, , i])) %*% expm::expm(exp1)
            
            C[(ii + nLay):(jj + nLay), (ii + 2 * nLay):(jj + 2 * nLay)] <-
                  +dj %*% expm::sqrtm(as.matrix(A[, , j])) %*% expm::expm(-exp2)
            C[(ii + nLay):(jj + nLay), (ii + 3 * nLay):(jj + 3 * nLay)] <-
                  -dj %*% expm::sqrtm(as.matrix(A[, , j])) %*% expm::expm(exp2)
            R[(ii + nLay):(jj + nLay)] <- Q[, j]
            incProgress(1)
      }
            
      }) # withProgress
      
      
      withProgress(message="Solving linear system of equations.", {
      # Mimic mldivide
      COEF <- solve(C[, (nLay + 1):(ncol(C) - nLay)], R) %>%
            c(rep(0, nLay), ., rep(0, nLay))
      }) # withProgress
      
      
      withProgress(message="Creating matrix with heads, lateral fluxes and seepage.", min=1, max=length(X), value=0, {
            
      phi <- matrix(0, nrow = nLay, ncol = length(X))
      q <- phi
      s <- q
      for (i in 1:length(X)) {
            iSec <- which(X[i] > x[1:(length(x) - 1)] & X[i] <= x[2:length(x)])
            k <- 2 * nLay * (iSec - 1) + 1
            l <- k + nLay - 1
            sqm <- expm::sqrtm(as.matrix(A[, , iSec]))
            d <- X[i] - xMidSec[iSec]
            C1 <-
                  expm::expm(-d * sqm) %*% COEF[k:l]
            C2 <-
                  expm::expm(d * sqm) %*% COEF[(k + nLay):(l + nLay)]
            C3 <- sqm %*% C1
            C4 <- sqm %*% C2
            phi[, i] <- C1 + C2 + H[, iSec]
            dia <- kD[, iSec]
            if (nLay == 1) {
                  dia %<>% as.matrix()
            }
            dia %<>% diag()
            q[, i] <- dia %*% (C3 - C4)
            sNet <- dia %*% sqm %*% (C3 + C4)
            s[nLay, i] <- sNet[nLay]
            if (nLay > 1) {
                  for (iLay in (nLay - 1):1) {
                        s[iLay, i] <- sNet[iLay] + s[(iLay + 1), i]
                  }
            }
      }
      
      })
      return(rbind(phi, q, s))
}

#' Plot results of function solve_mlay1d
#' @inheritParams solve_mlay1d
#' @param m Matrix with calculated heads; lateral fluxes; seepage
#' @layers number(s) of layers to plot (numeric)
#' @param ptype Type of plot to create ("phi"=default, "q", "s")
#' @param xvlines x-coordinates of (optional) vertical lines in plot (numeric vector)
#' return ggplot2 object
plot_mlay1d <- function(m, X, layers = 1, ptype = "phi", xvlines=NULL) {
      sel_names <- paste0(ptype, layers)
      m %<>% t() %>% as.data.frame()
      nlay <- ncol(m) / 3
      names(m) <-
            c(paste0("phi", 1:nlay),
              paste0("q", 1:nlay),
              paste0("s", 1:nlay))
      sel_names <- paste0(ptype, layers)
      m %<>% dplyr::select(all_of(sel_names))
      names(m) <- layers %>% as.character()
      m$X <- X
      m %<>% reshape2::melt(id.vars = "X", variable.name = "Aquifer")
      
      if (ptype == "phi") {
            ylab <- "Head (m+ref)"
            title <- "Head"
      } else if (ptype == "q") {
            ylab <- "Lateral flux (m2/d)"
            title <- "Lateral flux"
      } else {
            ylab <- "Seepage (m/d)"
            title <- "Seepage"
      }
      myplot <-
            ggplot(data = m, aes(
                  x = X,
                  y = value,
                  colour = Aquifer
            )) +
            geom_line(size=1) +
            labs(
                  colour = "Aquifer",
                  x = "X (m)",
                  y = ylab,
                  title = title
            ) +
            theme(plot.title = element_text(hjust = 0.5)) +
            geom_vline(xintercept = xvlines, linetype="dotted", 
                       color = "blue")
      return(myplot)
}

function(input, output, session) {
  nrw <-
    reactive({
      trunc(min(max(input$nrw, min_nr_sections), max_nr_sections))
    })
  ncl <-
    reactive({
      trunc(min(max(input$ncl, min_nr_aquifers), max_nr_aquifers))
    })
  
  observeEvent(input$Q, {
    Q <- limit_matrix(input$Q, min_Q, max_Q)
    updateMatrixInput(session, "Q", Q)
  })
  
  observeEvent(input$kD, {
    kD <- limit_matrix(input$kD, min_kD, max_kD)
    updateMatrixInput(session, "kD", kD)
  })
  
  observeEvent(input$c_, {
    c_ <- limit_matrix(input$c_, min_c, max_c)
    updateMatrixInput(session, "c_", c_)
  })
  
  observeEvent(input$h, {
    h <- limit_matrix(input$h, min_h, max_h)
    updateMatrixInput(session, "h", h)
  })
  
  observeEvent(input$x, {
    x <- sort(input$x) %>% as.matrix()
    updateMatrixInput(session, "x", x)
  })
  
  observeEvent(c(input$grid_min, input$grid_max), {
    if (input$grid_min > input$grid_max) {
      hlp <- input$grid_max
      updateNumericInput(session, "grid_max", value = input$grid_min)
      updateNumericInput(session, "grid_min", value = hlp)
    }
   
  })

  observeEvent(c(input$ncl, input$nrw), {
    updateNumericInput(session, "nrw", value = nrw())
    updateNumericInput(session, "ncl", value = ncl())
    
    # Update kD MatrixInput
    kD_resized <- resize_matrix(input$kD, matrix(0, nrw(), ncl()))
    colnames(kD_resized) <- paste0("kD", 1:ncl())
    rownames(kD_resized) <- paste("Section", 1:nrw())
    updateMatrixInput(session, "kD", kD_resized)
    
    # Update Q MatrixInput  
    Q_resized <- resize_matrix(input$Q, matrix(0, (nrw()-1), ncl()))
    colnames(Q_resized) <- paste0("Q", 1:ncl())
    rownames(Q_resized) <- paste("Section", 1:(nrw()-1))
    updateMatrixInput(session, "Q", Q_resized)    
    
    # Update c_ MatrixInput
    c_resized <- resize_matrix(input$c_, matrix(0, nrw(), ncl()))
    colnames(c_resized) <- paste0("c", 0:(ncl() - 1))
    rownames(c_resized) <- paste("Section", 1:nrw())
    updateMatrixInput(session, "c_", c_resized)
    
    # Update h MatrixInput
    h_resized <- resize_matrix(input$h, matrix(0, nrw(), 1))
    colnames(h_resized) <- "h" 
    rownames(h_resized) <- paste("Section", 1:nrw())
    updateMatrixInput(session, "h", h_resized)
    
    # Update x MatrixInput
    x_resized <- resize_matrix(input$x, matrix(0, (nrw() - 1), 1))
    updateMatrixInput(session, "x", x_resized)
    
  })

  X <-
    reactive({
      seq(input$grid_min, input$grid_max, length.out=input$nr_grid_points)
    })  
  
  m_ <- eventReactive(input$go, {
    h <- input$h %>% as.vector()
    x <- input$x %>% as.vector()
    X <- X() %>% as.vector()
    f <- input$f_ %>% as.vector()
    res <- solve_mlay1d(t(input$kD), t(input$c_), t(input$Q), h, x, X, f=f )
    return(res)
  })
  
  output$matrix <- renderDataTable({
    m <- m_()
    m %<>% t() %>% as.data.frame()
    nlay <- ncol(m) / 3
    names(m) <-
      c(paste0("phi", 1:nlay),
        paste0("q", 1:nlay),
        paste0("s", 1:nlay, "(x1000)"))
    m$X <- X()
    m[,1:nlay] %<>% round(.,2)
    m[,(2*nlay+1):(3*nlay)] <- 1000* m[,(2*nlay+1):(3*nlay)]
    m[,(nlay+1):(2*nlay)] %<>% round(.,4)
    m[,(2*nlay+1):(3*nlay)] %<>% round(.,1)
    m %<>% dplyr::relocate(X, .before = phi1)
    m
  })
  
  output$phi_plot <- renderPlot({
    m <- m_()
    nlay <- nrow(m) / 3
    x <- input$x %>% as.vector()
    plot_mlay1d(m, X(), layers = 1:nlay, ptype = "phi", xvlines=x)
  })
  
  output$q_plot <- renderPlot({
    m <- m_()
    nlay <- nrow(m) / 3
    x <- input$x %>% as.vector()
    plot_mlay1d(m, X(), layers = 1:nlay, ptype = "q", xvlines=x)
  })
  
  output$s_plot <- renderPlot({
    m <- m_()
    nlay <- nrow(m) / 3
    x <- input$x %>% as.vector()
    plot_mlay1d(m, X(), layers = 1:nlay, ptype = "s", xvlines=x)
  })
  
}