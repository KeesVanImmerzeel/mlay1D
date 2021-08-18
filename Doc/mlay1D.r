###_Install packages_###
packages.loading <-
      c( "expm",
         "installr",
         "rstudioapi",
         "pracma",
         "reshape2",
         "dplyr",
         "magrittr",
         "ggplot2"
      )
new.packages <-
      packages.loading[!(packages.loading %in% installed.packages()[, "Package"])]
if (length(new.packages))
      install.packages(new.packages, dependencies = TRUE)
lapply(packages.loading, require, character.only = TRUE)
installr::updateR() #If packages can't be installed update Rstudio

### Set work directory_###
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

### Functions ### 
#' Calculate the head and fluxes of an analytic steady-state, one-dimensional, 
#' leaky multi-aquifer model using n connected sections. 
#' Left and right most sections run to -inf and +inf resp.
#' @param kD Matrix of transmissivity values  ([L2/T], nLay by nSec)
#' @param c_ Matrix of vertical resistance values of overlaying aquitards ([T], nLay by nSec) 
#' @param Q  Matrix of nodal injections ([L2/T], nLay by nSec-1) 
#' @param h  Vector of given heads on top of each sections  ([L], 1 by nSec)
#' @param x  Vector of coordinates of intersection points (except +/-inf, [L], nNod by 1)
#' @param X  Vector of points where values will be computed ([L] =npoints)
#' @return   Matrix with calculated heads ([L] rows (1-nLay)); lateral fluxes ([L2/T] rows (nLay+1 - 2xnLay)); seepage ([L/T] rows (2xnLay+1 - 3xnLay))
solve_mlay1d <- function(kD, c_, Q, h, x, X) {
      nSec <- length(h)
      nNod <- nSec - 1
      nLay <- nrow(c_)
      
      Q %<>% cbind(0, ., 0)
      x %<>% c(-Inf, ., Inf)
      Nx <- x %>% length()
      H <- h %>% rep(nLay) %>% matrix(nrow = nLay, byrow = TRUE)
      
      #  are used to compute relative coordinates within sections
      xMidSec <- .5 * (x[-length(x)] + x[-1])
      xMidSec[1] <- x[2]
      xMidSec[length(xMidSec)] <- x[length(x) - 1]
      
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
      }
      
      C <- matrix(0, nrow = nLay * (2 * (Nx - 2)),
                  ncol = nLay * (2 * (Nx - 2) + 2))
      R <-  matrix(0, nrow = nLay * (2 * (Nx - 2)), ncol = 1)
      
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
            
            C[(ii + nLay):(jj + nLay), (ii:jj)] <-
                  -di %*% expm::sqrtm(as.matrix(A[, , i])) %*% expm::expm(-(x[j] - xMidSec[i]) * expm::sqrtm(as.matrix(A[, , i])))
            C[(ii + nLay):(jj + nLay), (ii + nLay):(jj + nLay)] <-
                  +di %*% expm::sqrtm(as.matrix(A[, , i])) %*% expm::expm(+(x[j] - xMidSec[i]) * expm::sqrtm(as.matrix(A[, , i])))
            
            C[(ii + nLay):(jj + nLay), (ii + 2 * nLay):(jj + 2 * nLay)] <-
                  +dj %*% expm::sqrtm(as.matrix(A[, , j])) %*% expm::expm(-(x[j] - xMidSec[j]) * expm::sqrtm(as.matrix(A[, , j])))
            C[(ii + nLay):(jj + nLay), (ii + 3 * nLay):(jj + 3 * nLay)] <-
                  -dj %*% expm::sqrtm(as.matrix(A[, , j])) %*% expm::expm(+(x[j] - xMidSec[j]) * expm::sqrtm(as.matrix(A[, , j])))
            R[(ii + nLay):(jj + nLay)] <- Q[, j]
      }
      
      # Mimic mldivide
      COEF <- solve(C[, (nLay + 1):(ncol(C) - nLay)], R) %>%
            c(rep(0, nLay), ., rep(0, nLay))
      
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
      return(rbind(phi, q, s))
}

#' Plot results of function solve_mlay1d
#' @inheritParams solve_mlay1d
#' @param m Matrix with calculated heads; lateral fluxes; seepage
#' @layers number(s) of layers to plot (numeric)
#' @param ptype Type of plot to create ("phi"=default, "q", "s")
#' return ggplot2 object
plot_mlay1d <- function(m, X, layers = 1, ptype = "phi") {
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
            geom_line() +
            labs(
                  colour = "Aquifer",
                  x = "X (m)",
                  y = ylab,
                  title = title
            ) +
            theme(plot.title = element_text(hjust = 0.5))
      return(myplot)
}

# Example 1
nLay <- 3
h <-
      c(-1.10,
        -3.85,
        -1.20,
        -1.00,
        -0.80,
        -0.40,
        -0.00,
        0.40,
        0.80,
        1.20,
        1.60)
nSec <- length(h)
kD <- matrix(rep(c(35 * 30, 80 * 30, (55 / 2) * 0.075), nSec), nrow = nLay)
c_ <- matrix(rep(c(50, 400, 85 / 0.075), nSec), nrow = nLay)
Q <- matrix(rep(0, nLay * (nSec - 1)), nrow = nLay)
x <- c(-1000, 1000, 3250, 4500, 5500, 6500, 7250, 8750, 9750, 10500)
X <- seq(-2500, 11000, 10)
m <- solve_mlay1d(kD, c_, Q, h, x, X)
m %>% plot_mlay1d(X, layers = c(1:nLay))
m %>% plot_mlay1d(X, layers = c(1:nLay), ptype = "q")
m %>% plot_mlay1d(X, layers = c(1:nLay), ptype = "s")

### Example 2
nLay <- 1
nSec <- 4
kD <- matrix(rep(2000, nSec), nrow = nLay)
c_ <- matrix(c(rep(1000, nSec - 1), 10), nrow = nLay)
Q <- matrix(rep(0, nLay * (nSec - 1)), nrow = nLay)
h <- c(-1.2, -0.2, 0.3, 0)
x <- c(830, 1185, 1300)
X <- seq(0, 1800, 10)
m <- solve_mlay1d(kD, c_, Q, h, x, X)
m %>% plot_mlay1d(X, layers = c(1:nLay))
m %>% plot_mlay1d(X, layers = c(1:nLay), ptype = "q")
m %>% plot_mlay1d(X, layers = c(1:nLay), ptype = "s")

