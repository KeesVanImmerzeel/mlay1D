###_Install packages_###
packages.loading <-
      c( "expm",
         "installr",
         "rstudioapi",
         "pracma",
         "reshape2",
         "dplyr",
         "magrittr",
         "ggplot2",
         "clipr",
         "udpipe"
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
      #      COEF <- limSolve::Solve.banded(C[, (nLay + 1):(ncol(C) - nLay)], R) %>%
      #                   c(rep(0, nLay), ., rep(0, nLay))
      
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
#' @param labls optional data-frame with x-coordinates of (optional) vertical lines in plot (numeric) and labels (character)
#' return ggplot2 object
plot_mlay1d <-
      function(m,
               X,
               layers = 1,
               ptype = "phi",
               labls = NULL) {
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
                  theme(plot.title = element_text(hjust = 0.5))
            if (!is.null(labls)) {
                  xintercept <- labls$xvlines[!is.na(labls$xvlines)]
                  myplot <-
                        myplot + geom_vline(
                              xintercept = xintercept,
                              linetype = "dotted",
                              color = "blue"
                        )
                  xrnge <- layer_scales(myplot)$x$range$range
                  xv <- c(xintercept, xrnge) %>% sort()
                  xv <- xv[1:(length(xv) - 1)] + diff(xv) / 2
                  myplot <-
                        myplot + annotate(
                              "text",
                              x = xv,
                              y = mean(layer_scales(myplot)$y$range$range),
                              label = labls$txt,
                              angle = 90
                        )
            }
            return(myplot)
      }

### Copy range B1:BW9 from spreadsheet 'Dwarsraai A.xlsx' (or 'Dwarsraai B.xlsx') to the clipboard.
### Then run the following code.
df <- clipr::read_clip_tbl(header = FALSE, dec = ".")
area <- df[1, ] %>% as.character()
df <- df[2:nrow(df), ] %>% mutate_if(is.character, as.numeric)
nLay <- (nrow(df) - 2) / 3
nSec <- ncol(df)
x <- df[1, 1:(nSec - 1)] %>% as.double() %>% as.vector()
h <- df[2, ] %>% as.double() %>% as.vector()
i <- seq(from = 3,
         to = nLay * 3,
         length.out = nLay)
c_ <- df[i, ] %>% data.matrix()
kD <- df[i + 1, ] %>% data.matrix()
Q <- df[i + 2, ] %>% data.matrix()
m <- solve_mlay1d(kD, c_, Q, h, x, X = x)
i <- c(1, which((area == udpipe::txt_previous(area, n = 1)) != TRUE))
labls <- data.frame(xvlines = x[i], txt = area[i])
labls$xvlines[1] <- NA
m %>% plot_mlay1d(X = x,
                  layers = c(1:nLay),
                  labls = labls)
m %>% plot_mlay1d(
  X = x,
  layers = c(1:nLay),
  ptype = "q",
  labls = labls
)
m %>% plot_mlay1d(
  X = x,
  layers = c(1:nLay),
  ptype = "s",
  labls = labls
)
