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
         "udpipe",
         "stats",
         "tidyr",
         "openxlsx"
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
#' @return   Matrix with X (row 1), calculated heads ([L] rows (2-nLay+1)); lateral fluxes ([L2/T] rows (nLay+2 - 2xnLay+1)); seepage ([L/T] rows (2xnLay+2 - 3xnLay+1))
solve_mlay1d <- function(kD, c_, Q, h, x, X, f=1) {
      if (f > 1) {
            # Create sub-sections per intersection interval to increase numerical stability
            print(paste("Refine modelgrid f=",as.character(f)))
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
            exp1 <- expm::sqrtm(as.matrix(A[, , i]))
            exp2 <- expm::sqrtm(as.matrix(A[, , j]))
            C[ii:jj, ii:jj] <-
                  +expm::expm(-(x[j] - xMidSec[i]) * exp1)

            C[ii:jj, (ii + nLay):(jj + nLay)] <-
                  +expm::expm(+(x[j] - xMidSec[i]) * exp1)

            C[ii:jj, (ii + 2 * nLay):(jj + 2 * nLay)] <-
                  -expm::expm(-(x[j] - xMidSec[j]) * exp2)

            C[ii:jj, (ii + 3 * nLay):(jj + 3 * nLay)] <-
                  -expm::expm(+(x[j] - xMidSec[j]) * exp2)
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
            #if (crit>25) {
            #      m <- solve_mlay1d(kD, c_, Q[,2:(ncol(Q)-1)], h, x[2:(length(x)-1)], X, f=3) 
            #      return( m )
            #}
            exp3 <- expm::sqrtm(as.matrix(A[, , i]))
            exp4 <- expm::sqrtm(as.matrix(A[, , j]))
            C[(ii + nLay):(jj + nLay), (ii:jj)] <-
                  -di %*% exp3 %*% expm::expm(-exp1)
            C[(ii + nLay):(jj + nLay), (ii + nLay):(jj + nLay)] <-
                  +di %*% exp3 %*% expm::expm(exp1)
            C[(ii + nLay):(jj + nLay), (ii + 2 * nLay):(jj + 2 * nLay)] <-
                  +dj %*% exp4 %*% expm::expm(-exp2)
            C[(ii + nLay):(jj + nLay), (ii + 3 * nLay):(jj + 3 * nLay)] <-
                  -dj %*% exp4 %*% expm::expm(exp2)
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
      
      return(rbind(X, phi, q, s))
}

#' Plot results of function solve_mlay1d
#' @inheritParams solve_mlay1d
#' @param m Matrix with X and calculated heads, lateral fluxes and seepage
#' @param layers number(s) of layers to plot (numeric)
#' @param ptype Type of plot to create ("phi"=default, "q", "s")
#' @param labls optional data-frame with x-coordinates of (optional) vertical lines in plot (numeric) and labels (character)
#' return ggplot2 object
plot_mlay1d <-
      function(m,
               layers = 1,
               ptype = "phi",
               labls = NULL) {
            m %<>% t() %>% as.data.frame()
            X <- m[,1]
            m <- m[,2:ncol(m)]
            nlay <- ncol(m) / 3
            
            sel_names <- paste0(ptype, layers)
            names(m) <-
                  c(paste0("phi", 1:nlay),
                    paste0("q", 1:nlay),
                    paste0("s", 1:nlay))
            sel_names <- paste0(ptype, layers)
            m %<>% dplyr::select(all_of(sel_names))
            names(m) <- layers %>% as.character()
            m$X <- X
            m %<>% reshape2::melt(id.vars = "X", variable.name = "Aquifer") %>% tidyr::drop_na()
            
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
                              y = layer_scales(myplot)$y$range$range[1],
                              label = labls$txt,
                              angle = 90
                        )
            }
            return(myplot)
      }

# Convert data that is read with openxlsx::read.xlsx to mLay1D input
xls2mlay2D <- function(df) {
  res <- list()
  
  res$area <- df[1, ] %>% as.character()
  df <- df[2:nrow(df), ]  %>% mutate_if(is.character, as.numeric)
  res$nLay <- (nrow(df) - 2) / 3
  res$nSec <- ncol(df)
  res$x <- df[1, 1:(res$nSec - 1)] %>% as.double() %>% as.vector()
  res$h <- df[2, ] %>% as.double() %>% as.vector()
  i <- seq(from = 3,
           to = res$nLay * 3,
           length.out = res$nLay)
  res$c_ <- df[i, ] %>% data.matrix()
  res$kD <- df[i + 1, ] %>% data.matrix()
  res$Q <- df[i + 2, ] %>% data.matrix()
  i <- c(1, which((res$area == udpipe::txt_previous(res$area, n = 1)) != TRUE))
  res$labls <- data.frame(xvlines = res$x[i], txt = res$area[i])
  res$labls$xvlines[1] <- NA
  return(res)
}

######################################################################

# Hulp functies om gemiddelde kwel (mm) te berekenen afhankelijk van de weerstand van de deklaag
.avseepage <-  function( cnew, res, xrange ){
  m <- solve_mlay1d(res$kD, c_=rbind(res$c[1,],cnew), res$Q, res$h, res$x, X = res$x)
  i <- which(res$x > xrange[1] & res$x < xrange[2])
  avseepage <- mean(m[7,i]) * 1000 %>% round(digits=2) 
  return(avseepage)
}
avseepage <- Vectorize(.avseepage, vectorize.args="cnew", USE.NAMES=FALSE)

######################################################################

# Bereken dikte van de neerslaglens 
#' @param x Afstand tot het midden van het perceel
#' @param L Slootafstand
#' @param kx horizontale doorlatendheid (m/d)
#' @param kz verticale doorlatendheid (m/d)
#' @param N Neerslagoverschot (m/d)
#' @param K Jaargemiddelde kwelintensiteit (m/d)
.lensdikte <- function(x, L, kx, kz, N, K){
  if ( x <= 0 | x >= L/2 ) {
    return(NA)
  }
  alpha <- N/(K+N)
  ksi <- x * pi / L
  dikte <- .5*(L/pi)*sqrt(kz/kx)*log(sin((1+alpha)*ksi)/sin((1-alpha)*ksi))
  return(dikte)
}
lensdikte <- Vectorize(.lensdikte, vectorize.args="x", USE.NAMES=FALSE)

### gemiddelde lensdikte
gem_lensdikte <- function( x, L, kx, kz, N, K) {
  mean(lensdikte(x, L, kx, kz, N, K))
}

.hlp_gem_lensdikte <- function( x, L, c_deklaag, N, K, Ddeklaag=12, ani_factor=10 ) {
  kz <- Ddeklaag / c_deklaag
  kx <- ani_factor * kz
  gem_lensdikte( x, L, kx, kz, N, K)
}

hlp_gem_lensdikte <- Vectorize(.hlp_gem_lensdikte, vectorize.args=c("c_deklaag","K"), USE.NAMES=FALSE)

### max lensdikte
max_lensdikte <- function( x, L, kx, kz, N, K) {
  max(lensdikte(x, L, kx, kz, N, K))
}

.hlp_max_lensdikte <- function( x, L, c_deklaag, N, K, Ddeklaag=12, ani_factor=10 ) {
  kz <- Ddeklaag / c_deklaag
  kx <- ani_factor * kz
  max_lensdikte( x, L, kx, kz, N, K)
}
hlp_max_lensdikte <- Vectorize(.hlp_max_lensdikte, vectorize.args=c("c_deklaag","K"), USE.NAMES=FALSE)

### min lensdikte
min_lensdikte <- function( x, L, kx, kz, N, K) {
  min(lensdikte(x, L, kx, kz, N, K))
}

.hlp_min_lensdikte <- function( x, L, c_deklaag, N, K, Ddeklaag=12, ani_factor=10 ) {
  kz <- Ddeklaag / c_deklaag
  kx <- ani_factor * kz
  min_lensdikte( x, L, kx, kz, N, K)
}
hlp_min_lensdikte <- Vectorize(.hlp_min_lensdikte, vectorize.args=c("c_deklaag","K"), USE.NAMES=FALSE)

###

# Schat fractie van het perceelsoppervlak waar de neerslaglens in de zomer verdwijnt
.P_lensdikte <- function(ymax, ycrit) {
  if (ymax > ycrit) {
      return(1-sqrt((1-ycrit/ymax)))
    } else {
      return(1)
    }
}
P_lensdikte <- Vectorize(.P_lensdikte, vectorize.args="ymax")


######################################################################

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
#kD <- matrix(rep(c(1, 1, 1), nSec), nrow = nLay)
c_ <- matrix(rep(c(50, 400, 85 / 0.075), nSec), nrow = nLay)
Q <- matrix(rep(0, nLay * (nSec - 1)), nrow = nLay)
x <- c(-1000, 1000, 3250, 4500, 5500, 6500, 7250, 8750, 9750, 10500)
X <- seq(-2500, 11000, 10)
m <- solve_mlay1d(kD, c_, Q, h, x, X, f=10)
m %>% plot_mlay1d(layers = c(1:nLay))
m %>% plot_mlay1d(layers = c(1:nLay), ptype = "q")
m %>% plot_mlay1d(layers = c(1:nLay), ptype = "s")

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
m %>% plot_mlay1d(layers = c(1:nLay))
m %>% plot_mlay1d(layers = c(1:nLay), ptype = "q")
m %>% plot_mlay1d(layers = c(1:nLay), ptype = "s")

### Example 3. Use spreadsheet as input.
### Copy range B1:BW9 from spreadsheet 'Example 3.xlsx' to the clipboard.
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
m %>% plot_mlay1d(
                  layers = c(1:nLay),
                  labls = labls)
m %>% plot_mlay1d(
                  layers = c(1:nLay),
                  ptype = "q",
                  labls = labls)
m %>% plot_mlay1d(
                  layers = c(1:nLay),
                  ptype = "s",
                  labls = labls)
fname <- "Dwarsraai a/voor/met zandbaan/met verhoogde weerstand tussen zandbaan en randsloot/Raai a voor met zandbaan"
fname <- "Dwarsraai a/voor/met zandbaan/Raai a voor met zandbaan"
fname <- "Dwarsraai a/na/met zandbaan/met verhoogde weerstand tussen zandbaan en randsloot/Raai a na met zandbaan"
fname <- "Dwarsraai a/na/met zandbaan/Raai a na met zandbaan"

fname <- "Dwarsraai b/voor/met zandbaan/met verhoogde weerstand tussen zandbaan en randsloot/Raai b voor met zandbaan"
fname <- "Dwarsraai b/voor/met zandbaan/Raai b voor met zandbaan"
fname <- "Dwarsraai b/na/met zandbaan/met verhoogde weerstand tussen zandbaan en randsloot/Raai b na met zandbaan"
fname <- "Dwarsraai b/na/met zandbaan/Raai b na met zandbaan"

m %>% saveRDS(fname)

# Example 4: Create sub-sections (50) per intersection interval to increase numerical stability.
nLay <- 4
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
kD <- matrix(rep(c(35 * 30, 80 * 30, (55 / 2) * 0.075, 1), nSec), nrow = nLay)
#kD <- matrix(rep(c(1, 1, 1), nSec), nrow = nLay)
c_ <- matrix(rep(c(50, 400, 85 / 0.075, 1), nSec), nrow = nLay)
Q <- matrix(rep(0, nLay * (nSec - 1)), nrow = nLay)
x <- c(-1000, 1000, 3250, 4500, 5500, 6500, 7250, 8750, 9750, 10500)
X <- seq(-2500, 11000, 10)
m <- solve_mlay1d(kD, c_, Q, h, x, X, f=50) # Results are produced.
m <- solve_mlay1d(kD, c_, Q, h, x, X) # No results produced.
m %>% plot_mlay1d(layers = c(1:nLay))
m %>% plot_mlay1d(layers = c(1:nLay), ptype = "q")
m %>% plot_mlay1d(layers = c(1:nLay), ptype = "s")

#####################################################################################

# Example 5: read directly from spreadsheet 
#selct <- c("A","voor")
#selct <- c("A","na")
#selct <- c("B","voor")
selct <- c("B","na")

# Dwarsraai A
if (selct[1] == "A") {
  xrange <-
    c(334, 625) # Afstanden in de raai waarop berekeningsresultaten worden gemiddeld etc.
  if (selct[2] == "voor") {
    fname1 <- "Dwarsraai A/voor/Voor - Dwarsraai A.xlsx"
    fname2 <- "Dwarsraai A/voor/met zandbaan/dikte_neerslaglens_voor.csv"
  } else {
    fname1 <- "Dwarsraai A/na/Na - Dwarsraai A.xlsx"
    fname2 <- "Dwarsraai A/na/met zandbaan/dikte_neerslaglens_na.csv"
  }
} else {
  # Dwarsraai B
  xrange <-
    c(370, 750) # Afstanden in de raai waarop berekeningsresultaten worden gemiddeld etc.
  if (selct[2] == "voor") {
    fname1 <- "Dwarsraai B/voor/Voor - Dwarsraai B.xlsx"
    fname2 <- "Dwarsraai B/voor/met zandbaan/dikte_neerslaglens_voor.csv"
  } else {
    fname1 <- "Dwarsraai B/na/Na - Dwarsraai B.xlsx"
    fname2 <- "Dwarsraai B/na/met zandbaan/dikte_neerslaglens_na.csv"
  }
}

df <- openxlsx::read.xlsx(fname1,rows=c(1:9),colNames=FALSE) %>% dplyr::select(-1)
df <- df[,2:ncol(df)]
inp <- df %>% xls2mlay2D()
#m <- solve_mlay1d(inp$kD, inp$c_, inp$Q, inp$h, inp$x, X = inp$x)
#m %>% plot_mlay1d(
#  layers = c(1:inp$nLay),
#  labls = inp$labls)
#m %>% plot_mlay1d(
#  layers = c(1:inp$nLay),
#  ptype = "q",
#  labls = inp$labls)
#m %>% plot_mlay1d(
#  layers = c(1:inp$nLay),
#  ptype = "s",
#  labls = inp$labls)

# Gevoeligheidsanalyse voor deklaagweerstand van kwel en dikte neerslaglens in randzone.
# i <- which(trimws(res$area)=="Hoofdwatergang")

cnew <- seq(1000,4000,100)
seepage <- avseepage( cnew, inp, xrange )
df <- data.frame("c_deklaag"=cnew, "kwelmm"=seepage)

# Gevoeligheidsanalyse voor deklaagweerstand en dikte neerslaglens in randzone.
L <- 7.5 # Slootafstand
x <- seq(.1,L/2-.1,length.out=25) # Afstand tot het midden van het perceel
N <- 0.3 # Gemiddelde jaarlijkse grondwateraanvulling gras (tabel 1.7 grondwaterzakboekje)
ycrit <- 1.25 # 0.1 m neerslagtekort gedeeld door freatische bergingscoefficient lichte klei 0.08 (-)
df$min_lensdikte <- hlp_min_lensdikte( x, L, c_deklaag=df$c_deklaag, N, K=df$kwelmm/1000 )
df$gem_lensdikte <- hlp_gem_lensdikte( x, L, c_deklaag=df$c_deklaag, N, K=df$kwelmm/1000) 
df$max_lensdikte <- hlp_max_lensdikte( x, L, c_deklaag=df$c_deklaag, N, K=df$kwelmm/1000 )
df$P <- P_lensdikte(ymax=df$max_lensdikte, ycrit)
df %>% write.table(fname2, sep="\t", row.names=FALSE)



