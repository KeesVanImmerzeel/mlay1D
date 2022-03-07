###_Install packages_###
packages.loading <-
  c(
    "dplyr",
    "installr",
    "rstudioapi",
    "dplyr",
    "magrittr"
  )
new.packages <-
  packages.loading[!(packages.loading %in% installed.packages()[, "Package"])]
if (length(new.packages))
  install.packages(new.packages, dependencies = TRUE)
lapply(packages.loading, require, character.only = TRUE)
installr::updateR() #If packages can't be installed update Rstudio

### Set work directory_###
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# seepage row
rownr <- 7

na_metzb <- readRDS('na/met zandbaan/Raai A na met zandbaan')
x <- 650
i <- which(na_metzb[1,]==x)
na_metzb_verhwrst <- readRDS('na/met zandbaan/met verhoogde weerstand tussen zandbaan en randsloot/Raai A na met zandbaan')

voor_metzb <- readRDS('voor/met zandbaan/Raai a voor met zandbaan')
voor_metzb_verhwrst <- readRDS('voor/met zandbaan/met verhoogde weerstand tussen zandbaan en randsloot/Raai A voor met zandbaan')

na_kwel_metzb <- na_metzb[rownr,i] * 1000 # mm/d
voor_kwel_metzb <- voor_metzb[rownr,i] * 1000 # mm/d
toename_kwel_metzb <- na_kwel_metzb -  voor_kwel_metzb

na_kwel_metzb_verhwrst <- na_metzb_verhwrst[rownr,i] * 1000 # mm/d
voor_kwel_metzb_verhwrst <- voor_metzb_verhwrst[rownr,i] * 1000 # mm/d
toename_kwel_metzb_verhwrst <- na_kwel_metzb_verhwrst - voor_kwel_metzb_verhwrst

m <- c(voor_kwel_metzb,na_kwel_metzb,voor_kwel_metzb_verhwrst,na_kwel_metzb_verhwrst) %>% matrix(nrow=2,ncol=2)
m %<>% rbind( m[2,] - m[1,] ) %>% round(digits=2)
