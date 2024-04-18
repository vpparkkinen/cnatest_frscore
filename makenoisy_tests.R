if (is.na(Sys.getenv("RSTUDIO", unset = NA))) {
  setwd(system2("pwd", stdout = TRUE))
} else {
  path <- rstudioapi::getActiveDocumentContext()$path
  Encoding(path) <- "UTF-8"
  setwd(dirname(path))
}

library(cna)
source("helper_funcs.R")

g <- expand.grid(seq(0.3,0.8,by=0.1),seq(0.05,0.8,by=0.1), 8:20)

#ndats <- replicate(100, makenoisy_asf(prevalence = 0.2, noiselevel = .05, N=30))
ndats <- mapply(makenoisy_asf, 
                prevalence = g[,1], 
                noiselevel = g[,2], 
                N = g[,3], 
                SIMPLIFY = FALSE)
noiserows <- lapply(ndats, function(x) attributes(x)$noise)
input_dat <- lapply(ndats, function(x) attributes(x)$makenoisy_asf.input.data)
diffs_noise <- mapply(dplyr::setdiff, noiserows, input_dat, SIMPLIFY = FALSE)
allg <- mapply(dplyr::setequal, noiserows, diffs_noise)
all(allg)

