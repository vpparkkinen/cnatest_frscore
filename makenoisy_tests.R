if (is.na(Sys.getenv("RSTUDIO", unset = NA))) {
  setwd(system2("pwd", stdout = TRUE))
} else {
  path <- rstudioapi::getActiveDocumentContext()$path
  Encoding(path) <- "UTF-8"
  setwd(dirname(path))
}

library(cna)
source("helper_funcs.R")

g <- expand.grid(seq(0.3,0.8,by=0.1),seq(0.05,0.8,by=0.05), 8:20)

#ndats <- replicate(100, makenoisy_asf(prevalence = 0.2, noiselevel = .05, N=30))
ndats <- mapply(makenoisy_asf, 
                prevalence = g[,1], 
                noiselevel = g[,2], 
                N = g[,3], 
                SIMPLIFY = FALSE)

# are the noise rows really noise?

noiserows <- lapply(ndats, function(x) attributes(x)$noise)
input_dat <- lapply(ndats, function(x) attributes(x)$makenoisy_asf.input.data)
diffs_noise <- mapply(dplyr::setdiff, noiserows, input_dat, SIMPLIFY = FALSE)
allg <- mapply(dplyr::setequal, noiserows, diffs_noise)
all(allg)

# double check
models <- lapply(ndats, function(x) attributes(x)$makenoisy_asf.info$model)
cdat <- lapply(models, function(x) ct2df(selectCases(x))) 
diff_cdat_ndat <- mapply(dplyr::setdiff, x = ndats, y = cdat) 
all_gd <- mapply(dplyr::setequal, diff_cdat_ndat, noiserows)
all(all_gd)

prevs <- lapply(ndats, function(x) attributes(x)$makenoisy_asf.info$prevalence)
outs <- lapply(ndats, 
               \(x) gsub("=1","", attributes(x)$makenoisy_asf.info$outcome))

# are prevalences (nearly) what they should be?

r_prev <- mapply(\(x,y) mean(x[,y]), x = ndats, y = outs)
r_prev - g[,1]

# are noise levels what they should be?

r_noiseprop <- mapply(\(x,y) nrow(x) / nrow(y), x = noiserows, y = ndats)
r_noiseprop - g[,2]

rrnoise_prop <- mapply(\(x, y) nrow(x) / nrow(y), diff_cdat_ndat, ndats)
rrnoise_prop - g[,2]

