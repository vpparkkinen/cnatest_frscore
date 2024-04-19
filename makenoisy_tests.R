if (is.na(Sys.getenv("RSTUDIO", unset = NA))) {
  setwd(system2("pwd", stdout = TRUE))
} else {
  path <- rstudioapi::getActiveDocumentContext()$path
  Encoding(path) <- "UTF-8"
  setwd(dirname(path))
}

library(cna)
library(doParallel)
source("helper_funcs.R")
n_cores <- detectCores() - 2
options(mc.cores = n_cores)

# set up test
g <- expand.grid(seq(0.05, 0.9, by = 0.05),
                 seq(0.05, 0.9, by = 0.05),
                 8:128)


ndats <- mcmapply(makenoisy_asf,
                  prevalence = g[, 1],
                  noiselevel = g[, 2],
                  N = g[, 3],
                  SIMPLIFY = FALSE)

# are the noise rows really noise?
noiserows <- lapply(ndats, function(x) attributes(x)$noise)
input_dat <- lapply(ndats,
  function(x) {
    attributes(x)$makenoisy_asf.input.data
  }
)
diffs_noise <- mcmapply(dplyr::setdiff,
                        noiserows,
                        input_dat,
                        SIMPLIFY = FALSE)
allg <- mcmapply(
  \(x, y) dplyr::setequal(x, unique(y))
  , x = diffs_noise, y = noiserows
)
all(allg)

# double check
models <- mclapply(ndats,
                   \(x) attributes(x)$makenoisy_asf.info$model)
cdat <- mclapply(models, \(x) ct2df(selectCases(x)))
diff_cdat_ndat <- mcmapply(dplyr::setdiff, x = ndats, y = cdat)
all_gd <- mcmapply(
  \(x, y) dplyr::setequal(x, unique(y)),
  x = diff_cdat_ndat,
  y = noiserows
)
all(all_gd)


# are prevalences (nearly) what they should be?
prevs <- mclapply(ndats,
  \(x) {
    attributes(x)$makenoisy_asf.info$prevalence
  }
)
outs <- mclapply(ndats,
  \(x) gsub("=1", "", attributes(x)$makenoisy_asf.info$outcome)
)
r_prev <- mcmapply(\(x, y) mean(x[, y]), x = ndats, y = outs)
prevdiff <- r_prev - g[, 1]
mean(prevdiff)
max(prevdiff)
which(prevdiff == max(prevdiff))


# are noise levels (nearly) what they should be?
r_noiseprop <- mcmapply(\(x, y) nrow(x) / nrow(y), x = noiserows, y = ndats)
noisediff <- r_noiseprop - g[, 2]
mean(noisediff)
max(abs(noisediff))
which(noisediff == max(noisediff))

# are sample sizes what they should be?
ns <- mclapply(ndats,
               \(x) attributes(x)$makenoisy_asf.info$N)
r_ns <- mclapply(ndats, nrow)
identical(ns, r_ns)
