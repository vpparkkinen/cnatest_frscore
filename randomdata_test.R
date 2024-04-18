if (is.na(Sys.getenv("RSTUDIO", unset = NA))) {
  setwd(system2("pwd", stdout = TRUE))
} else {
  path <- rstudioapi::getActiveDocumentContext()$path
  Encoding(path) <- "UTF-8"
  setwd(dirname(path))
}

if(!exists("cnaTest")) source("cnaTest.R")
if(!exists("flipout")) source("helper_funcs.R")
library(doParallel)
library(frscore)
n_cores <- detectCores() - 2
options(mc.cores = n_cores)

r_datasets <- bs_dat_create(Nsets = 500, size = 8)

outcomes <- lapply(r_datasets, \(x) sample(names(x), 1))

pvals_temp <- mcmapply(cnaTest,
                       d = r_datasets,
                       outcome = outcomes,
                       SIMPLIFY = FALSE)

pvals <- unlist(lapply(pvals_temp, `[[`, 8))

sig <- which(pvals <= 0.05)
length(sig) / length(r_datasets)
sig_dsets <- r_datasets[sig]
sig_outcomes <- outcomes[sig]

nonsig_dsets <- r_datasets[-sig]
nonsig_outcomes <- outcomes[-sig]






