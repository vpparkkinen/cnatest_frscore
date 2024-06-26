if (is.na(Sys.getenv("RSTUDIO", unset = NA))) {
  setwd(system2("pwd", stdout = TRUE))
} else {
  path <- rstudioapi::getActiveDocumentContext()$path
  Encoding(path) <- "UTF-8"
  setwd(dirname(path))
}

source("cnaTest.R")
source("helper_funcs.R")
library(doParallel)
library(frscore)
n_cores <- detectCores() - 2
options(mc.cores = n_cores)

set.seed(25)
n_models <- 1000
n_factors <- 7
noise_prop <- 0.2
tquantile <- 1
N <- 40
rand_outcome <- TRUE



targets <- replicate(n_models, randomAsf(n_factors), simplify = FALSE)

outcomes <- lapply(targets, rhs)


clean_dsets <- mclapply(targets,
                      function(x) cna::some(ct2df(selectCases(x)), N))

for (i in seq_along(clean_dsets)) {
  clean_dsets[[i]] <- cbind(clean_dsets[[i]],
        data.frame(X = rbinom(nrow(clean_dsets[[i]]), 1, 0.5)))
}

if (rand_outcome){
  cat("!!! outcome selected at random from DGS non-outcomes !!!\n!!! is this what you want? !!!")
  outcomes <- mapply(\(x, y) {sample(names(x)[names(x) != y], 1)},
                     x = clean_dsets,
                     y = outcomes,
                     SIMPLIFY = FALSE)
}



noisy_dsets <- mcmapply(flipout,
                      data = clean_dsets,
                      outcome = outcomes,
                      MoreArgs = list(proportion = noise_prop),
                      SIMPLIFY = FALSE)

pvals_temp <- mcmapply(cnaTest,
                     d = noisy_dsets,
                     outcome = outcomes,
                     SIMPLIFY = FALSE)

pvals <- unlist(lapply(pvals_temp, `[[`, 8))

sig <- which(pvals < 0.05)

frscore_results <- mcmapply(frscored_cna,
                            x = noisy_dsets,
                            outcome = outcomes,
                            MoreArgs = list(output = "asf",
                            comp.method = "is.submodel",
                            fit.range = c(0.9, 0.6)),
                            SIMPLIFY = FALSE)

top_frscore <- lapply(
  frscore_results,
  function(x) x[[1]][x[[1]]$score >= quantile(x[[1]]$score, tquantile) ,2]
  )

correctness <- mcmapply(any_submodel,
                        x = top_frscore,
                        y = targets,
                        SIMPLIFY = FALSE)

correctness <- unlist(correctness)
all_cor <- sum(correctness, na.rm = TRUE)
cor_all_avg <- all_cor / n_models

sig_all_cor <- sum(correctness[sig], na.rm = TRUE)
sig_correctness <- sig_all_cor / length(sig)

cor_all_avg
sig_correctness
length(sig)
length(correctness)
nresults <- lapply(top_frscore, function(x) length(x))
mean(unlist(nresults))
max(unlist(nresults))

# average complexity of models when correct model among top results, both in all results, and in those for sig data sets
correctness[is.na(correctness)] <- FALSE
topf_vec <- unlist(top_frscore[correctness])
sig_cor_idx <- correctness[sig]


complx_allcor <- unlist(lapply(topf_vec, cna::getComplexity))
sig_cor_complx <- unlist(lapply(topf_vec[sig_cor_idx], cna::getComplexity))

mean(complx_allcor)
mean(sig_cor_complx)

# same for false results

false_topf_vec <- unlist(top_frscore[!correctness])
complx_allfalse <- unlist(lapply(false_topf_vec, cna::getComplexity))
mean(complx_allfalse)

sig_incor_idx <- !correctness[sig]
incorrects <- unlist(top_frscore[!correctness[sig]])
sig_incor_complx <- unlist(lapply(topf_vec[sig_incor_idx], cna::getComplexity))
mean(sig_incor_complx)

# complexity of models for significant vs insignificant data sets, regardless of correctness

sig_results <- unlist(top_frscore[sig])
mean(unlist(lapply(sig_results, cna::getComplexity)))
mean(unlist(lapply(top_frscore[-sig], cna::getComplexity)))
