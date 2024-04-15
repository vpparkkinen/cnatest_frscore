if(is.na(Sys.getenv("RSTUDIO", unset = NA))){
  setwd(system2("pwd", stdout = TRUE)) # if not in RStudio, assume R runs in
} else {                               # a shell. otherwise assume RStudio
  path <- rstudioapi::getActiveDocumentContext()$path
  Encoding(path) <- "UTF-8"
  setwd(dirname(path))
}

source("cnaTest.R")
library(doParallel)
library(frscore)
n_cores <- detectCores() - 2

flipout <- function(data, outcome, proportion){
  N <- nrow(data)
  out_col <- which(names(data) == outcome)
  range_o <- min(data[,out_col]):max(data[,out_col])
  n_to_flip <- round(N*proportion)
  n_to_flip <- if(n_to_flip == 0L) 1L else n_to_flip
  w_rows <- sample(1:N, n_to_flip)
  for(row in w_rows){
    ov <- data[row, out_col]
    not_ov <- range_o[range_o != ov]
    change_ov_to <- sample(not_ov, 1)
    data[row, out_col] <- change_ov_to
  }
  return(data)
}

n_models <- 10
n_factors <- 6
noise_prop <- 0.3
targets <- replicate(n_models, randomAsf(n_factors), simplify = FALSE)
outcomes <- lapply(targets, rhs)

clean_dsets <- lapply(targets, function(x) ct2df(selectCases(x)))
noisy_dsets <- mapply(flipout, 
                      data = clean_dsets, 
                      outcome = outcomes,
                      proportion = noise_prop, 
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
                            output = "asf",
                            comp.method = "is.submodel",
                            SIMPLIFY = FALSE)

top_frscore <- lapply(frscore_results, function(x) x[[1]][1,2])
correctness <- mcmapply(is.submodel, 
                        x = top_frscore, 
                        y = targets,
                        SIMPLIFY = FALSE)

correctness <- unlist(correctness)
cor_all_avg <- sum(correctness) / n_models

sig_correctness <- sum(correctness[sig]) / length(sig)

cor_all_avg
sig_correctness


