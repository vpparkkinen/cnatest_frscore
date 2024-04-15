
# Functions related to permutation test
# M. Ambuehl, March 2024

suppressPackageStartupMessages(suppressWarnings({
  require(dplyr)
  require(cna)
  require(collapse)
  library(ggplot2)
  loadNamespace("matrixStats")
}))

# ------------------------------------------------------------------------------

#' @param fast Logical: FALSE = implementation using dplyr (slower); TRUE = using collapse (fast).
minErr1 <- function(d, outcome, nfact, fast = TRUE){
  if (inherits(d, "configTable")) 
    stop("Doesn't currently accept 'configTable's - use ct2df() to create a pure data.frame.")
  stopifnot(
    is.data.frame(d), d>=0, d<=1, 
    length(outcome) == 1, outcome %in% names(d), 
    nfact>0, nfact <= length(d)-1)
  factors <- setdiff(names(d), outcome)
  if (!all(d %in% 0:1)){
    d[factors] <- round(d[factors])
  }
  n <- nrow(d)
  ifactors <- match(factors, names(d))
  factCombs <- t(combn(ifactors, nfact))
  ii <- cbind(rep(seq_len(n), length(factCombs)), 
              rep(c(factCombs), each = n))
  x <- as.vector(matrix(d[ii], ncol = nfact) %*% 2^(seq(nfact)-1))
  df <- data.frame(x, 
                   y = rep(d[[outcome]], nrow(factCombs)),
                   g = rep(1:nrow(factCombs), each = n))
  if (fast){
    grps <- GRP(df, by = c("g", "x"))
    interim1 <- grps$groups
    yy <- fsum(array(c(df$y, 1-df$y), c(nrow(df), 2)), 
               grps, use.g.names = FALSE)
    interim1$err <- matrixStats::rowMins(yy) 
    calc <- fsum(interim1$err, interim1$g, use.g.names = FALSE)/n
  } else {
    calc <- df %>% 
      group_by(x, g) %>% 
      summarize(err = min(sum(y), sum(1-y)), 
                .groups = "drop") %>% 
      group_by(g) %>% 
      summarize(err = sum(err)/n) %>% 
      pull(err)
  }
  out <- data.frame(t(combn(factors, nfact)))
  names(out) <- paste0("fact", seq_len(nfact))
  out$err <- calc
  out[order(out$err), , drop = FALSE]
}

minErr2 <- function(d, outcome, nfact = 1:min(length(d)-1, 6), ...){
  minErrList <- lapply(nfact, minErr1, d = d, outcome = outcome, ...)
  minErrList <- lapply(minErrList, function(x) x[which.min(x$err), ])
  out <- data.frame(
    nfact, 
    factors = NA,
    err = sapply(minErrList, "[[", "err"))
  out$factors <- split.default(unlist(mapply(function(d, j) d[, 1:j, drop = TRUE], minErrList, nfact), use.names = FALSE), 
                               rep(seq_along(nfact), nfact))
  out
}

# ==============================================================================

minErr1sim <- function(d, outcome, nfact, nsim = 100, 
                       ysim = replicate(nsim, sample(d[[outcome]]), simplify = TRUE), 
                       fast = TRUE){
  factors <- setdiff(names(d), outcome)
  if (!all(d %in% 0:1)){
    d[factors] <- round(d[factors])
  }
  n <- nrow(d)
  ifactors <- match(factors, names(d))
  factCombs <- t(combn(ifactors, nfact))
  ii <- array(c(rep(seq_len(n), length(factCombs)), 
                rep(c(factCombs), each = n)), 
              dim = c(n*length(factCombs), 2))
                
  x <- as.vector(matrix(d[ii], ncol = nfact) %*% 2^(seq(nfact)-1))
  df0 <- data.frame(
    x, 
    y = rep(d[[outcome]], nrow(factCombs)),
    g = rep(1:nrow(factCombs), each = n))
  # ----- identical with simErr1() up to here -----
  df <- rep_rows(df0, rep(seq_len(nrow(df0)), nsim))
  df$y <- as.vector(ysim[rep(seq_len(nrow(ysim)), nrow(factCombs)), ])
  df$sim <- rep(seq_len(nsim), each = nrow(df0))
  if (fast){
    grps <- GRP(df, by = c("sim", "g", "x"))
    interim1 <- grps$groups
    yy <- fsum(array(c(df$y, 1-df$y), c(nrow(df), 2)), 
               grps, use.g.names = FALSE)
    interim1$err <- matrixStats::rowMins(yy) 
    grps <- GRP(interim1, c("sim", "g"))
    interim2 <- grps$groups
    interim2$err <- fsum(interim1$err, grps, use.g.names = FALSE)/n
    fmin(interim2$err, interim2$sim, use.g.names = FALSE)
  } else {
  df %>%
    group_by(sim, g, x) %>%
    summarize(err = min(sum(y), sum(1-y)),
              .groups = "drop") %>%
    group_by(sim, g) %>%
    summarize(err = sum(err)/n, .groups = "drop") %>%
    group_by(sim) %>%
    summarize(err = min(err)) %>%
    pull(err)
  }
}

rep_rows <- function(x, i) {
  out <- lapply(x, function(z) if (length(dim(z)) != 2L) z[i] else z[i, , drop = FALSE])
  attr(out, "row.names") <- .set_row_names(length(i))
  class(out) <- "data.frame"
  out
}

minErr2sim <- function(d, outcome, nfact = 1:min(length(d)-1, 6), nsim = 100, 
                       ysim = replicate(nsim, sample(d[[outcome]]), simplify = TRUE), 
                       ...){
  stopifnot(
    is.data.frame(d), d>=0, d<=1, 
    length(outcome) == 1, outcome %in% names(d), 
    nfact>0, nfact <= length(d)-1, 
    identical(dim(ysim), c(nrow(d), as.integer(nsim))))
  lapply(nfact, 
         function(i) minErr1sim(d, outcome, nfact = i, nsim = nsim, ysim = ysim, ...))
}

#' cnaTest
#' 
cnaTest <- function(d, outcome, nfact = 1:min(ncol(d)-1, 6), nsim = 2000, 
                    ysim = replicate(nsim, sample(d[[outcome]]), simplify = TRUE), 
                    ...){
  tm0 <- proc.time()[[3]]
  err_data_detailed <- minErr2(d, outcome, nfact = nfact, ...)$err
  err_data <- mean(err_data_detailed)
  err_ref_detailed <- do.call(rbind, minErr2sim(d, outcome, nfact = nfact, 
                                                nsim = nsim, ysim = ysim, ...))
  err_ref <- colMeans(err_ref_detailed)
  eps <- .Machine$double.eps * 10
  out <- list(call = match.call(), nfact = nfact, nsim = nsim, 
              err_data = err_data, err_data_detailed = err_data_detailed,
              err_ref = err_ref, err_ref_detailed = err_ref_detailed, 
              p_value = mean(err_data >= err_ref - eps))
  out$execution_time <- proc.time()[[3]] - tm0
  class(out) <- "cnaTest"
  out
}

print.cnaTest <- function(x, ...){
  x0 <- x
  x$err_ref <- x$err_ref_detailed <- x$err_data_detailed <- NULL
  print(unclass(x, ...))
  invisible(x0)
}
`[.cnaTest` <- function(x, i){
  match_i <- match(i, x$nfact, 0)
  x$nfact <- x$nfact[match_i]
  x$call$nfact <- x$nfact
  x$err_data_detailed <- x$err_data_detailed[match_i]
  x$err_data <- mean(x$err_data_detailed)
  x$err_ref_detailed <- x$err_ref_detailed[match_i, , drop = FALSE]
  x$err_ref <- colMeans(x$err_ref_detailed)
  x$p_value <- mean(x$err_data >= x$err_ref)
  x
}
summary.cnaTest <- function(object, ...){
  nf <- object$nfact
  smry <- as.data.frame(t(sapply(nf, function(i){
    unclass(object[i])[c("nfact", "nsim", "err_data", "p_value")]
  })))
  list(call = object$call, nsim = object$nsim, summary = smry, 
       execution_time = object$execution_time)
}
    


plot_violin <- function(x, nfact = c(0, x$nfact), width = .5, ...){
  stopifnot(inherits(x, "cnaTest"))
  nfact_keep <- ifelse(nfact == 0, "Overall", as.character(nfact))
  d0 <- data.frame(
    nfact = factor(c("Overall", x$nfact), levels = c("Overall", x$nfact)), 
    err = c(x$err_data, x$err_data_detailed)) %>% 
    filter(nfact %in% nfact_keep)
  df <- data.frame(
    err = as.vector(rbind(x$err_ref, x$err_ref_detailed)), 
    nfact = rep(c("Overall", as.character(x$nfact)), x$nsim)) %>% 
    mutate(nfact = factor(nfact, levels = c("Overall", x$nfact))) %>% 
    filter(nfact %in% nfact_keep) %>% 
    droplevels
  if ("Overall" %in% nfact_keep){
    levels(df$nfact)[[1]] <- levels(d0$nfact)[[1]] <- 
      paste0("{", paste0(x$nfact, collapse = ","), "}")
  }
  i_nfact <- seq_along(levels(df$nfact))
  out <- df %>% 
    ggplot(aes(nfact, err, group = nfact)) +
    geom_violin(color = 4, fill = 4, alpha = .25, width = .6) +
    geom_point(data = d0, col = 2, size = 3) +
    geom_segment(data = d0, 
                 aes(x = i_nfact - width/2, xend = i_nfact + width/2,
                     y = err, yend = err), 
                 col = 2, linewidth = 1.5)
    
  if (length(nfact_keep)>1 && "Overall" %in% nfact_keep){
    out <- out + geom_vline(xintercept = 1.5)
  }
  out
}

plot_bars <- function(x, nfact = c(0, x$nfact), width = 3, ...){
  stopifnot(inherits(x, "cnaTest"))
  nfact_keep <- ifelse(nfact == 0, "Overall", as.character(nfact))
  d0 <- data.frame(
    nfact = factor(c("Overall", x$nfact), levels = c("Overall", x$nfact)), 
    err = c(x$err_data, x$err_data_detailed)) %>% 
    filter(nfact %in% nfact_keep)
  df <- data.frame(
    err = as.vector(rbind(x$err_ref, x$err_ref_detailed)), 
    nfact = rep(c("Overall", as.character(x$nfact)), x$nsim)) %>% 
    mutate(nfact = factor(nfact, levels = c("Overall", x$nfact))) %>% 
    filter(nfact %in% nfact_keep) %>% 
    droplevels
  if ("Overall" %in% nfact_keep){
    levels(df$nfact)[[1]] <- levels(d0$nfact)[[1]] <- 
      paste0("{", paste0(x$nfact, collapse = ","), "}")
  }
  tbl <- df %>% 
    mutate(err = round(err, 12)) %>% 
    count(nfact, err) %>% 
    group_by(nfact) %>% 
    mutate(prop = n / sum(n))
  max_distinct <- tbl %>% count(nfact) %>% pull(n) %>% max
  if (max_distinct>30) warning("bars might not be a suitable type of plot.")
  tbl %>% 
    ggplot(aes(err, prop)) + 
    facet_wrap(~ paste0("nfact = ", nfact), scales = "free_y") +
    geom_segment(aes(xend = err, yend = 0), 
                 col = 4, linewidth = width, alpha = .7) +
    geom_vline(data = d0, mapping = aes(xintercept = err), 
               col = 2, alpha = 0.6, linewidth = 1) +
    theme(panel.spacing = unit(20, "pt"))
}

plot_ecdf <- function(x, nfact = c(0, x$nfact), width = 1, ...){  
  stopifnot(inherits(x, "cnaTest"))
  nfact_keep <- ifelse(nfact == 0, "Overall", as.character(nfact))
  d0 <- data.frame(
    nfact = factor(c("Overall", x$nfact), levels = c("Overall", x$nfact)), 
    err = c(x$err_data, x$err_data_detailed)) %>% 
    filter(nfact %in% nfact_keep)
  df <- data.frame(
    err = 
      as.vector(rbind(x$err_ref, x$err_ref_detailed)), 
    nfact = rep(c("Overall", as.character(x$nfact)), x$nsim)) %>% 
    mutate(nfact = factor(nfact, levels = c("Overall", x$nfact))) %>% 
    filter(nfact %in% nfact_keep) %>% 
    droplevels
  subTests <- c(list(x), lapply(x$nfact, function(i) x[i]))
  subTests <- subTests[c("Overall", x$nfact) %in% nfact_keep]
  d0$err_data <- sapply(subTests, "[[", "err_data")
  d0$p_value <- sapply(subTests, "[[", "p_value")
  if ("Overall" %in% nfact_keep){
    levels(df$nfact)[[1]] <- levels(d0$nfact)[[1]] <- 
      paste0("{", paste0(x$nfact, collapse = ","), "}")
  }
  df %>% 
    ggplot(aes(err)) + 
    facet_wrap(~ paste0("nfact = ", nfact), scales = "free_y") +
    geom_hline(yintercept = c(0:1), linewidth = .8, col = 8) +
    geom_hline(yintercept = .05, col = 8) +
    geom_step(stat = "ecdf", col = 4, linewidth = width, alpha = .8) +
    geom_vline(data = d0, mapping = aes(xintercept = err_data), col = 2, linewidth = 1, alpha = .5) +
    geom_hline(data = d0, mapping = aes(yintercept = p_value), col = 2, linetype = 2, alpha = .5) 
}

plot.cnaTest <- function(x, type = c("ecdf", "violins", "bars"), ...){
  type <- match.arg(type)
  switch(type, 
         ecdf = plot_ecdf(x, ...),
         violins = plot_violin(x, ...), 
         bars = plot_bars(x, ...))
}

