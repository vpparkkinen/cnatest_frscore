flipout <- function(data, outcome, proportion) {
  N <- nrow(data)
  out_col <- which(names(data) == outcome)
  range_o <- min(data[, out_col]):max(data[, out_col])
  n_to_flip <- round(N * proportion)
  n_to_flip <- if (n_to_flip == 0L) 1L else n_to_flip
  w_rows <- sample(1:N, n_to_flip)
  for(row in w_rows){
    ov <- data[row, out_col]
    not_ov <- range_o[range_o != ov]
    change_ov_to <- sample(not_ov, 1)
    data[row, out_col] <- change_ov_to
  }
  return(data)
}

any_submodel <- function(x, y) {
  if(is.null(x)) return(NA)
  as <- lapply(x, function(z) bis_submodel(z,y))
  as <- unlist(as)
  if(length(as) == 1 && is.na(as)){
    return(NA)
  } else {
    return(any(as))
  }
  
}


bis_submodel <- function(x,y){
  if (is.na(x) || is.null(x)){
    NA
  } else {
    is.submodel(x, y)
  }
}

bs_dat_create <- function(Nsets = 1e3, 
                          size = 30, 
                          varnum = 7,
                          type = c("cs", "fs"),
                          varnames = LETTERS[1:varnum]){
  type = match.arg(type)
  if (type == "fs"){
    c <- quote(runif(size, min = 0, max = 1))
  } 
  if (type == "cs"){
    c <- quote(rbinom(n = size, size = 1, prob = 0.5))
  }
  dsets <- vector("list", Nsets)
  for(i in 1:Nsets){
    dsets[[i]] <- data.frame(setNames(
      replicate(varnum, eval(c), simplify = FALSE), varnames))
  }
  return(dsets)
}


prevalence_compliant_noisify <- function(model, data, outcome, noiselevel){
  if((noiselevel * nrow(data)) %% 1 != 0L) warning("noiselevel is not a fraction of N")
  N <- nrow(data)
  n_noise <- round(noiselevel * N)
  cdat <- ct2df(selectCases(model, full.ct(data)))
  ndat <- data.frame(setdiff(ct2df(full.ct(data)), cdat))
  ndatsplit <- split(ndat, ndat[,outcome])
  datasplit <- split(data, data[,outcome])
  data_ps <- lapply(datasplit, function(x) round((nrow(x) / N) * n_noise))
  if(sum(unlist(data_ps)) == 0L){
    idx <- sample(1:length(datasplit),1)
    tc <- datasplit[[idx]][-sample(1:nrow(datasplit[[idx]]),1), ]
    tn <- ndatsplit[[idx]][sample(1:nrow(ndatsplit[[idx]]),1), ]
    temp_ndata <- tn
    temp_cdata <- datasplit
    temp_cdata[[idx]] <- tc
  } else {
    temp_cdata <- mapply(\(x,y){if (y == 0L) x else x[-sample(1:nrow(x), y),]}, 
                         datasplit, 
                         data_ps,
                         SIMPLIFY = FALSE)
    temp_ndata <- mapply(\(x,y){
      if (y == 0L) {x[0L,]} else {x[sample(1:nrow(x), y, replace = TRUE),]}
      }, 
                         ndatsplit, 
                         data_ps,
                         SIMPLIFY = FALSE)  
  }
  
  ndata_all <- if(class(temp_ndata) == "data.frame") {
    temp_ndata
    } else {do.call(rbind, temp_ndata)}
  cdata_all <- do.call(rbind, temp_cdata)
  out <- rbind(cdata_all, ndata_all)
  attr(out, "noise") <- ndata_all
  return(out)
}


prevalence_fixer <- function(data, outcome, prevalence, N){
  if((prevalence * N) %% 1 != 0L) warning("prevalence is not a fraction of N")
  if(class(substitute(outcome, parent.frame())) == "call"){o <- outcome} else{
    o <- substitute(outcome)  
  }
  
  o_idx <- eval(o, data, parent.frame())
  o_present <- data[o_idx,]
  opnr <- nrow(o_present)
  if (nrow(data) == N && prevalence == opnr / N){
    return(data)
  }
  o_absent <- data[!o_idx,]
  oanr <- nrow(o_absent)
  n_preval_rows <- round(N * prevalence)
  n_nonpreval_rows <- N - n_preval_rows
  if (opnr > n_preval_rows){
    prev_out <- o_present[sample(1:opnr, n_preval_rows),]
  } else {
    op_sample <- sample(1:opnr, n_preval_rows - opnr, replace = T)
    prev_out <- rbind(o_present, o_present[op_sample,])
  }
  if (oanr > n_nonpreval_rows){
    nonprev_out <- o_absent[sample(1:oanr, n_nonpreval_rows),]
  } else {
    oa_sample <- sample(1:oanr, n_nonpreval_rows - oanr, replace = T)
    nonprev_out <- rbind(o_absent, o_absent[oa_sample,])
  }
   out <- rbind(prev_out, nonprev_out)
   return(out)
}


makenoisy_asf <- function(model, 
                        data = ct2df(selectCases(model)), 
                        outcome, 
                        prevalence, 
                        noiselevel,
                        N = nrow(data)){
  oc_temp <- substitute(outcome)
  oc <- deparse(oc_temp)
  oc <- gsub(" ", "", strsplit(oc, "==")[[1]][1])
  prev_dat <- prevalence_fixer(data, oc_temp, prevalence, N)
  out <- prevalence_compliant_noisify(model, prev_dat, oc, noiselevel)
  return(out)
}
