## this file contains:
##    global variables,
##    PMF functions for discrete distributions,
##    some help functions

 
getDensity <- function(distr_abbr, x, params) {
  distr_type <- getDistrInfo(distr_abbr)$type
  if(distr_type == 'discrete') x = round(x)
  dfunc <- ifelse(distr_type == 'discrete',
                  sprintf('xm%s', distr_abbr),
                  sprintf('d%s', distr_abbr))
  if(distr_abbr == 'weibull' && is.null(params$scale)) {
    params$scale <- 1/params$lambda
    params$lambda <- NULL
  }
  if(distr_abbr == 'erlang') {
    dfunc <- 'dgamma'
  }
  
  den <- tryCatch(do.call(dfunc, c(list(x = x), params)), 
                  error = function(e) NULL)
  return(den)
  
}

getBreaks <- function(data, min.bins = 20, max.bins = 100) {
  if (all(data == round(data))) {
    return(seq(min(data) - 0.5, max(data) + 0.5, by = 1))
    
  } else {
    nb <- sqrt(length(data)) %>% 
      round() %>% 
      min(., max.bins) %>% 
      max(., min.bins)
    return(seq(min(data), max(data), length.out = nb + 1))
  }
}

# Helper to format multiple datasets of unequal length into a tab-separated table string
list2tsv <- function(ds_list) {
  if (length(ds_list) == 0) return("")
  
  max_len <- max(sapply(ds_list, length))
  df <- lapply(ds_list, `length<-`, max_len) %>% as.data.frame()
  
  tc <- textConnection("out_str", "w", local = TRUE)
  write.table(df, file = tc, sep = "\t", row.names = FALSE, col.names = TRUE, na = "", quote = FALSE)
  close(tc)
  
  paste(out_str, collapse = "\n")
}

makeSafeId <- function(id, algo = 'md5') {
    paste0("ds_", digest::digest(id, algo = algo))
}

## format parameter list for better display
param2str <- function(param, collapse = '\t', pad = 10, indent = 4) {
  if(is.list(param)) param <- unlist(param)
  if(is.null(names(param))) names(param) <- str_c('parameter', 1:length(param))
  
  pad <- max(pad, nchar(names(param)))
  pad <- ifelse(collapse == "\n", pad, 0)
  s <- sapply(names(param), function(pname) {
    sprintf("%-*s%-*s %.4f",indent," ", pad, paste0(pname, ":"), param[[pname]])
  }) %>% str_c(., collapse = collapse)

  return(s)
}


## to format fit results for display
## input: 
##    fit.res: a list of fitting result from distrEstimate
## return: 
##    a vector of strings
##    each element is a formatted string of one fitting result
fit2str <- function(fit.res, collapse = '\n', pad = 12) {
  dname <- fit.res$distr.name %>% rev() %>% `[`(1) %>% 
    sprintf("%s: %s", "Distribution", .)
  pval <- fit.res$p.value #%>% 
  pval_str <- ifelse(pval < 1e-4, sprintf("%.4e", pval), sprintf("%.4f",pval)) %>% 
    sprintf("  %-*s %s", pad, "p-value:", .)
  chisq <- fit.res$chisq.statistics %>% round(4) %>% 
    sprintf("  %-*s %s", pad, "Chisq:", .)
  param_str <- param2str(fit.res$parameters, collapse = "\n", pad = pad-2) %>% 
    sprintf("  %s:\n%s", "Parameters", .)
  
  str <- c(dname, pval_str, chisq, param_str) %>% str_c(., collapse = "\n")
 return(str) 
}

## to get information for a given distribution
getDistrInfo <- function(distribution) {
  distr.idx <- which(distribution[1] == DISTR.NAMES[,2] | distribution[1] == DISTR.NAMES[, 1])
  if (!(length(distr.idx) == 1)) {
    stop(sprintf('Ambigous distribution name.\n Check `DISTR.NAMES` for distributions supported.', 
                 distribution))
  }
  
  n.params <- DISTR.NAMES$n.params[distr.idx]
  param.names <- DISTR.NAMES[distr.idx, str_c('param', 1:n.params, sep = '.')] %>% 
    unlist() %>% unname()
  param.default <- DISTR.NAMES[distr.idx, str_c('default', 1:n.params, sep = '.')] %>% 
    unlist() %>% unname() %>% as.list() %>% `names<-`(param.names)
  
  ## get lower and upper boundary of supported range 
  suppressWarnings( {
    x.min <- DISTR.NAMES[distr.idx, 'support.min'] %>% as.numeric
    x.max <- DISTR.NAMES[distr.idx, 'support.max'] %>% as.numeric
  })
  
  if(is.na(x.min)) x.min <- DISTR.NAMES[distr.idx, 'support.min']
  if(is.na(x.max)) x.max <- DISTR.NAMES[distr.idx, 'support.max']
  
  distrInfo <- list(idx = distr.idx,
                    name = DISTR.NAMES$distr.abbr[distr.idx],
                    name2 = DISTR.NAMES$distr.name[distr.idx],
                    type = DISTR.NAMES$distr.type[distr.idx],
                    n.params = DISTR.NAMES$n.params[distr.idx],
                    param.names = param.names,
                    param.defaults = param.default,
                    x.min = x.min,
                    x.max = x.max
  )
  
  return(distrInfo)
}

## to generate random parameters for a given distribution
makeParams <- function(distribution = 'unif') {
  distr <- getDistrInfo(distribution)
  dnm <- distr$name
  params <- rep(NA, distr$n.params)
  
  if (dnm == 'beta') { 
    u <- xunif(2)
    s1 <- case_when(u[1]<0.3 ~ xunif(1, 0.1, 1),
                    u[1]<0.8 ~ xunif(1, 1, 20),
                   .default = 1)
    s2 <- case_when(u[2] < 0.3 ~ xunif(1, 0.1, 1),
                    u[2] < 0.8 ~ xunif(1, 1, 20), 
                    .default = 1)
    params <- c(s1, s2)
  } 
  else if (dnm == 'cauchy') {
    u <- xunif(1)
    loc <- xnorm(1, mean = 0, sd = 20)
    s2<- case_when(u < 0.3 ~ xunif(1, 0.1, 1),
                   u < 0.8 ~ xunif(1, 1, 20), 
                   .default = 1)
    
    params <- c(loc, s2)
  }
  else if (dnm == 'exp') {
    params <- xunif(1, min = 0.1, max = 20)
  }
  else if(dnm == 'gamma') {
    u <- xunif(2)
    s1 <- case_when(u[1]<0.3 ~ xunif(1, 0.1, 1),
                    u[1]<0.8 ~ xunif(1, 1, 20),
                    .default = 1)
    s2 <- case_when(u[2] < 0.3 ~ xunif(1, 0.1, 1),
                    u[2] < 0.8 ~ xunif(1, 1, 20), 
                    .default = 1)
    params <- c(s1, s2)
  }
  else if(dnm == 'norm') {
    params[1] <- xnorm(1, mean = 0, sd = 20)
    params[2] <- xunif(1, min = 1, max = 20)
  }
  else if(dnm == 'tri') {
    abc <- xtri(3, -30,0, 30) %>% sort()
    params <- abc
  }
  else if(dnm == 'unif') {
    ab <- xunif(2, -20, 20) %>% sort()
    params <- ab
  } 
  else if(dnm == 'weibull') {
    u <- xunif(2)
    k <- case_when(u[1]<0.3 ~ xunif(1, 0.1, 1),
                   u[1]<0.8 ~ xunif(1, 1, 20),
                   .default = 1)
    lamb <- case_when(u[2] < 0.3 ~ xunif(1, 0.1, 1),
                      u[2] < 0.8 ~ xunif(1, 1, 20), 
                      .default = 1)
    params <- c(k, lamb)
  }
  else if(dnm == 'erlang') {
    u <- xunif(2)
    k <- case_when(u[1]<0.3 ~ 1,
                   .default = xrandint(1, 1, 20))
    lamb <- case_when(u[2] < 0.3 ~ xunif(1, 0.1, 1),
                      u[2] < 0.8 ~ xunif(1, 1, 20), 
                      .default = 1)
    params <- c(k, lamb)
  } 
  
  ## discrete distributions ============
  else if(dnm == 'dunif') {
    params <- xrandint(-30, 30) %>% sort()
  }
  else if(dnm == 'bern') {
    params <- xunif(1)
  }
  else if(dnm == 'binom') {
    params[1] <- xrandint(1, 6, 30)
    params[2] <- xunif(1, 0.2, 0.8)
  }
  else if(dnm == 'geom') {
    params <- xunif(1, 0.2, 0.8)
  }
  else if(dnm == 'nbinom') {
    params[1] <- xrandint(1, 1, 10)
    params[2] <- xunif(1, 0.2, 0.8)
  }
  else if(dnm == 'pois') {
    u <- xunif(1, 0.1, 1)
    lamb <- case_when(u < 0.3 ~ xunif(1),
                      u < 0.8 ~ xunif(1, 1, 20), 
                      .default = 1)
    params <- lamb
  }
  names(params) <- distr$param.names
  return(as.list(params))
}

## to validate if the range of data x fits a given distribution
## x: numeric vector
## distr: target distribution; see DISTRIBUTION.NAMES for supported distributions
## 
## return: 
dataRangeValidation <- function(x, distr, params = NULL) {
  
  if (!all(is.numeric(x))) stop('x has to be numeric.')
  
  LB <- F
  UP <- F
  
  xrange <- range(x)
  distrInfo <- getDistrInfo(distr)
  d.min <- distrInfo$x.min
  d.max <- distrInfo$x.max
  
  if(is.null(params)) {
    if(!is.numeric(as.numeric(d.min))) {
      LB <- T
    } else {
      LB <- ifelse(as.numeric(d.min) < xrange[1], T, F)
    }
    
    if(!is.numeric(d.max)) {
      UP <- F
    } else {
      UP <- ifelse(d.max > xrange[2], T, F)
    }
  }
  
  if(!is.numeric(d.min)) {
    if(is.null(params)) {
      LB <- T
    } else {
      
    }
  }
  
  
  
}


validate_data_in_distribution_range <- function(data, distr_abbr, dist_table, params = list()) {
  # 1. Type check: Ensure data is a numeric vector
  if (!is.numeric(data)) {
    stop("The 'data' argument must be a numeric vector.")
  }
  
  # 2. Lookup distribution specifications
  dist_row <- dist_table[dist_table$distr.abbr == distr_abbr, ]
  if (nrow(dist_row) == 0) {
    stop(paste0("Distribution '", distr_abbr, "' not found in the reference table."))
  }
  
  # Map default parameter values
  param_values <- list()
  for (i in 1:3) {
    p_name <- as.character(dist_row[[paste0("param.", i)]])
    p_def <- dist_row[[paste0("default.", i)]]
    
    if (!is.na(p_name) && p_name != "") {
      param_values[[p_name]] <- p_def
    }
  }
  
  # Override defaults with user-supplied custom parameters
  if (length(params) > 0) {
    for (name in names(params)) {
      param_values[[name]] <- params[[name]]
    }
  }
  
  # 3. Helper to resolve dynamic support bounds (e.g., 'min', 'max', 'size', 'Inf', '-Inf')
  resolve_bound <- function(bound_val) {
    bound_str <- trimws(as.character(bound_val))
    
    if (is.na(bound_str) || bound_str == "") return(NA)
    if (bound_str == "-Inf") return(-Inf)
    if (bound_str == "Inf") return(Inf)
    
    # If the bound refers to a parameter name (e.g., 'size', 'min', 'max')
    if (bound_str %in% names(param_values)) {
      val <- param_values[[bound_str]]
      if (is.na(val) || is.null(val)) {
        stop(paste0("Parameter '", bound_str, "' required for support bound is not defined."))
      }
      return(as.numeric(val))
    }
    
    return(as.numeric(bound_str))
  }
  
  min_val <- resolve_bound(dist_row[["support.min"]])
  max_val <- resolve_bound(dist_row[["support.max"]])
  
  # 4. Check boundaries
  out_of_bounds_mask <- (data < min_val) | (data > max_val)
  
  list(
    is_valid = !any(out_of_bounds_mask, na.rm = TRUE),
    distribution = distr_abbr,
    support_min = min_val,
    support_max = max_val,
    out_of_bounds_count = sum(out_of_bounds_mask, na.rm = TRUE),
    total_elements = length(data)
  )
}
## === continuous PDF functions ===========================================
## f(x) = a*r*(a*x)^(r-1)*exp(-(a*x)^r), r= shape, a = lambda
xdweibull <- function(x, shape = 1, lambda = 1, scale = NULL) {
  if(is.null(lambda) && !is.null(scale)) lambda = 1/scale
  if(!is.numeric(shape) || shape < 0) stop('shape cannot be negative')
  if(!is.numeric(lambda) || lambda < 0) stop('lambda cannot be negative.')
  a <- lambda
  r <- shape
  
  fx <- a*r*(a*x)^(r-1)*exp(-(a*x)^r)
  return(fx)
  
}

## === discrete PMF functions ==============================================
## bernoulli pmf
## x = {0,1}
xmbern <- function(x, prob = 0.5) {
  if (!all(x %in% c(0,1))) {
    stop('x have to be 0 or 1')
  }
  
  fx <- ifelse(x, prob, 1-prob)
  return(fx)
}

## binomal pmf
## x = {0,1,2,3,...}
xmbinom <- function(x, size, prob = 0.5) {
  if(!is.numeric(x) || !all(x == round(x))|| any(x<0)) {
    stop('x have to be non-negative integers')
  }
  pmf <- choose(size, x) * prob^x * (1-prob)^(size-x)
  pmf <- ifelse(is.na(pmf), 0, pmf)
  return(pmf)
}

## geom pmf
## number of trials till success
## x = {1,2,3,...}
xmgeom <- function(x, prob = 0.5) {
  if(!is.numeric(x) || !all(x == round(x))|| any(x<1)) {
    stop('x have to be positive integers')
  }
  prob * (1-prob)^(x-1)
}

## negative binomial pmf
## number of trials till r-th success.
## x ~ c(1,2,3,...)
xmnbinom <- function(x, size, prob = 0.5) {
  if(!is.numeric(x) || !all(x == round(x))|| any(x<size)) {
    stop('x have to be positive integers')
  }
  pmf <- choose(x-1, size-1) * prob^size * (1-prob)^(x-size)
  pmf <- ifelse(is.na(pmf), 0, pmf)
  return(pmf)
}

## pois distribution pmf
## x = c(0, 1, 2, 3, ...)
xmpois <- function(x, lambda) {
  if(!is.numeric(x) || !all(x == round(x))|| any(x<0)) {
    stop('x have to be non-negative integers')
  }
  lambda^x * exp(-lambda) / factorial(x)
}


