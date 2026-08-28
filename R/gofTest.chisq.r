## given a set of observed data, fit it with probability distributions
## a wrap for gofTest.chisqC and gofTest.chisqD.
## 
## returns a list
##    estimated parameters
##    chisq statistics
##    p.values
##    cut.points
##    n.obs
##    n.exp
distrEstimate <- function(x, 
                         distribution = NULL, 
                         paramList = NULL,
                         estimate.parameters = T,
                         alpha = 0.05,
                         return.best = T,
                         lower.tail = F) {
  
  x.raw <- x
  if(!is.numeric(x)) {
    message('data is coerced to numeric')
    x <- as.numeric(x)
  }
  
  bad.obs <- sum(is.na(x) | is.infinite(x))
  x <- x[!(is.na(x) | is.infinite(x))]
  if(length(x) < MIN_OBS) {
    stop(sprintf('Need at least %d data points for the estimate.\n\t%d supplied', MIN_OBS, length(x)))
  }

  
  if(!is.null(distribution)) {
    distrs <- distribution
  } else {
    if (!all(x == round(x))) {
      distrs <- DISTR.NAMES$distr.abbr[DISTR.NAMES$distr.type == 'continuous']
    } else if (mean(x==0) > 0.1) {
      distrs <- DISTR.NAMES$distr.abbr[DISTR.NAMES$distr.type == 'discrete']
    } else if (all(x == round(x)) && length(table(x))<10) {
      distrs <- DISTR.NAMES$distr.abbr[DISTR.NAMES$distr.type == 'discrete']
    } else {
      distrs <- DISTR.NAMES$distr.abbr
    }
  }
  
  chisq.results <- vector(mode = 'list', length = length(distrs))
  names(chisq.results) <- distrs
  for(distr_ in distrs){
    distr <- getDistrInfo(distr_)
    type <- distr$type
    
    fit_ <- tryCatch({
      if(type == 'discrete') {
        gofTest.chisqD(x = x, 
                       distribution = distr$name,
                       paramList = paramList,
                       estimate.parameters=estimate.parameters,
                       # expected.prob = expected.prob,
                       alpha = alpha,
                       lower.tail = lower.tail)
      } else {
        gofTest.chisqC(x = x, 
                       distribution = distr$name,
                       paramList = paramList,
                       estimate.parameters=estimate.parameters,
                       # expected.prob = expected.prob,
                       alpha = alpha,
                       lower.tail = lower.tail)
      }
    }, error = function(e) {
      message(sprintf('skipping %s: %s', distr_, conditionMessage(e)))
      NULL
    })
    chisq.results[[distr_]] <- fit_
  
  } 
   
  chisq.results[sapply(chisq.results, is.null)] <- NULL   
  if(return.best) {
    pvals <- lapply(chisq.results, function(chi) chi$p.value) %>% unlist()
    chival <- lapply(chisq.results, function(chi) chi$chisq.statistics)
    idx <- which(pvals == max(pvals))
  } else {
    idx <- 1:length(chisq.results)
  }

  result <- lapply(idx, function(i) {c(chisq.results[[i]], 
                                        list(data = x.raw, 
                                             n.obs = length(x.raw), 
                                             bad.obs = bad.obs))
                                    })
  
  return(result)
}

## for chisq goodness of fit test
gofTest.chisqD <- function(x,
                           distribution,
                           paramList = NULL,
                           estimate.parameters = T,
                           alpha = 0.05,
                           lower.tail = F){
  
  distr <- getDistrInfo(distribution)

  if(is.null(paramList)) {
    estimate.parameters <- T
  }

  if(estimate.parameters) {
    funcFit <- sprintf('x%sFit', distr$name)
    paramList <- do.call(funcFit, list(x = x, param.names = distr$param.names))
    if(!is.list(paramList) && (all(unlist(paramList) == -1))) {
      message(sprintf('data does not fit %s', distr$name))
      return(NULL)
    }
  } else {
    if(!all(names(paramList) %in% distr$param.names)) {
      invalid_params <- names(paramList)[!(names(paramList) %in% distr$params.names)] %>% 
        str_c(., collapse = ';')
      stop(sprintf('parameter %s not recognized.\nPlease check DISTR.NAMES for parameter names.',
                   invalid_params))
    }
  }
  
  if (is.character(distr$x.min)) {
    min.name <- distr$x.min 
    distr$x.min <- paramList[[min.name]]
  }
  
  if (is.character(distr$x.max)) {
    max.name <- distr$x.max
    distr$x.max <- paramList[[max.name]]
  }
  n <- length(x)
  max.obs <- max(x)
  min.obs <- min(x)
  
  Ev_range <- switch(distr$name,
                     'bern' = c(0,1),
                     'binom' = 0:paramList$size,
                     'geom' = 1:max.obs,
                     'nbinom' = paramList$size:max.obs,
                     'pois' = 0:max.obs)
  
  mfunc <- sprintf('xm%s', distr$name)
  args <- c(list(x = Ev_range[-length(Ev_range)]), paramList)
  prob_xs <- do.call(mfunc, args)
  prob_tail <- 1-sum(prob_xs)
  prob_xs <- c(prob_xs, prob_tail)
  Eii <- prob_xs * length(x)
  names(Eii) <- Ev_range
 
  Oii <- table(x)
  ## fill observation that are missing in Oii
  Oii <- sapply(names(Eii), function(nm) ifelse(is.na(Oii[nm]), 0, Oii[nm]))  
  
  if(!(length(Oii) == length(Eii))) {
    stop("Expected value and Observed values in different length")
  }
  
  ## combine lower and upper tails that are less than 5
  ## won't check values in the middle
  k <- length(Eii)
  if (Eii[1] < 5) {
    csum <- cumsum(Eii)
    idx1 <- max(which(csum < 5))
  } else {
    idx1 <- 1
  }
  
  if (Eii[k] < 5) {
    rcsum <- rev(cumsum(rev(Eii)))
    idx2 <- tryCatch(min(which(rcsum < 5)), error = k+1)
  } else {
    idx2 <- k
  }
  
  Eii <- c(sum(Eii[1:idx1]), Eii[(idx1+1):(idx2-1)], sum(Eii[idx2:k]))
  Oii <- c(sum(Oii[1:idx1]), Oii[(idx1+1):(idx2-1)], sum(Oii[idx2:k]))
    
  kk <- length(Eii)
  names(Eii)[1] <- ifelse(idx1 == 1, 
                          Ev_range[1],
                          sprintf('%s-', as.numeric(names(Eii)[2])-1)
                          )
  names(Eii)[kk] <- ifelse(idx2 == k, 
                           Ev_range[k],
                           sprintf('%s+', as.numeric(names(Eii)[kk-1])+1))
  #}
  names(Oii) <- names(Eii)
  
  n.param.est <- ifelse(estimate.parameters, length(paramList), 0)
  n.param.est <- ifelse(distr$name == 'bern', 0, n.param.est)
  
  df <- length(Eii) - 1 - n.param.est
  chisq <- sum((Oii-Eii)^2/Eii)
  
  if(df <=0 ) {
    message(sprintf('df = %d\ndata is not sufficuent for chisq test', df))
    chisq.a <- NA
    p.val <- 0
  } else {
    chisq.a <- qchisq(alpha, df = df, lower.tail = lower.tail)
    p.val <- pchisq(chisq, df = df, lower.tail = lower.tail)
  }
  
  return(list(distr.name = unique(c(distr$name2,distr$name)),
              parameters = unlist(paramList),
              n.valid.observes = length(x),
              chisq.statistics = chisq,
              p.value = p.val,
              df = df,
              alpha = alpha,
              chisq.critical = chisq.a,
              bins = names(Eii),
              number.observed = Oii,
              number.expected = Eii
  ))
}  

## for continuous distributions
gofTest.chisqC <- function(x, 
                           distribution,
                           paramList = NULL,
                           estimate.parameters = T,
                           alpha = 0.05,
                           lower.tail = F) {
  
  distr <- getDistrInfo(distribution)

  if(is.null(paramList)) {
    estimate.parameters <- T
  }
  
  if(estimate.parameters) {
    funcFit <- sprintf('x%sFit', distr$name)
    paramList <- do.call(funcFit, list(x = x, param.names = distr$param.names))
    if(!is.list(paramList) && (all(unlist(paramList) == -1))) {
      message(sprintf('data does not fit %s', distr$name))
      return(NULL)
    }
  } else {
    if(!all(names(paramList) %in% distr$param.names)) {
      invalid_params <- names(paramList)[!(names(paramList) %in% distr$params.names)] %>% 
        str_c(., collapse = ';')
      stop(sprintf('parameter %s not recognized.\nPlease check DISTR.NAMES for parameter names.',
                   invalid_params))
    }
  }
  
  if (is.character(distr$x.min)) {
    min.name <- distr$x.min 
    distr$x.min <- paramList[[min.name]]
  }
  
  if (is.character(distr$x.max)) {
    max.name <- distr$x.max
    distr$x.max <- paramList[[max.name]]
  }
  
  n <- length(x)
  k <- min(100, ceiling(2*n^(1/3))-1)
  while(T) {
    Ei <- rep(n/k, k)
    if (k <= 2) break
    if(any(Ei < 5)) {
      k <- k-1
    } else {
      break
    }
  }
  ps <- (1:k)/k
  
  qfunc <- sprintf('q%s', distr$name)
  
  if (distr$name == 'weibull') {
    scale <- 1/paramList$lambda
    names(paramList) <- c('shape', 'scale')
    paramList$scale <- scale
  } else if (distr$name == 'erlang') {
    ## special case of gamma when shape is integer
    if(!(round(paramList$shape) == paramList$shape)) {
      stop('shape should be an integer for erlang distribution')
    }
    qfunc <- 'qgamma'
 # )  
  }
  
  cut.points <- sapply(ps, function(p) {
                            args <- c(list(p = p), paramList)
                            do.call(qfunc, args)})
  
  ## in some extreme cases, the tail of cut.points will be duplicated
  ## merge Ei of duplicated cells if have duplicated cut.points
  if(any(duplicated(cut.points))) {
    dup.idx <- duplicated(cut.points)
    Ei[dup.idx[1]-1] <- sum(Ei[(dup.idx[1]-1):max(dup.idx)])
    Ei <- Ei[1:(dup.idx[1]-1)]
  }
  
  cut.points <- unique(c(distr$x.min, cut.points))
  
  Oi <- cut(x, breaks = cut.points) %>% table()
  
  chisq <- sum((Oi-Ei)^2/Ei)
  n.param.est <- ifelse(estimate.parameters, length(paramList), 0)
  n.param.est <- ifelse(distr$name == 'bern', 0, n.param.est)
  df <- length(Ei) - 1 - n.param.est
  
  if(df <=0 ) {
    message(sprintf('df = %d\ndata is not sufficuent for chisq test', df))
    chisq.a <- NA
    p.val <- 0
  } else {
    chisq.a <- qchisq(alpha, df = df, lower.tail = lower.tail)
    p.val <- pchisq(chisq, df = df, lower.tail = lower.tail)
  }
  
  if (distr$name == 'weibull') {
    lambda <- 1/paramList$scale
    names(paramList) <- c('shape', 'lambda')
    paramList$lambda <- lambda
  }
  
  return(list(distr.name = unique(c(distr$name2,distr$name)),
              parameters = unlist(paramList),
              n.valid.observes = length(x),
              chisq.statistics = chisq,
              p.value = p.val,
              df = df,
              alpha = alpha,
              chisq.critical = chisq.a,
              bins = cut.points,
              number.observed = unname(Oi),
              number.expected = Ei
              ))
}
