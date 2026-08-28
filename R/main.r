source('./R/utilities.r')
source('./R/prn_generator.r')
source('./R/prob_fit.r')
source('./R/gofTest.chisq.r')


## given a set of observed data, fit it with probability distributions
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
