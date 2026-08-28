
## Norm distribution estimators, continuous distribution ==================================================
## x ~ (-Inf, Inf)
xnormFit <- function(x, method = c('mom', 'mle'), param.names = c('mean', 'sd')) {
  cat('fitting normal distribution...\n') 
  if(!is.numeric(x)) {
    message('non numeric numbers.')
    return(-1)
  }
   method = tolower(method[1])
   if (!method %in% c('mom', 'mle')) stop(sprintf('method not supported: %s', method))
   
   params <- switch(method,
                    'mom' = xnormFitMoM(x),
                    'mle' = xnormFitMLE(x)) %>% 
     setNames(., param.names) %>% 
     as.list()
   
   return(params)
   
}

## MLE estimate of norm distribution 
xnormFitMLE <- function(x, param.names = c('mean', 'sd')){
  cat('fitting normal distribution...\n') 
  if(!is.numeric(x)) {
    message('non numeric numbers.')
    return(-1)
  }
  
  L <- function(param, data) {
    mean <- param[1]
    sd <- param[2]
    nllh <- 0.5 * sum((data-mean)^2)/sd^2 + 0.5 * length(data) * log(2*pi*sd^2)
    return(nllh)
  }
  
  result <- optim(c(0,1), fn = L, data = x, lower = c(-Inf, 1e-6), method = 'L-BFGS-B')
  return(result$par)
}

## MOM estimate of norm distribution 
## x ~ (-Inf, Inf)
xnormFitMoM <- function(x, param.names = c('mean', 'sd')){
  cat('fitting normal distribution...\n') 
  if(!is.numeric(x)) {
    message('non numeric numbers.')
    return(-1)
  }
  
  n <- length(x)
  mu_hat <- mean(x)
  sigma_hat <- sqrt(sum((x-mu_hat)^2)/n)

  return(c(mu_hat, sigma_hat))
}


## MLE estimate of uniform distribution
## x ~ (-Inf, Inf)
xunifFit <- function(x, param.names = c('min', 'max')) {
  cat('fitting uniform distribution...\n')   
  if(!is.numeric(x)) {
    message('non numeric numbers.')
    return(-1)
  }
  
  n <- length(x)
  params = list(min = min(x), max = max(x))
  
  return(params)
}

## MLE estimate of exponential distribution
## x ~ [0, Inf)
xexpFit <- function(x, param.names = c('rate')){
  cat('fitting exponential distribution...\n')
  if(!is.numeric(x) || any(x < 0)) {
    message('Not exponential distribution\n\tNo estimate returned...')
    return(-1)
  }
  
  ## log likelihood
  L <- function(lambda) { 
    n <- length(x)
    llh <- n * log(lambda) - lambda * sum(x) 
    return(llh)
  } 
  
  lambda <- optimize(L, 
                     interval = c(1e-6, 1e23), 
                     maximum = TRUE)$maximum 
  
  params <- list(rate = lambda)
  
  return(params)
  
}

## MLE estimate of triangular distribution
## x ~ (-Inf, Inf)
xtriFit <- function(x, param.names = c('min', 'mode', 'max')) {
  cat('fitting triangular distribution...\n') 
  if(!is.numeric(x)) {
    message('non numeric numbers.')
    return(-1)
  }
  
  ## method of moments to estimate initial parameters
  a0 <- min(x)
  c0 <- max(x)
  b0 <- 3 * mean(x) - a0 - c0
  if(b0 < a0) b0 <- a0 + 1e-6
  if(b0 > c0) b0 <- c0 - 1e-6
  
  ## logliklihood for more accurate estimate.
  L <- function(param, data) {
    a <- param[1]
    b <- param[2] ## b as mode
    c <- param[3]
    if(!(a<=b && b<=c)) return(1e23)
    
    idx_left <- which(data < b)
    idx_right <- which(data >= b)
    
    rng_left <- b-a + 1e-6
    rng_right <- c-b + 1e-6
    rng <- c-a + 1e-6
    
    llh_1  <- -sum(log(2) + log(data[idx_left] - a+1e-6) - log(rng_left) - log(rng))
    llh_2 <- -sum(log(2) + log(c - data[idx_right]+1e-6) - log(rng_right) - log(rng))
    return(llh_1 + llh_2)
    
  }
  
  result <- optim(c(a0,b0,c0), 
                  fn = L, 
                  data = x, 
                  method = 'L-BFGS-B',
                  lower = c(-Inf, min(x), max(x)+1e-6),
                  upper = c(min(x)-1e-6, max(x), Inf))
  params <- result$par
  names(params) <- param.names
  return(as.list(params))
}

## Estimate of erlang distribution parameters
xerlangFit <- function(x, method = c('mom', 'mle'), 
                       param.names = c('shape', 'rate'), 
                       shape.upper = 500) {
  cat('fitting erlang distribution...\n') 
  if(!is.numeric(x) || any(x<0)) {
    message('not erlang distribution\n\tNo estimate returned...')
    return(-1)
  }
  
  method = method[1]
  if (!(method %in% c('mom', 'mle'))) {
    stop(sprintf('method not supported: %s', method))
  }
  
  params <- switch(method,
                   'mom' = xerlangFitMoM(x),
                   'mle' = xerlangFitMLE(x, shape.upper)) %>% 
    setNames(., param.names) %>% 
    as.list()
  
  return(params)
  }

## MLE estimate of erlang distribution
## x  ~ [0, Inf)
xerlangFitMLE <- function(x, shape.upper = 500) {
  ## MOM of lambda
  rate0 <- mean(x)/var(x)
  
  L <- function(rate, k, data) {
    n <- length(data)
    llh <- n * k * log(rate) + (k-1)*sum(log(data+1e-6)) - rate*sum(data) - n*lgamma(k)
    return(-llh) 
  }
  
  ## loop through integers to find best shape parameter
  ks <- 1:11
  while(T) {
    objs <- rep(0, length(ks))
    rates <- rep(0, length(ks))
    for (i in 1:length(ks)) {
      k <- ks[i]
      para <- try(optim(rate0,fn=L,data = x,k=k, lower = 1e-6,upper = 1e6, method = 'L-BFGS-B'))
      if (inherits(para, "try-error")) {
        return(-1)
      } else {
        rates[i] <- para$par
        objs[i] <- para$value
      }

    }
    
    ## if objective functon minimize at the largest k,
    ## assume it will keep rising. 
    ## set one overlap value; may not be necessary
    min.idx <- which(objs == min(objs)[1])
    if(min.idx == length(ks)) {
      if (max(ks) > shape.upper) {
        message('shape excess upper bound.\nreturned values are not necessary optimal')
        k_hat <- ks[min.idx]
        rate_hat <- rates[min.idx]
        break
      }
      ks <- ks + 10
      ks <- ks[ks<shape.upper]
    } else {
      rate_hat <- rates[min.idx]
      k_hat <- ks[min.idx]
      break
    }
    
  } ## end of while loop
  
  return(c(k_hat, rate_hat))
}

## MOM estimate of erlang distribution
## x  ~ [0, Inf)
xerlangFitMoM <- function(x){
  #n <- length(x)
  mu <- mean(x)
  #sigma2 <- sum((x-mu)^2)/n
  sv <- var(x)
  rate <- mu/sv
  k <- round(mu * rate)
  
  return(c(k, rate))
  
}

## MLE estimate of weibull distribution
## x ~ [0, Inf)
xweibullFit <- function(x, param.names = c('shape', 'lambda')) {
  cat('fitting weibull distribution...\n')
  if(!is.numeric(x) || any(x<0)) {
    message('Not weibull distribution\n\tNo estimate returned...')
    return(-1)
  }
  
  nL <- function(params, data) {
    ## log likelihood function
    ## return negative log-likihood for optimization.
    
    n <- length(x)
    r <- params[1]
    lambda <- params[2]
    llh <- n*log(r) + n*r*log(lambda) + (r-1)*sum(log(data+1e-6)) - (lambda^r)*sum(data^r)
    return(-llh)
  }
  result <- optim(c(1,1), fn = nL, data = x, lower = c(1e-6, 1e-6), method = 'L-BFGS-B')
  params <- result$par
  names(params) <- param.names
  
  return(as.list(params))
}

## MLE estimate of gamma distribution
## x ~ [0, Inf)
xgammaFit <- function(x, param.names = c('shape', 'rate')) {
  cat('fitting gamma distribution...\n') 
  if(!is.numeric(x) || any(x<0)) {
    message('Not gamma distribution\n\tNo estimate returned...')
    return(-1)
  }
  n <- length(x)
  
  nL <- function(params, data) {
    ## log likihood function, return negative for minimization
    r <- params[1]
    lambda <- params[2]
    
    llh <- n * r * log(lambda) + (r-1) * sum(log(data)) - lambda * sum(data) - n * log(gamma(r))
    return(-llh)
  }
  
  result <- optim(c(1,1), fn = nL, data = x, method = 'L-BFGS-B', lower = c(1e-6, 1e-6))
  params = result$par
  names(params) <- param.names
  
  return(as.list(params))
}

## MOM estimate of beta distribution
## x ~ [0,1]
xbetaFit <- function(x,param.names = c('shape1', 'shape2')) {
  cat('fitting beta distribution...\n') 
  if(!is.numeric(x) || any(x>1) || any(x<0)) {
    message('Not beta distribution\n\tNo estimate returned...')
    return(-1)
  }
  
  n <- length(x)
  mu <- mean(x)
  s2 <- sum((x-mu)^2)/(n-1)
  
  beta <- (1-mu)^2*mu/s2 - 1 + mu
  alpha <- beta * mu/(1-mu)
  params <- c(alpha, beta)
  names(params) <- param.names
  return(as.list(params))
}

## MLE estimate of cauchy distribution
## x ~ (-Inf, Inf)
xcauchyFit <- function(x, param.names = c('location', 'scale')) {
  cat('fitting cauchy distribution...\n') 
  if(!is.numeric(x)) {
    message('Not numeric...')
    return(-1)
  }
  
  nL <- function(param, data) {
    x0 <- param[1]
    gamma <- param[2]
    n <- length(data)
    if(gamma<=0) return(Inf)
    
    nllh <- n*log(pi * gamma) + sum(log(1+((data-x0)/gamma)^2))
    return(nllh)
  }
  
  x0_start <- median(x)
  scale_start <- IQR(x) / 2
  result <- optim(c(x0_start, scale_start), 
                  fn = nL, 
                  data = x, 
                  method = 'L-BFGS-B',
                  lower = c(-Inf, 1e-6))
  
  params <- result$par
  names(params) <- param.names
  
  return(as.list(params))
}


## discrete distributions ======================================================

## MLE estimate of bernoulli distribution
## x in c(0, 1)
xbernFit <- function(x, param.names = 'prob'){
  cat('fitting bernoulli distribution...\n') 
  if(!all(unique(x) %in% c(0,1))) {
    message('Not bernoulli distribution\n\tNo estimation returned...')
    return(-1)
  }

  return(list(prob = mean(x)))
  
}

## Estimate of binomial distribution
## x in c(0, 1, 2, ...)
xbinomFit <- function(x, 
                      method = c('mom', 'mle'), 
                      param.names = c('size', 'prob'), 
                      size.upper = 1000) {
  cat('fitting binomial distribution...\n') 
  if(!(is.numeric(x))||(!all(x == round(x))) || any(x<0)) {
    message('Not geometric distribution\n\tNo estimation returned...')
    return(-1)
  }
  
  method <- tolower(method[1])
  if (!method %in% c('mom', 'mle')) {
    stop(sprintf('method not supported: %s', method))
  }
  
  params <- switch(method,
                   'mom' = xbinomFitMoM(x),
                   'mle' = xbinomFitMLE(x, size.upper)) %>% 
    setNames(., param.names) %>% 
    as.list()
  
  return(params)
}

## MoM estimate of binomial distribution
xbinomFitMoM <- function(x) {
  s_mean <- mean(x)
  s_var <- var(x)
  prob <- 1 - s_var/s_mean
  size <- round(s_mean/prob)
  
  return(c(size, prob))
}

## MLE estimate of binomial distribution
xbinomFitMLE <- function(x, size.upper = 1000){

  L <- function(prob, size, data) { 
    n <- length(data)
    llh <- sum(lchoose(size, data) + data * log(prob) + (size-data)*log(1-prob))
    return(-llh)
  } 
  
  ## loop through integers to find best shape parameter
  ## since k starts at max(x), it should be pretty close to the optimal size
  ks <- (0:10) + max(x)
  while(T) {
    objs <- rep(0, length(ks))
    probs <- rep(0, length(ks))
    for (i in 1:length(ks)) {
      k <- ks[i]
      result <- optim(0.5, L, size = k, data = x, 
                      lower = 1e-6, upper = 1-1e-6,  
                      method = 'L-BFGS-B') 
      probs[i] <- result$par
      objs[i] <- result$value
    }
    
    min.idx <- which(objs == min(objs)[1])
    if(min.idx == length(ks)) {
        
        if (max(ks) >= size.upper) {
          ## return estimates if upper size reached.
          message('shape excess upper bound.\nreturned values are not necessary optimal')
          k_hat <- ks[min.idx]
          prob_hat <- probs[min.idx]
          break
        }
        
        ks <- ks + 10
        ks <- ks[ks <= size.upper]
    } else {
        prob_hat <- probs[min.idx]
        k_hat <- ks[min.idx]
        break
    }

    
  } ## end of while loop
  
  #params <- c(k_hat, prob_hat)
  #names(params) <- param.names
  return(c(k_hat, prob_hat))
  
}


## MLE estimate of geometric distribution
## x in c(1, 2, 3, ...)
xgeomFit <- function(x, param.names = 'prob'){
  cat('fitting geometric distribution...\n') 
  if(!(is.numeric(x))||(!all(x == round(x))) || any(x<=0)) {
    message('not geom distribution...')
    return(-1)
  }
  
  if(all(x == 1)) {
    return(list(prob = 1))
  }
  
  L <- function(p, data) { 
    p <- p[1]
    n <- length(data)
    llh <- (sum(data) - n) * log(1-p) + n * log(p)
    return(-llh)
  }
  result <- optim(0.5, L, data = x, method = 'L-BFGS-B', lower = 1e-6, upper = 1 - 1e-6)
  params <- result$par
  names(params) <- param.names
  
  return(as.list(params))
  
}


## Estimate of negative binomial distribution
xnbinomFit <- function(x, method=c('mom', 'mle'), param.names = c('size', 'prob')) {
  cat('fitting negative binomial distribution...\n') 
  if(!(is.numeric(x)) || (!all(x == round(x))) || any(x<0)) {
    message('Not negative binomial distribution\n\tNo estimate returned...')
    return(-1)
  }
  
  method <- method[1]
  if(!(method %in% c('mom', 'mle'))) {
    stop(sprintf('method not supported: %s', method))
  }
  
  params <- switch(method,
                   'mom'=xnbinomFitMoM(x),
                   'mle'=xnbinomFitMLE(x)) %>% 
    setNames(., param.names) %>% 
    as.list()
  return(params)
  
}

## MoM estimate of negative binomial distribution
## modeled as # of total trials until r-th (size) success 
## prob: prob of success.
## x in c(1, 2, 3, ...)
xnbinomFitMoM <- function(x) {
  s_mean <- mean(x)
  s_var <- var(x)
  
  prob <- s_mean/(s_var + s_mean)
  size <- round(prob * s_mean)
  return(c(size, prob))
}

## MLE estimate of negative binomial distribution
## modeled as # of total trials until r-th (size) success 
## prob: prob of success.
## x in c(1, 2, 3, ...)
## works but slow due to looping through size parameter
xnbinomFitMLE <- function(x, param.names = c('size', 'prob')){
  L <- function(p, size, data) {
    if( p<0 || p > 1) return(Inf)
    if(size <= 0 || size > min(data)) return (Inf)
    llh <- sum(lchoose(data-1, size-1) + (data-size) * log(1-p) + size * log(p))
    return(-llh)
  }
  
  
  ks <- 1:min(11, min(x))
  while(T) {
    objs <- rep(0, length(ks))
    probs <- rep(0, length(ks))
    for (i in 1:length(ks)) {
      k <- ks[i]
      result <- optim(0.5, L, size = k, data = x, 
                      lower = 1e-6, upper = 1-1e-6,  
                      method = 'L-BFGS-B') 
      probs[i] <- result$par
      objs[i] <- result$value
    }
    
    min.idx <- which(objs == min(objs)[1])
    if(min.idx == length(ks)) {
        if (max(ks) >= min(x) ) {
          message('shape excess upper bound.\nreturned values are not necessary optimal')
          k_hat <- ks[min.idx]
          prob_hat <- probs[min.idx]
          break
        }
      
        ks <- ks + 10
        ks <- ks[ks <= min(x)]
    } else {
      prob_hat <- probs[min.idx]
      k_hat <- ks[min.idx]
      break
    }
    
  } ## end of while loop
  
  return(c(k_hat, prob_hat))
}

## MLE estimate of poisson distribution
## x in c(0, 1, 2, 3, ...)
xpoisFit <- function(x, param.names = 'lambda',  test='chisq'){
  cat('fitting poisson distribution...\n') 
  if(!(is.numeric(x)) || (!all(x == round(x))) || any(x<0)) {
    message('Not poisson distribution\n\tNo estimate returned...')
    return(-1)
  }
  
  L <- function(lamb,data) {
    n <- length(data)
    llh <- log(lamb) * sum(data) - n*lamb - sum(lfactorial(data))
    return(-llh)
  }
  result <- optim(1, L, data = x, lower = 1e-6, method = "L-BFGS-B")
  params <- result$par
  names(params) <- param.names
  return(as.list(params))
  
}














