
## to set global variable SEED
## if seed not supplied,  reset 
set_seed <- function(seed = NULL) {
  if(!(is.numeric(seed) || is.null(seed))) stop('seed has to be a number.')
  if(!is.null(seed)) SEED <<- floor(seed)
  else SEED <<- NULL
}

## make random seed from system time
## and the memory address of the time variable
.makeSeed <- function(){
  time <- as.numeric(Sys.time())
  addr <- obj_address(time) %>% 
    substr(13, 18) %>% 
    strtoi(base = 16L)
  
  return(time %% addr)
}

## Desert island algorithm
## get integers from 0 to 2^31 - 1
## the first number is always small when use a small seed
## first number is discarded
.desertIsland <- function(n, seed = SEED) {
  if(is.null(seed)) {
    seed <- .makeSeed()
  } else {
    seed = as.integer(seed)
    if (is.na(seed)) stop ("seed has to be an integer")
  }
  
  n <- n+5
  a <- 16807
  c <- 0
  
  X <- vector('integer', n)
  if (seed == 0) {
    X[1] <- (seed * a + 127773) %% M
  } else {
    X[1] <- seed
  }

  for (i in 2:n) {
    K <- floor(X[i-1] * a * M_INV)
    xx <- a*(X[i-1] - 127773*K) - 2836*K + 13
    X[i] <- ifelse(xx < 0, xx + M, xx)
  }
  
  if(!is.null(SEED))  SEED <<- X[i]
  
  return(X[6:n])
}


## generate uniform distribution
xunif <- function(n, min = 0, max = 1) {

  n = as.integer(n)
  if (is.na(n)) {
    stop("n has to be an integer.")
  }
  
  R <- .desertIsland(n)
  U <- R * M_INV
  return(min(min,max) + abs(max-min) * U) ## to accommodate user errors where b < a
}


## to generate n random integers in the range of (a,b]
## parameters:
##   n, number of random integers
##   a, lower bound
##   b, upper bound
## return:
##   random integers range from (a, b].

## the function can also deal with situation where b < a

xrandint <- function(n, min = 0, max = 10) {
  U <- xunif(n)
  min_ <- as.integer(min)
  max_ <- as.integer(max)
  
  if (any(is.na(c(min_, max_)))) {
    stop(sprintf('Parameters "min" and "max" has to be integers. Got min = %s; max = %s', min, max))
  }
  
  return(ceiling(min(min_,max_) + abs(max_ - min_)*U))
}

## bern distribution
xbern <- function(n, prob = 0.5) {
  if (prob < 0 || prob > 1) {
    stop(sprint('prob should be within [0, 1], got %s', prob))
  }
  
  U <- xunif(n)
  x <- ifelse(U <= prob, 1, 0)
  
  return(x)
  
}

## generate binomial distribution
xbinom <- function(n, size = 10, prob = 0.5) {
  if(!is.numeric(size) || size %% 1 != 0 || size < 1) {
    stop('size must be a positive integer.')
  }
  if(prob > 1 || prob < 0) stop('p must be in the range of [0,1]')
  
  B <- xbern(n * size, prob = prob)
  Bx <- matrix(B, nrow = size, ncol = n)
  
  return(colSums(Bx))
  
}

## generate poisson distribution
xpois <- function(n, lambda = 1){
  if(!is.numeric(lambda) || lambda <= 0) stop('lambda has to be a positive number.')
  
  xfunc <- function() {
      a <- exp(-lambda)
      p <- 1
      X <- -1
      while(a < p) {
        U <- xunif(1)
        p <- p*U
        X <- X+1
      }
      return(X)
  }
  
  if(lambda > 20) {
    R <- sapply(1:n, function(x) xfunc())
  } else {
    Z <- xnorm(n)
    R <- lambda + Z * sqrt(lambda) + 0.5
    R <- ifelse(R<0, 0, floor(R))
  }
  return(R)
}

## generate negative binomial distribution
## generates as sum of geom distribution
## modeled as # of total trials until r-th (size) success 
## prob is the probability of success
xnbinom <- function(n, size = 10, prob = 0.5) {
  
  if(!is.numeric(size) || size %% 1 != 0 || size < 1) {
    stop('r must be a positive integer.')
  }
  if(prob > 1 || prob < 0) stop('p must be in the range of [0,1]')
  
  G <- xgeom(n * size, prob = prob)
  Gx <- matrix(G, nrow = size, ncol = n)
  
  return(colSums(Gx))
  
}

## geometric distribution
## number of trials before success
## f(x|p) = p*(1-p)^(x-1), x ~ c(1,2,3,...), p ~ [0,1]
xgeom <- function(n, prob = 0.5) {
  
  if(prob > 1 || prob < 0) stop('p must be in the range of [0,1]')
  
  if(prob == 1) {
    return(rep(1, n))
  }
  U <- xunif(n)
  X <- ceiling(log(U)/log(1-prob))
  
  return(X)
  
}

## generate exponential distribution
## rate is the same as lambda
xexp <- function(n, rate = 1) {
  U <- xunif(n)
  return(-log(1-U)/rate)

}

## generate normal distribution
xnorm <- function(n, mean = 0, sd = 1, method = c('box-muller', 'polar')) {
  if(!is.numeric(sd) || as.numeric(sd) < 0) {
    stop('sd must be a positive number')
  }
  
   method = method[1]
   if (!(method %in% c('polar', 'box-muller'))) {
     message('Only support polar and box-muller methods; set to polar method')
     method = 'polar'
   }
   
   if(method == 'polar') {
    z <- xnormPL(n)
   } else {
    z <- xnormBM(n)
   }
   
   return(z*sd + mean)
}

## generate Standard Normal Distribution using Polar method
xnormPL <- function(n) {
  polar_func <- function() {
    while(TRUE) {
      U <- xunif(2)
      V <- 2*U-1
      W <- sum(V**2)
      if(W <= 1) {
        Y <- sqrt(-2*log(W)/W)
        return(Y * V)
        }
     }
  }
  
  Z <- replicate(n, polar_func())
  Z <- as.vector(Z)[xrandint(n, 1, 2*n)]
  return(Z)
  
}

## generate Standard Normal Distribution using Box-muller method
xnormBM <- function(n) {
  nx <- 2*n
  U <- xunif(nx)
  U1 <- U[(1:nx) %% 2 == 1]
  U2 <- U[(1:nx) %% 2 == 0]
  z1 <- sqrt(-2*log(U1)) * cos(2*pi*U2)
  z2 <- sqrt(-2*log(U1)) * sin(2*pi*U2)
  z <- c(z1, z2)[xrandint(n, 0, 2*n)]
  return(z)
  
}

## generate triangle distribution
## min = a, max = c, mode = b
xtri <- function(n, min = 0,mode = 1, max = 2) {
  if (any(!is.numeric(c(min, max, mode)))) stop('all parameters must be numeric')
  if (!(min < mode)) stop('min must less than mode') 
  if (!(mode < max)) stop('mode must less than max')  
  
  f <- (mode - min)/(max - min)
  U <- xunif(n)
  U <- ifelse(U<=f, 
              min + sqrt(U*(mode-min)*(max-min)), 
              max - sqrt((1-U)*(max-mode)*(max-min)))
  
  return(U)
}

## generate Erlang distribution
## summation of log(U) * (-1/lambda)
## become inefficient when n is big
## Erlang is a special case of gamma
## therefore can use gamma function to generate Erlang
xerlang <- function(n, shape = 1, rate = 1) {
  if(!is.numeric(shape) || shape %% 1 != 0 || shape < 1) {
    stop('shape must be a positive integer.')
  }
  if(!is.numeric(rate) || rate <= 0) stop('rate must be positive.')
  
  U <- xunif(n * shape)
  Ux <- matrix(log(U), nrow = shape, ncol = n)
  
  return(-colSums(Ux)/rate)
}


## generate gamma distribution
## shape = alpha, rate = beta,  or lambda 
## acceptance - rejection
## From Law2015 book page 453 - 455
xgamma <- function(n, shape = 1, rate = 1) {
  if(!is.numeric(shape) || shape <= 0) stop('shape has to be positive.')
  if(!is.numeric(rate) || rate <= 0) stop('rate  has to be positive.')
  
  GS_gamma <- function(shape, b) {
    while(TRUE) {
      U <- xunif(2)
      P <- b*U[1]
      
      if(P > 1) {
        ## step 3
        X <- -log((b-P)/shape)
        if(U[2] <= X**(shape-1)) return(X)
      } else {
        ## step 2
        X <- P**(1/shape)
        if(U[2] <= exp(-X)) return (X)
      }
    }
  }
  
  GB_gamma <- function(a, b, d, q, shape, theta) {
    while(TRUE){
      U <- xunif(2)
      rU <- U[1]/(1-U[1])
      V <- a*log(rU)
      Z <- U[2]*U[1]*U[1]
      
      Y <- shape*exp(V)
      W <- b + q*V-Y
      
      if((W+d-theta*Z) >= 0) {
        return(Y)
      } else {
        if(W >= log(Z)) return(Y)
      }
    }
  }
  
  if (shape < 1) {
    b <- (exp(1) + shape)/exp(1)
    G <- replicate(n, GS_gamma(shape, b))
    
  } else if (shape > 1) {
    a <- 1/sqrt(2*shape - 1)
    b <- shape - log(4)
    q <- shape + 1/a
    theta <- 4.5
    d <- 1 + log(theta)
    G <- replicate(n, GB_gamma(a, b, d, q, shape, theta))
    
  } else {
    ## shape == 1
    G <- xexp(n, rate = 1)
  }
  return(rate * G)
}

## generate beta distribution
## B(a,b) ~ X/(X+Y) where X ~ gamma(a, 1) and Y ~ gamma(b, 1)
## shape1 = alpha, shape2 = beta
xbeta <- function(n, shape1 = 1, shape2 = 1) {
  if(!is.numeric(shape1) || shape1 <= 0) stop('shape1 has to be positive.')
  if(!is.numeric(shape2) || shape2 <= 0) stop('shape2 has to be positive.')
  
  X <- xgamma(n, shape = shape1, rate = 1)
  Y <- xgamma(n, shape = shape2, rate = 1)
  return(X/(X+Y))
  
}

## generate weibull distribution
## PDF: f(x) = a*r*(a*x)^(r-1)*exp(-(a*x)^r), a = lambda
## CDF: F(x) = 1 - e**(-(lambda*x)**r)
## here shape = r; lambda = 1/scale of rweibull 
## inverse transformation
xweibull <- function(n, shape = 1, lambda = 1, scale = NULL) { 
  
  if(is.null(lambda) && !is.null(scale)) lambda = 1/scale
  if(!is.numeric(shape) || shape < 0) stop('shape cannot be negative')
  if(!is.numeric(lambda) || lambda < 0) stop('lambda cannot be negative.')
  
  U <- xunif(n)
  W <- (-log(1-U))**(1/shape)/lambda
  
  return(W)
}

## cauchy distribution
## inverse transformation.
xcauchy <- function(n, location = 0, scale = 1) {
  
  U <- xunif(n)
  C <- location + scale * tan(pi * (U-0.5))
  
  return(C)
}











