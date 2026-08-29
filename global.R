suppressMessages({
  library(shiny)
  library(dplyr)
  library(stringr)
  library(rlang)
  library(randtests)
  library(EnvStats)
  library(rclipboard)
  library(base64enc)
})

home <- "./DistriFit"
home <- "./"
source(file.path(home, 'R/utilities.r'))
source(file.path(home, 'R/prn_generator.r'))
source(file.path(home, 'R/prob_fit.r'))
source(file.path(home, 'R/gofTest.chisq.r'))

## Global variables
M <- 2**31 - 1
M_INV <- 1/M
SEED <- NULL

MIN_OBS <- 30

DISTR.NAMES <- read.delim(file.path(home, 'R/distribution.names.txt'), header = T, sep = '\t') %>% 
  arrange(distr.type, distr.abbr)

DEFAULT_GEN_PARAMS <- setNames(
  lapply(seq_len(nrow(DISTR.NAMES)), function(i) {
    row <- DISTR.NAMES[i, ]
    np  <- row$n.params
    pnames <- unlist(row[paste0('param.', seq_len(np))])
    pvals  <- as.numeric(unlist(row[paste0('default.', seq_len(np))]))
    setNames(as.list(pvals), pnames)
  }),
  DISTR.NAMES$distr.abbr
)

DISTR_CHOICES <- setNames(
  DISTR.NAMES$distr.abbr,
  sprintf('%s (%s)', DISTR.NAMES$distr.name, DISTR.NAMES$distr.abbr)
)

getTopFits <- function(x, pcutoff, alpha, selected_distrs) {
    distrs <- if (is.null(selected_distrs) || "All applicable" %in% selected_distrs) NULL else selected_distrs
    n.top <- ifelse(!is.null(distrs), length(distrs), 3)
    fits <- tryCatch(
      distrEstimate(x = x, distribution = distrs, alpha = alpha, return.best = FALSE),
      error = function(e) NULL
    )
    if (is.null(fits) || length(fits) == 0) return(NULL)

    pvals      <- sapply(fits, function(r) r$p.value)
    res.sorted <- fits[order(pvals, decreasing = TRUE)]
    topFit       <- res.sorted[sapply(res.sorted, function(r) r$p.value >= pcutoff)]

    if (length(topFit) > 0) {
      topFit[seq_len(min(n.top, length(topFit)))]
    } else {
      res.sorted[seq_len(min(n.top, length(res.sorted)))]
    }
}


## the main function to plot data distribution and fitted density function
plotDistr <- function(data, top.results = NULL, is.fitted = FALSE, title = "") {
  n <- length(data)
  brks <- getBreaks(data)
  xr <- range(data)
  x_pad <- diff(xr) * 0.08
  if (x_pad == 0) x_pad <- 1
  xlims <- c(xr[1] - x_pad, xr[2] + x_pad)
  
  layout(matrix(c(1, 2), nrow = 2, ncol = 1), heights = c(1, 2.5))
  
  # --- SUBPLOT 1: Horizontal Boxplot ---
  par(mar = c(0, 4, 2, 4))
  
  mn  <- mean(data)
  med <- median(data)
  std   <- sd(data)
  min_val <- min(data)
  max_val <- max(data)
  
  mean_sd_col <- "salmon"
  med_col <- "gray30"
  min_col <- "orange"
  max_col <- "orange"
  hist_col0 <- 'steelblue2'
  pt_col <- col2rgb(hist_col0) %>% 
    as.list() %>% 
    c(., alpha = 0.5*255, maxColorValue = 255) %>% 
    do.call('rgb', .)
  
  plot(1, type = "n", xlim = xlims, ylim = c(0.2, 1.8),
       axes = FALSE, xlab = "", ylab = "", main = title, cex.main = 1.0)
  
  boxp <- boxplot(data, 
                  horizontal = T, 
                  axes = F,
                  xlim = c(0.2, 1.8), ylim = xlims,
                  xlab = '', ylab = '', main = title, cex.main = 1.0,
                  col = NA, border = rgb(0.4,0.4,0.4, alpha = 0.5), 
                  outline = F,
                  add = T)
  set.seed(64)
  jitter_y <- jitter(rep(1, length(data)-2), amount = 0.2)
  points(sort(data)[2:(n-1)], 
         jitter_y, 
         col = rep(pt_col, n-2), 
         pch = 16, 
         cex = 0.8)
  points(c(min_val, max_val),
         c(1,1),
         col = min_col,
         pch = 16, 
         cex=1)

  
  segments(mn, 0.8, mn, 1.2, col = mean_sd_col, lty = 1, lwd = 1.2)
  arrows(x0 = max(xlims[1], mn - std), y0 = 1, x1 = min(xlims[2], mn + std), y1 = 1, 
         code = 3, angle = 90, length = 0.06, col = mean_sd_col, lwd = 1.2)
  
  text(min_val, 1.45, labels = sprintf("Min: %.2f", min_val), col = min_col, cex = 0.75, font = 2)
  text(max_val, 1.45, labels = sprintf("Max: %.2f", max_val), col = max_col, cex = 0.75, font = 2)
  text(med, 1.55, labels = sprintf("Med: %.2f", med), col = med_col, cex = 0.75, font = 2)
  text(mn, 1.3, labels = sprintf("Mean \u00B1 sd: %.2f \u00B1 %.2f", mn, std), col = mean_sd_col, cex = 0.75, font = 2)
  
  # --- SUBPLOT 2: Histogram ---
  par(mar = c(4, 4, 0, 4))
  
  h <- hist(data, breaks = brks, plot = FALSE)
  hist_col <- if (is.fitted) "gray80" else hist_col0
  
  plot(h, freq = TRUE, col = hist_col, border = "white",
       xlim = xlims, main = "", xlab = "x", ylab = "Frequency", cex.axis = 0.85)
  
  if (is.fitted && !is.null(top.results) && length(top.results) > 0) {
    pdf_cols <- RColorBrewer::brewer.pal(3, 'Set2') 
    #xs <- seq(xlims[1], xlims[2], length.out = 500)
    xs <- seq(xr[1], xr[2], length.out = 500)
    max_dens <- 0
    curves_data <- list()
    
    for (i in seq_along(top.results)) {
      res   <- top.results[[i]]
      dabbr <- rev(res$distr.name)[1]
      params <- as.list(res$parameters)
      dens   <- tryCatch(getDensity(dabbr, xs, params), error = function(e) NULL)
      print(head(dens))
      
      if (!is.null(dens) && all(is.finite(dens))) {
        max_dens <- max(max_dens, max(dens, na.rm = TRUE))
        curves_data[[i]] <- list(dabbr = dabbr, dens = dens, res = res)
      }
    }
    
    if (length(curves_data) > 0) {
      max_dens <- max(max_dens * 1.15, 1e-4)
      par(new = TRUE)
      plot(xs, 
           rep(0, length(xs)), 
           type = "n", 
           xlim = xlims, 
           ylim = c(0, max_dens),
           axes = FALSE, xlab = "", ylab = "")
      axis(4, col = "#D62728", col.axis = "#D62728", cex.axis = 0.85)
      mtext("Density / PMF", side = 4, line = 2.5, col = "#D62728", cex = 0.8)
      
      legend.labels <- character(0)
      legend.cols   <- character(0)
      
      for (i in seq_along(curves_data)) {
        cd <- curves_data[[i]]
        if (is.null(cd)) next
        lines(xs, cd$dens, col = pdf_cols[i], lwd = 2.5)
        
        pval_str <- if (cd$res$p.value < 0.0001) sprintf("%.3e", cd$res$p.value) else sprintf("%.4f", cd$res$p.value)
        legend.labels <- c(legend.labels, sprintf("%s (p = %s)", cd$dabbr, pval_str))
        legend.cols   <- c(legend.cols, pdf_cols[i])
      }
      
      if (length(legend.labels) > 0) {
        legend("topright", legend = legend.labels, col = legend.cols, lwd = 2, bty = "n", cex = 0.8)
      }
    }
  }
}
