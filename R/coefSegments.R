#' @title Coefficients from segments
#' @describe Computes the trend and amplitudes of segments derived from BFAST-type models
#' 
#' @param x Numeric vector
#' @param dates Vector of dates corresponding to \code{x}. Can be in format \code{"\%Y-\%m-\%d"} or as decimal years.
#' @param type Character. Time series type (see \code{\link{bfastts}})
#' @param bp Vector of dates. Can be in format \code{"\%Y-\%m-\%d"} or as decimal years. Omit to compute breakpoints automatically (see \code{\link{breakpoints}}). Set to \code{0} for a single segment model (no breakpoints).
#' @param model Character. Either a OLS model fit (\code{'lm'}) or a robust linear model (\code{'rlm'}; see package \code{\link{MASS}})
#' @param formula Either \code{response ~ harmon} or \code{response ~ harmon + trend}
#' @param order Numeric. Order of the harmonic component (see \code{\link{bfastpp}}). Only works with \code{order=1} for now.
#' @param plot Logical. Produce a plot?
#' @param dataLabel Character. Label of the y-axis of the data plot, if applicable.
#' @param ... Additional arguments to be passed to \code{\link{breakpoints}} (if \code{breakpoints} argument is omitted).
#'
#' @return \code{matrix} with 2 columns (for amplitude and trend, respectively) and 1 row for each chronological segment detected in \code{x}
#' 
#' @details If \code{bp} is left as \code{NULL}, breakpoints will be automatically detected using \code{\link{breakpoints}}. However, the model supplied in \code{model} is not used in this step - these are only fitted to each segment after breakpoints are determined.
#'
#' @author Ben DeVries
#'   
#' @import bfast
#' @import zoo
#' @import strucchange
#' @import MASS
#' @export
#' 

coefSegments <- function(x, dates, type = 'irregular', bp = NULL, model = c('lm', 'rlm'), formula = response ~ harmon + trend, order = 1, plot = FALSE, dataLabel = 'data', ...) {
  
  
  #### Pre-processing
  
  # strip x to a vector (e.g. in case a zoo is supplied)
  x <- as.numeric(x)
  
  # reformat dates if needed (should be in yyyy-mm-dd for bfastts)
  if(is.numeric(dates)) {
    year <- floor(dates)
    jday <- round((dates - year) * 365, 0)
    dates <- as.Date(paste(year, jday, sep = ''), format = '%Y%j')
  }
  
  # make regular bfast-type time series
  bts <- bfastts(x, dates, type = type)
  
  # bfast pre-processing
  bpp <- bfastpp(bts, order=order)
  
  
  #### Breakpoints
  
  # compute breakpoints if necessary and assign segments
  if(is.null(bp)) {
    
    bp <- breakpoints(response ~ time, data = bpp, ...)
    bpp$segment <- breakfactor(bp)
    bd <- bpp$time[bp$breakpoints]
    
  } else if(bp[1] != 0) {
    
    # reformat bp to numeric if needed
    if(!is.null(bp) & class(bp) != 'numeric') {
      bp <- as.numeric(format(bp, format = '%Y')) + as.numeric(format(bp, format = "%j"))/365
      bd <- bp
    } else if(class(bp) == 'numeric' & bp[1] != 0) {
      bd <- bp
    }
    
    seg <- paste('segment', c(1:(length(bp) + 1)), sep = '')
    bpp$segment <- seg[1]
    for(i in 2:length(seg)) {
      bpp$segment[bpp$time > bp[i-1]] <- seg[i]
    }
    
  } else if(bp[1] == 0) {
    
    bpp$segment <- 'segment1'
    bp <- NULL
    bd <- NULL
    
  }
  
  
  #### Models
  
  # OLS or robust lm
  model <- model[1]
  
  if(length(unique(bpp$segment)) > 1) {
    
    # check formula (only response ~ harmon or response ~ harmon + trend allowed)
    if(all(as.character(formula) == c("~", "response", "harmon"))) {
      formula <- response ~ segment/harmon
    } else if(all(as.character(formula) == c("~", "response", "harmon + trend")) | all(as.character(formula) == c("~", "response", "trend + harmon"))) {
      formula <- response ~ segment/(harmon+trend)
    } else {
      stop('only harmon and harmon+trend models supported at this time.')
    }
    
    # check fit (only lm or rlm allowed)
    if(model == 'lm') {
      m <- lm(formula, data=bpp)
    } else if(model == 'rlm') {
      m <- rlm(formula, data=bpp)
    } else {
      stop("model must be either \'lm\' or \'rlm\'.")
    }
    
  } else {
    
    if(model == 'lm') {
      m <- lm(formula, data=bpp)
    } else if(model == 'rlm') {
      m <- rlm(formula, data=bpp)
    } else {
      stop("model must be either \'lm\' or \'rlm\'.")
    }
    
  }
  
  
  #### Amplitudes
  
  # compute amplitude
  coef <- coefficients(m)
  cosines <- coef[which(grepl('harmoncos', names(coef)))]
  sines <- coef[which(grepl('harmonsin', names(coef)))]
  amps <- sqrt(cosines^2 + sines^2)
  names(amps) <- unique(bpp$segment)
  
  
  #### Trends
  
  if(any(grepl('trend', as.character(formula)))) {
    trends <- coef[which(grepl('trend', names(coef)))]
  } else {
    trends <- NULL
  }
  
  
  
  #### Lengths
  if(length(unique(bpp$segment) > 1)) {
    len <- c(max(bpp$time[bpp$segment == 'segment1']) - min(bpp$time),
             max(bpp$time) - min(bpp$time[bpp$segment == 'segment2']))
  } else {
    len <- max(bpp$time) - min(bpp$time)
  }
  
  #### plots
  
  # plot
  if(plot) {
    
    # predict values
    bpp$prediction <- predict(m, newdata = bpp)
    
    y <- zoo(x, dates)
    time(y) <- as.numeric(format(time(y), "%Y")) + as.numeric(format(time(y), "%j"))/365
    
    lo <- matrix(c(1:2), nr=2, nc=1)
    layout(lo)
    op <- par(mar = c(0, 5, 0, 5), oma = c(3, 0, 2, 0))
    
    # plot limits
    xmax <- ceiling(max(time(y)))
    xmin <- ceiling(min(time(y)))
    ymax <- max(bpp$response)
    ymin <- min(bpp$response)
    ymax <- ymax + ((ymax - ymin) / 10)
    ymin <- ymin - ((ymax - ymin) / 10)
    
    plot(y, ylab = dataLabel, xlab = "", type = 'b', pch = '*', xlim = c(xmin, xmax), xaxt = 'n', ylim = c(ymin, ymax))
    
    ## trend plot
    if(any(grepl('trend', as.character(formula)))) { 
      for(i in 1:length(unique(bpp$segment))) {
        seg <- unique(bpp$segment)[i]
        
        trend <- subset(bpp, segment == seg)
        trend$harmon[] <- 0
        trend$prediction <- predict(m, newdata = trend)
        
        # plot trend
        lines(zoo(trend$prediction[trend$segment == seg], trend$time[trend$segment == seg]), col = 4)
      }
    }
    
    abline(v = bd, lty = 3)
    
    ## season plot
    ymax <- max(coef[1] + amps)
    ymin <- min(coef[1] - amps)
    ymax <- ymax + ((ymax - ymin) / 10)
    ymin <- ymin - ((ymax - ymin) / 10)
    plot(y, type = 'n', xlim = c(xmin, xmax), ylab = 'season', ylim = c(ymin, ymax), yaxt = 'n', xlab = 'Time', xaxt = 'n')
    axis(4)
    
    # plot intercept
    lines(zoo(coef[1], bpp$time), col = 4, lty = 2)
    
    for(i in 1:length(unique(bpp$segment))) {
      seg <- unique(bpp$segment)[i]
      
      seas <- subset(bpp, segment == seg)
      seas$trend <- 0
      seas$prediction <- predict(m, newdata = seas)
      
      if(i > 1) {
        seas$prediction[seas$segment == seg] <- seas$prediction[seas$segment == seg] - coefficients(m)[i]
      }
      
      # plot season model
      lines(zoo(seas$prediction[seas$segment == seg], seas$time[seas$segment == seg]), col = 4, xaxt = 'n', xlim = c(xmin, xmax))
      
      # plot amplitude
      lines(zoo(coef[1] + amps[i], seas$time[seas$segment == seg]), lty = 4, col = 4)
      lines(zoo(coef[1] - amps[i], seas$time[seas$segment == seg]), lty = 4, col = 4)
    }
    
    abline(v = bd, lty = 3)
    
    axis(1, at = c(xmin:xmax), labels = c(xmin:xmax))
    par(op)
    layout(1)
  }
  
  
  if(!is.null(trends)) {
    coefs <- cbind(amps, trends, len)
    colnames(coefs) <- c('amplitude', 'trend', 'length')
  } else {
    coefs <- cbind(amps, len)
    colnames(coefs) <- c('amplitude', 'length')
  }
  
  return(coefs)
}

