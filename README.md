# Integrated Forest Monitoring using LTS and CBM

Data and R code accompanying the paper:

DeVries, B., Pratihast, A.K., Verbesselt, J., Kooistra, L. and Herold, M. 2016. Characterizing forest change using community-based monitoring data and Landsat time series. *PLOS ONE*, in press.

## Introduction

This paper presents an approach to integrate community-based monitoring (CBM) data with dense Landsat times series (LTS) to characterize subtle forest changes. Here, we focus on deforestation and degradation processes, but this approach can be used to measure and monitor a broader range of forest change variables.

This repository holds the local data and the R code used to carry out this study. Because of size restrictions, the LTS data are not included here (except for a time series for an example location). Landsat image data can be downloaded for free via the USGS.

## Community-Based Monitoring data

In this study, we worked with local forest rangers to demonstrate the potential of CBM data when integrated with LTS. A series of observations made by the rangers over a period of three years is included here. These observations were labelled as either deforestation ("DEF") or degradation ("DEG") events. We also included a set of no-change ("NOCH") observations based on visual interpretation of a random sample using high-resolution imagery to complete the training dataset.

```R
library(rgdal)
library(raster)
obs <- readOGR('data', 'local_observations')

plot(obs, col = 'white')
plot(subset(obs, label == "NOCH"), col = '#4daf4a', pch = '*', add = TRUE)
plot(subset(obs, label == "DEG"), col = '#377eb8', pch = '+', add = TRUE)
plot(subset(obs, label == "DEF"), col = '#e41a1c', pch = 'X', cex = 0.7, add = TRUE)
scalebar(25000, label = '25km')
legend('topright', legend = c('DEF', 'DEG', 'NOCH'), pch = c('X', '+', '*'), 
       pt.cex = c(0.7, 1, 1), col = c('#e41a1c', '#377eb8', '#4daf4a'), bty = 'n')
```

<div style="text-align:center">
<img src ="figs/local_observations.png" />
</div>


## Spectral bands and indices

In this paper we used a range of indices in addition to the original six Landsat bands (excluding the thermal band). The original six bands are included here for an example observation (where deforestation ("DEF") was observed).

```R
LTS <- read.csv('data/LTS.csv')
head(LTS)
```

Each column represents a spectral band, and can be neatly visualized using the ```zoo``` package:

```R
library(zoo)
swir1 <- zoo(LTS$SWIR1, as.Date(LTS$date))
plot(swir1, type = 'b', pch = '*', cex = 0.65)
```

We computed a range of indices based on these spectral bands:

```R
LTS$NDVI <- (LTS$NIR - LTS$R) / (LTS$NIR + LTS$R)
LTS$NDMI <- (LTS$NIR - LTS$SWIR1) / (LTS$NIR + LTS$SWIR1)
LTS$NBR <- (LTS$NIR - LTS$SWIR2) / (LTS$NIR + LTS$SWIR2)
LTS$NBR2 <- (LTS$SWIR1 - LTS$SWIR2) / (LTS$SWIR1 + LTS$SWIR2)
```

For tasseled cap indices, we used a single set of coefficients, assuming that the atmospheric correction also corrected for sensor-dependent differences.

```R
brightness <- c(0.2043, 0.4158, 0.5524, 0.5741, 0.3124, 0.2303)
greenness <- c(-0.1603, -0.2819, -0.4934, 0.7940, -0.0002, -0.1446)
wetness <- c(0.0315, 0.2021, 0.3102, 0.1594, -0.6806, -0.6109)

LTS$TCB <- apply(LTS[, c(2:7)], 1, FUN=function(x) sum(x * brightness))
LTS$TCG <- apply(LTS[, c(2:7)], 1, FUN=function(x) sum(x * greenness))
LTS$TCW <- apply(LTS[, c(2:7)], 1, FUN=function(x) sum(x * wetness))
```

And finally, the tasseled cap angle (TCA) was computed as the angle between the TCB and TCG vectors.

```R
LTS$TCA <- atan(LTS$TCG / LTS$TCB)
```

## Random Forest covariates

We used several data reduction approaches to derive random forest covariates for *each* of spectral bands and indices (hereafter simply referred to as "bands"). The covariates used for the random forest models in this paper were based on two types of LTS reductions: (1) full time series descriptors; and (2) segment-based descriptors.

### (1) Full time series descriptors

We used the robust linear model (rlm) from the ```MASS``` R package to fit a linear model over the entire time series for each of the spectral bands. We first converted the irregularly sampled LTS to a regular "daily" time series, filling in missing observations with ```NA``` using the ```bfast``` package, as in this example using the ```TCW``` band.

```R
library(bfast)
bts <- bfastts(LTS$TCW, dates = as.Date(LTS$date), type = "irregular")
bpp <- bfastpp(bts, order = 1)
```

```bfastts``` produces a ```ts``` object, with a time index expressed as a numeric (ie. decimal year) and ```bfastpp``` produces a ```data.frame``` from this time series that we can use to easily produce fitted models.

```R
library(MASS)
m <- rlm(response ~ as.numeric(time), data = bpp)
rlmIntercept <- coef(m)[1] + coef(m)[2] * min(bpp$time)
rlmTrend <- coef(m)[2]
bpp$pred_rlm <- predict(m, newdata = bpp)

plot(zoo(bpp$response, bpp$time), type = 'b', pch = '*', cex = 0.65, ylab = "TCW", xlab = "Time")
lines(zoo(bpp$pred_rlm, bpp$time), col = 'red', lty = 2)
```

<div style="text-align:center">
<img src ="figs/rlm_int_trend.png" />
</div>

### (2) Segment-based descriptors

The full time series descriptor may be a good measure of gradual changes over time, but as we can see in the plot above, it does not fully represent abrupt changes. For this reason, we used the BFAST method to derive temporal segments, and fit RLM models to each of these segements individuals. For simplicity (ie. to keep the total number of covariates manageble), we limited the segments to a maximum of one segment.

The function ```coefSegments()``` included in this repository can derive and plot these segments for a time sereis. This function is based on the ```breakpoints()``` function in the ```strucchange``` package. To limit to a binary output (ie. maximum 1 breakpoint), set ```breaks = 1```.

Here, we use the TCW time series as an example, but this was performed for all spectral bands.

```R
library(strucchange)
source('R/coefSegments.R')
segs <- coefSegments(LTS$TCW, dates = as.Date(LTS$date), model = 'rlm', breaks = 1, plot = TRUE, dataLabel = 'TCW')
print(segs)
```

<div style="text-align:center">
<img src ="figs/segs.png" />
</div>

