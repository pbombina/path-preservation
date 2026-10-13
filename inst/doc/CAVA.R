## ----opts, echo=FALSE---------------------------------------------------------
knitr::opts_chunk$set(fig.width=7, fig.height=5)

## ----packs--------------------------------------------------------------------
library(Preservation)

## ----sampleData---------------------------------------------------------------
data("prot5")
dim(prot5)
table(colvec)

## ----drs----------------------------------------------------------------------
data("drs")
sapply(drs, dim)

## ----fig01, fig.width=8, fig.height=8, fig.cap="Four views of the CAVA subset after dimension reduction."----
opar <- par(mfrow = c(2,2))
for (tag in names(drs)) {
  plot(drs[[tag]], col = colvec, pch = 16, main = tag, xlab = "X1", ylab = "X2")
}
par(opar)

## ----PWmetrics----------------------------------------------------------------
hiD <- dist(prot5)
loD <- lapply(drs, dist)
suppressWarnings(
  scores <- sapply(c("MilnorDistortion", "SigmaDistortion", "Stress", 
                     "EarthMover"), function(metric) {
    cat(metric, "\n", file = stderr())
    FUN <- get(metric)
    sapply(names(drs), function(tag) {
      FUN(hiD, loD[[tag]])
    })
  })
)
m1d <- sapply(names(drs), function(tag) {
  M1Distortion(prot5, drs[[tag]])
})
scores <- cbind(scores, M1Distortion = m1d)
round(scores, 3)
Preservation:::PWBest(colnames(scores))

## ----neigh--------------------------------------------------------------------
ajd <-sapply(names(drs), function(tag) {
  AverageJaccardDistance(prot5, drs[[tag]], k = 50)
})
scores <- cbind(scores, ajd)
# too slow: qnx <- QNX(hiD, loD[["MDS"]])
quad <- c(MDS = "#B00068", TSNE = "#1C8356", UMAP = "#325A9B", DM = "#FEAF16" )
ct <- lapply(names(drs), function(tag) {
  Preservation:::compute_ct_list(prot5, drs[[tag]])
})
avtrust <- sapply(ct, function(X) X[,3])
colnames(avtrust) <- names(drs)

avcont <- sapply(ct, function(X) X[,6])
colnames(avcont) <- names(drs)

## ----fig.width = 9, fig.height=5----------------------------------------------
opar <- par(mfrow = c(1,2))
plot(avtrust[,1], type = "n", ylim = c(0.1, 0.3), xlab = "Neighbors", ylab = "Trustworthiness")
for (tag in names(drs)) lines(avtrust[, tag], lwd = 2, col = quad[tag])
legend("topright", names(drs), col = quad, lwd = 2)

plot(avcont[,1], type = "n", ylim = c(0, 0.25), xlab = "Neighbors", ylab = "Continuity")
for (tag in names(drs)) lines(avcont[, tag], lwd = 2, col = quad[tag])
legend("topleft", names(drs), col = quad, lwd = 2)
par(opar)

## ----idex---------------------------------------------------------------------
data("path-idex")
length(idex)
highD <- prot5[idex,]
dim(highD)

## ----fig02, fig.width=8, fig.height=8, fig.cap="Four views of the path of the first principal components of the CAVA subset after dimension reduction."----
opar <- par(mfrow = c(2,2))
for (tag in names(drs)) {
  lowD <- drs[[tag]][idex,]
  plot(drs[[tag]], col = colvec, pch = 16, main = tag, xlab = "X1", ylab = "X2")
  lines(lowD, lwd = 2)
}
par(opar)

## ----LD-----------------------------------------------------------------------
scores <- sapply(c("LengthDistortion", "SegmentVariance", "Curvature", "SpatialSimilarity"),
                 function(metric) {
  FUN <- get(metric)
  sapply(names(drs), function(tag) {
    lowD <- drs[[tag]][idex,]
    FUN(highD, lowD)
  })
})
round(scores, 3)

Preservation:::PathBest(colnames(scores))

## ----rk-----------------------------------------------------------------------
rscores <- sapply(c("RankLengthDistortion", "RankCurvature", "RankSpatialSimilarity"),
                 function(metric) {
  FUN <- get(metric)
  sapply(names(drs), function(tag) {
    lowD <- drs[[tag]][idex,]
    FUN(highD, lowD)
  })
})
round(rscores, 3)
scores <- cbind(scores, rscores)

## ----fig03, fig.width=8, fig.height=8, fig.cap="Four views of the smoothed PC1 path of the CAVA subset after dimension reduction."----
smoosh <- function( coords) {
  ticks <- 1:nrow(coords) # arbitrary parametrization
  ## fit smooth curve to each coordinate as a function of the parameter
  fitted <- apply(coords, 2, function(X) predict(loess(X ~ ticks, span = 0.1)))
  delta <- fitted - coords # difference between actual path and smooth fitted curve
  val <- list(coords = coords, fitted = fitted, delta = delta)
  class(val) <- "smoosh"
  val
}
fits <- lapply(names(drs), function(tag) smoosh(drs[[tag]][idex,]))
names(fits) <- names(drs)
plot.smoosh <- function(obj, main = "", col = "purple", ...) {
  plot(obj$coords, main = main, type = "b", col = "grey80", ...)
  lines(obj$fitted, col = col, lwd = 2)
  invisible(main)
}
lines.smoosh <- function(obj, col = "purple") {
  lines(obj$fitted, col = col, lwd = 2)
}
opar <- par(mfrow = c(2, 2))
for (tag in names(fits)) {
  plot(fits[[tag]], main=tag, col = "black", xlab = "X1", ylab = "X2")
  points(drs[[tag]], col = colvec, pch = ".", cex=3)
  lines(fits[[tag]], col = "black")
}
par(opar)

## ----morePathhMetrics---------------------------------------------------------
mscores <- sapply(c("Smoothness", "Coil", "Contact", "EndpointDisplacement"),
                 function(metric) {
  cat(metric, "\n", file = stderr())
  FUN <- get(metric)
  sapply(names(drs), function(tag) {
    lowD <- drs[[tag]][idex,]
    FUN(highD, lowD)
  })
})
round(mscores, 3)
scores <- cbind(scores, mscores)

## ----twist--------------------------------------------------------------------
round(t(scores), 3)

## ----si-----------------------------------------------------------------------
sessionInfo()

