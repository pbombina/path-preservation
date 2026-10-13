### (C) Copyright 2026 Polina Bombina and Kevin R. Coombes
###
### Functions to compute path metrics
### Simplified into its own "source-able" file by KRC

###################################
### Smoothness
###
comp_smooth <- function( coords) {
  ticks <- 1:nrow(coords) # arbitrary parametrization
  ## Adaptive loess span: keep >= ~10 points per local window so the fit stays
  ## defined for short paths. At nrow(coords) >= 100 this equals the original
  ## span = 0.1, so longer paths are unaffected; below 100 it widens the window
  ## (e.g. a 15-point loop uses span = 10/15 ~= 0.67) instead of returning NaN.
  sp <- max(0.1, 10 / nrow(coords))
  ## fit smooth curve to each coordinate as a function of the parameter
  fitted <- apply(coords, 2, function(X) predict(loess(X ~ ticks, span = sp)))
  delta <- fitted - coords # difference between actual path and smooth fitted curve
  SE <- apply(delta^2, 1, sum) # Euclidean distance from path to fitted curve
  sqrt(SE)
}
Smoothness <- function(M1, M2) {
  ## RMSE aggregation (sqrt(mean(SE))), matching the metric's definition in
  ## the manuscript ("root mean squared deviation from the fitted curve") and
  ## the original exploratory script -- not a plain mean of per-point
  ## deviations. comp_smooth() already returns sqrt(SE) per point, so squaring
  ## it back before averaging recovers SE.
  log(sqrt(mean(comp_smooth(M1)^2)) / sqrt(mean(comp_smooth(M2)^2)))
}
RankSmoothness <- function(M1, M2) {
  cor(comp_smooth(M1), comp_smooth(M2), method = "spearman")
}


###################################
### Coil
###
comp_coil <- function(coords) {
  troid <- apply(coords, 2, mean)
  temp <- rbind(coords, troid)
  ed <- dist(temp)
  L <- nrow(temp)
  dd <- as.matrix(ed)[L,][-L]
  diff(range(dd))
}
Coil <- function(M1, M2) {
  log(comp_coil(M1) / comp_coil(M2))
}

###################################
### Contact
###
comp_contact <- function(coords, q) {
  ed <- as.matrix(dist(coords)) # small values here mean close contact
  q10 <- quantile(ed, q)        # use a quantile to decide what small means
  ai <- which(ed < q10, arr.ind = TRUE) # get matrix coordinates of close cells
  spread <- abs(apply(ai, 1, diff)) # number of steps apart in the sequence
  mean(spread) # not sure if this is the best statistical summary
}
Contact <- function(M1, M2, q = 0.1) {
  log(comp_contact(M1, q) / comp_contact(M2, q))
}


###################################
### EndpointDisplacement
###
endpt_disp <- function(coords, k) {
  start <- apply(coords[1:k,], 2, mean)
  L <- nrow(coords)
  end <- apply(coords[(L-k+1):L, ], 2, mean)
  sqrt(sum((start-end)^2))
}
EndpointDisplacement <- function(M1, M2, k=5) {
  log(endpt_disp(M1, k) / endpt_disp(M2, k))
}

###################################
### Pairwise Distance Rank
### (Spearman correlation of all pairwise distances, not just along the path;
### used alongside the cell-cycle/cyclic metrics below)
###
RankPairwiseDistance <- function(M1, M2) {
  cor(as.vector(dist(M1)), as.vector(dist(M2)), method = "spearman")
}

###################################
### Radial Spread
### (no analogue for an open path; distance of each point from the path's own
### centroid is only meaningful when the path is a closed loop, as in the
### cell-cycle data)
###
comp_radial <- function(coords) {
  troid <- apply(coords, 2, mean)
  sqrt(rowSums(sweep(coords, 2, troid, "-")^2))
}
RadialSpreadPreservation <- function(M1, M2) {
  cor(comp_radial(M1), comp_radial(M2), method = "spearman")
}

###################################
### Cyclic Length Distortion, Curvature, and Segment Variance
### Closed-loop analogues of LengthDistortion, Curvature, and SegmentVariance
### that wrap the last point in the path around to the first, as needed for
### a cell-cycle loop instead of an open path.
###
comp_cyclic_edges <- function(coords) {
  coords <- as.matrix(coords)
  n <- nrow(coords)
  nxt <- c(2:n, 1)
  sqrt(rowSums((coords[nxt, , drop = FALSE] - coords)^2))
}
CyclicLengthDistortion <- function(M1, M2, eps = 1e-10) {
  log((sum(comp_cyclic_edges(M2)) + eps) / (sum(comp_cyclic_edges(M1)) + eps))
}

comp_cyclic_curvature <- function(coords) {
  coords <- as.matrix(coords)
  n <- nrow(coords)
  prv <- c(n, 1:(n - 1))
  nxt <- c(2:n, 1)
  v1 <- coords - coords[prv, , drop = FALSE]
  v2 <- coords[nxt, , drop = FALSE] - coords
  len1 <- sqrt(rowSums(v1^2))
  len2 <- sqrt(rowSums(v2^2))
  cosang <- rowSums(v1 * v2) / (len1 * len2)
  cosang[len1 == 0 | len2 == 0] <- NA_real_
  acos(pmin(pmax(cosang, -1), 1))
}
CyclicCurvature <- function(M1, M2, eps = 1e-10) {
  log((mean(comp_cyclic_curvature(M2), na.rm = TRUE) + eps) /
      (mean(comp_cyclic_curvature(M1), na.rm = TRUE) + eps))
}

CyclicSegmentVariance <- function(M1, M2, eps = 1e-10) {
  log((var(comp_cyclic_edges(M2)) + eps) / (var(comp_cyclic_edges(M1)) + eps))
}

###################################
### Crossing Number
### (self-intersection count of a single 2-D embedded path; unlike the other
### metrics in this file, it has no high-dimensional counterpart to compare
### against, so it takes one matrix rather than a pair)
###
segments_intersect <- function(p1, p2, p3, p4, eps = 1e-12) {
  orient <- function(a, b, c)
    (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])
  o1 <- orient(p1, p2, p3); o2 <- orient(p1, p2, p4)
  o3 <- orient(p3, p4, p1); o4 <- orient(p3, p4, p2)
  (o1 * o2 < -eps) && (o3 * o4 < -eps)
}
comp_crossings <- function(coords) {
  coords <- as.matrix(coords)
  n <- nrow(coords)
  if (n < 4) return(0)
  count <- 0
  for (i in 1:(n - 3)) {
    for (j in (i + 2):(n - 1)) {
      if (j == i + 1) next
      if (segments_intersect(coords[i, ], coords[i + 1, ],
                              coords[j, ], coords[j + 1, ])) count <- count + 1
    }
  }
  count
}
CrossingNumber <- function(M) {
  comp_crossings(M)
}
CrossingNorm <- function(M, eps = 1e-8) {
  comp_crossings(M) / (nrow(as.matrix(M))^2 + eps)
}
