\name{More Path Preservation Metrics}
\alias{Smoothness}
\alias{RankSmoothness}
\alias{Coil}
\alias{Contact}
\alias{EndpointDisplacement}
\alias{RankPairwiseDistance}
\alias{RadialSpreadPreservation}
\alias{CyclicLengthDistortion}
\alias{CyclicCurvature}
\alias{CyclicSegmentVariance}
\alias{CrossingNumber}
\alias{CrossingNorm}

\title{ More Metrics for Path Preservation}
\description{
  Metrics to estimate how well or poorly a dimension reduction algorithm
  preserves paths, including closed-loop (cyclic) paths such as a
  cell-cycle trajectory.
}
\usage{
Smoothness(M1, M2)
RankSmoothness(M1, M2)
Coil(M1, M2)
Contact(M1, M2, q = 0.1)
EndpointDisplacement(M1, M2, k = 5)
RankPairwiseDistance(M1, M2)
RadialSpreadPreservation(M1, M2)
CyclicLengthDistortion(M1, M2, eps = 1e-10)
CyclicCurvature(M1, M2, eps = 1e-10)
CyclicSegmentVariance(M1, M2, eps = 1e-10)
CrossingNumber(M)
CrossingNorm(M, eps = 1e-8)
}
\arguments{
  \item{M1}{ a data matrix whose rows are the points in the original space.}
  \item{M2}{ a data matrix whose rows are the same points after
    dimension reduction.}
  \item{M}{ a data matrix whose rows are the points of a single embedded
    path (used by the metrics that have no high-dimensional counterpart).}
  \item{q}{a real number; the quantile of pairwise distances used to define "close".}
  \item{k}{an integer; the number of points to average at each end of the path.}
  \item{eps}{a small real number added to numerator and denominator of a
    log-ratio metric to keep it defined when a cyclic path has a
    degenerate (zero-length or zero-variance) segment, or to a denominator
    in \code{CrossingNorm} to avoid division by zero.}
}
\value{
  Each of the functions returns a real number representing the
  corresponding scores comparing the spaces before and after dimension reduction.
}

\author{Kevin R. Coombes <krc@silicovore.com>}
\examples{
data(prot5) # original space
data(drs) # spaces after dimension reduction
Smoothness(prot5, drs[["TSNE"]])
}
\keyword{htest}
\keyword{models}
