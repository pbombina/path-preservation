# Preservation

Metrics for assessing how well dimension reduction (DR) algorithms preserve
the structure of high-dimensional data: pairwise distances, local
neighborhoods, and paths (open trajectories and closed loops such as a
cell-cycle trajectory).

Authors: Kevin R. Coombes (maintainer) and Polina Bombina.
The package is also distributed through the
[OOMPA project on R-Forge](http://oompa.r-forge.r-project.org/).

## Installation

```r
# install.packages("remotes")
remotes::install_github("pbombina/path-preservation")
```

## Metrics

| Family | Functions |
|---|---|
| Pairwise distances | `MilnorDistortion`, `SigmaDistortion`, `Stress`, `M1Distortion`, `SpearmanRho`, `EarthMover` |
| Neighborhoods | `AverageJaccardDistance`, `QNX`, `AvgContinuity`, `AvgTrustworthiness` |
| Paths (log-ratio) | `LengthDistortion`, `SegmentVariance`, `Curvature`, `SpatialSimilarity`, `Smoothness`, `Coil`, `Contact`, `EndpointDisplacement` |
| Paths (rank-based) | `RankLengthDistortion`, `RankCurvature`, `RankSpatialSimilarity`, `RankSmoothness`, `RankPairwiseDistance` |
| Closed loops | `CyclicLengthDistortion`, `CyclicCurvature`, `CyclicSegmentVariance`, `RadialSpreadPreservation` |
| Single embedded path | `CrossingNumber`, `CrossingNorm` |

## Example

```r
library(Preservation)
data("prot5")      # 3096 cells x 51 proteins (CAVA CITE-Seq subset)
data("drs")        # 2-D embeddings: MDS, TSNE, UMAP, DM
data("path-idex")  # idex: ordered row indices defining a path

hd <- prot5[idex, ]
ld <- drs[["UMAP"]][idex, ]

LengthDistortion(hd, ld)
RankCurvature(hd, ld)
CrossingNorm(ld)
```

See the vignette for a full walk-through: `vignette("CAVA", package = "Preservation")`.

## Reference

For the distance and neighborhood metrics, see
[doi:10.1101/2025.04.28.651069](https://doi.org/10.1101/2025.04.28.651069).

## License

Apache License 2.0
