
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DRIP

<!-- badges: start -->

[![R-CMD-check](https://github.com/kangy10/DRIP/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kangy10/DRIP/actions/workflows/R-CMD-check.yaml)
[![Codecov](https://codecov.io/gh/kangy10/DRIP/branch/main/graph/badge.svg)](https://app.codecov.io/gh/kangy10/DRIP)
<!-- badges: end -->

DRIP is a collection of functions for edge detection, image denoising
and image deblurring. Below we categorize the functions based on their
purpose.

## Function overview

### Edge detection

| Function            | Description                             |
|---------------------|-----------------------------------------|
| `stepDiff()`        | Step edge detection statistics          |
| `stepEdge()`        | Step edge detection                     |
| `stepEdgeParSel()`  | Step edge detection parameter selection |
| `roofDiff()`        | Roof edge detection statistics          |
| `roofEdge()`        | Roof edge detection                     |
| `roofEdgeParSel()`  | Roof edge detection parameter selection |
| `edgeModify1()`     | Edge modification                       |
| `edgeModify2()`     | Edge modification                       |
| `edgeParSelPilot()` | Edge detection parameter selection      |
| `dKQ()`             | Edge detection performance measure      |

### Image denoising

| Function | Description |
|----|----|
| `restore3Stage()` | Three-stage denoising |
| `restore3StageParSel()` | Three-stage denoising parameter selection |
| `JPLLK_surface()` | Gradient-based denoising |
| `surfaceCluster()` | Clustering-based denoising |
| `surfaceCluster_bandwidth()` | Clustering-based denoising parameter selection |

### Image deblurring

| Function    | Description                         |
|-------------|-------------------------------------|
| `jpex()`    | JPEX deblurring                     |
| `cv.jpex()` | JPEX deblurring parameter selection |

## Installation

You can install the CRAN version of DRIP with:

``` r
install.packages("DRIP")
```
