## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(DRIP)

## ----step, fig.width=6--------------------------------------------------------
stepedge <- stepEdge(image = sar, bandwidth = 10, thresh = 17, 
                     degree = 1)
par(mfrow = c(1, 2), mar = c(3, 1, 1, 1), xaxt = "n", yaxt = "n")
image(sar, col = gray(c(0:255)/255))
image(1 - stepedge, col = gray(c(0:255)/255))

## ----step-par-pilot, fig.width=7----------------------------------------------
edgeParSelPilot(sar, edgeType = "step", degree = 1, bandwidth = c(6, 8, 10), probs = c(0.75, 0.85, 0.95))

## ----step-par-----------------------------------------------------------------
set.seed(24)
parSel <- stepEdgeParSel(image = sar, bandwidth = c(9, 10), degree = 1, 
                          thresh = c(17, 21), nboot = 10)
print(parSel, type = "all")

## ----three-stage--------------------------------------------------------------
fit <- restore3Stage(image = sar, bandwidth = 4, step_edge = stepedge, 
                  roof_edge = array(0, dim(sar)))
par(mfrow = c(1, 1), mar = c(3, 1, 1, 1), xaxt = "n", yaxt = "n")
image(fit, col = gray(c(0:255)/255))

## -----------------------------------------------------------------------------
bw_3stage <- restore3StageParSel(image = sar, bandwidth = 4:5, 
                              step_edge = stepedge, 
                              roof_edge = array(0, dim(sar)))
print(bw_3stage, type = "all")

## ----jpex, fig.width=6--------------------------------------------------------
deblur <- jpex(image = stopsign, bandwidth = 2, sigma = 0.00623, 
               alpha = 0.001)
names(deblur)
par(mfrow = c(1, 2), mar = c(3, 1, 1, 1), xaxt = "n", yaxt = "n")
image(stopsign, col = gray(c(0:255)/255))
image(deblur$deblurred, col = gray(c(0:255)/255))

## ----cv-jpex------------------------------------------------------------------
cv.out <- cv.jpex(image = stopsign, bandwidths = c(2, 3), ncpus = 1)
print(cv.out, type = "all")

