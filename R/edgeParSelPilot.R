# This is R source code for function 'edgeParSelPilot' in the R
# package 'DRIP'.
# Creator: Yicheng Kang
# Date: December 23, 2025

edgeParSelPilot <- function(image, edgeType = "step", degree = 0, blur = FALSE, bandwidth = c(3, 5, 7), probs = c(0.8, 0.9, 0.95)) {
  stopifnot(is.character(edgeType) && edgeType %in% c("step", "roof") && length(edgeType) == 1)
  nband <- length(bandwidth)
  nprobs <- length(probs)
  matout <- matrix(NA, nrow = nband, ncol = nprobs)
  row.names(matout) <- paste0("bandwidth=", bandwidth)
  colnames(matout) <- paste0("probs=", probs)
  par(mfrow = c(1, 3), mar = c(3, 1, 1, 1), xaxt = "n", yaxt = "n")
  if(edgeType == "step") {
    for (k in seq_along(bandwidth)) {
      diffout <- stepDiff(image = image, bandwidth = bandwidth[k], degree = degree, blur = blur, plot = TRUE)
      title(sub = paste0("bandwidth=", bandwidth[k]), line = 1.2)
      matout[k, ] <- quantile(diffout, probs = probs)
    }
  } else { # roof edge
    for (k in seq_along(bandwidth)) {
      diffout <- roofDiff(image = image, bandwidth = bandwidth[k], blur = blur) # roof edge detection statistics do not have degree or plot
      image(diffout, col = gray(c(0:255)/255))
      title(sub = paste0("bandwidth=", bandwidth[k]), line = 1.2)
      matout[k, ] <- quantile(diffout, probs = probs)
    }
  }
  return(matout)
}
