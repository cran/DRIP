# This is R source code for function 'stepEdg' in the DRIP
# package.
# Date: April 30, 2024
# Creator: Yicheng Kang

stepEdge <- function(image, bandwidth, thresh, degree = 1, 
                     blur = FALSE, plot = FALSE){
  degree <- round(degree)
  stopifnot((degree == 0 || degree == 1))
  if (degree == 0) {
    if (blur) {
      out <- stepEdgeLC2K(image = image, bandwidth = bandwidth, 
                      thresh = thresh, plot = plot)
    } else { # no blur
      out <- stepEdgeLCK(image = image, bandwidth = bandwidth, 
                     thresh = thresh, plot = plot)
    }
  } else { # degree = 1
    if (blur) {
      out <- stepEdgeLL2K(image = image, bandwidth = bandwidth, 
                      thresh = thresh, plot = plot)
    } else { # no blur
      out <- stepEdgeLLK(image = image, bandwidth = bandwidth, 
                     thresh = thresh, plot = plot)
    }
  }
  return(out)
}