# This is the R souce code for function ``stepDiff'' in the DRIP 
# package.
# Creator: Yicheng Kang
# Date: 04/30/2024.

stepDiff <- function(image, bandwidth, degree = 1, blur = FALSE, 
                     plot = FALSE) {
  degree <- round(degree)
  stopifnot((degree == 0 || degree == 1))
  if (degree == 0) {
    if (blur) {
      out <- diffLC2K(image = image, bandwidth = bandwidth, 
                      plot = plot)
    } else { # no blur
      out <- diffLCK(image = image, bandwidth = bandwidth, 
                     plot = plot)
    }
  } else { # degree = 1
    if (blur) {
      out <- diffLL2K(image = image, bandwidth = bandwidth, 
                      plot = plot)
    } else { # no blur
      out <- diffLLK(image = image, bandwidth = bandwidth, 
                     plot = plot)
    }
  }
  return(out)
}
