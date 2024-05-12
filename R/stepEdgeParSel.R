# This is R source code for function 'stepEdgeParSel' in the DRIP
# package.
# Date: April 30, 2024
# Creator: Yicheng Kang

stepEdgeParSel <- function(image, bandwidth, thresh, nboot, 
                           degree = 1, blur = FALSE){
  degree <- round(degree)
  stopifnot((degree == 0 || degree == 1))
  if (degree == 0) {
    if (blur) {
      out <- stepEdgeParSelLC2K(image = image, bandwidth = bandwidth, 
                      thresh = thresh, nboot = nboot)
    } else { # no blur
      out <- stepEdgeParSelLCK(image = image, bandwidth = bandwidth, 
                     thresh = thresh, nboot = nboot)
    }
  } else { # degree = 1
    if (blur) {
      out <- stepEdgeParSelLL2K(image = image, bandwidth = bandwidth, 
                      thresh = thresh, nboot = nboot)
    } else { # no blur
      out <- stepEdgeParSelLLK(image = image, bandwidth = bandwidth, 
                     thresh = thresh, nboot = nboot)
    }
  }
  return(out)
}
