# This is R source code for function 'stepEdgeParSel' in the DRIP
# package.
# Date: May 30, 2024
# Creator: Yicheng Kang

# S3 Class

foo <- list(dKQ = matrix(0.0, nrow = 2, ncol = 2), bandwidth = as.integer(3),
            threshold = 1.2)
class(foo) <- "Edge_Parameters"

# The print method for S3 class

print.Edge_Parameters <- function(x, type = "all", ...) {
  stopifnot(class(x) == "Edge_Parameters")
  if (!(type %in% c("matrix", "parameters", "all"))) {
    stop("Wrong type value")
  } else {
    if (type == "matrix") {
      cat("The bootstrap matrix:\n")
      print.table(x$dKQ)
    }
    if (type == "parameters") {
      cat("The selected bandwidth: ", x$bandwidth, "\n")
      cat("The selected threshold: ", x$threshold, "\n")
    }
    if (type == "all") {
      cat("The bootstrap matrix:\n")
      print.table(x$dKQ)
      cat("The selected bandwidth: ", x$bandwidth, "\n")
      cat("The selected threshold: ", x$threshold, "\n")
    }
  }
}

# The summary method for S3 class

summary.Edge_Parameters <- function(object, ...) {
  stopifnot(class(object) == "Edge_Parameters")
  cat("The bootstrap matrix:\n")
  print.table(object$dKQ)
}

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
  out1 <- list(dKQ = out$output_matrix,
               bandwidth = out$selected_bandwidth,
               threshold = out$selected_threshold)
  class(out1) <- "Edge_Parameters"
  return(out1)
}
