# This is R source code for function 'dKQ' in the R package
# 'image'.

dKQ <- function(set1, set2){
  if (!is.matrix(set1) || !is.matrix(set2))
    stop('Both arguments must be matrices.')
  if ((dim(set1)[1] != dim(set2)[1]) ||
      (dim(set1)[2] != dim(set2)[2]))
    stop('Two matrices must be of the same size.')
  if (dim(set1)[1] != dim(set1)[2])
    stop('Argument must be a square matrix.')
  if ((length(set1[(set1 != 0) & (set1 != 1)]) >= 1) ||
      (length(set2[(set2 != 0) & (set2 != 1)]) >= 1))
    stop('Both matrices can only have entries of 0 or 1.')
  n <- dim(set1)[1]
  out <- .Fortran(C_d_kq, n = as.integer(n - 1), edge1 =
    matrix(as.integer(set1), ncol = n), edge2 =
    matrix(as.integer(set2), ncol = n), dKQ = as.double(100))
  return(out$dKQ)
}

