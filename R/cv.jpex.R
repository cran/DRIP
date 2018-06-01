cv.jpex = function(image, bandwidth) {
    if (!is.matrix(image)) {
        stop("image data must be a matrix")
    }
    else {
        n1 = dim(image)[1]
        n2 = dim(image)[2]
    }
    if (n1 != n2)
        stop("image data must be a square matrix")
    if (!is.numeric(bandwidth))
        stop("bandwidth must be numeric")
    n1 = dim(image)[1]
    z = matrix(as.double(image), ncol = n1)
    LLK = matrix(as.double(0), n1, n1)
    nband = length(bandwidth)
    cvs = rep(as.double(0), nband)
    for (iband in 1:nband) {
        cv = cvs[iband]
        k = as.integer(bandwidth[iband])
        out = .C("LOOCV", Zin=z, nin=n1, kin=k, LLKin=LLK, cv=cv, PACKAGE="DRIP")
        cvs[iband] = out$cv
    }
    band.min = as.integer(bandwidth[which.min(cvs)])
    out = .C("LOOCV", Zin=z, nin=n1, kin=band.min, LLKin=LLK, cv=cv, PACKAGE="DRIP")
    LLK = out$LLKin
    sigma = sqrt(mean((z - LLK)^2))
    return(list(LLK=LLK, sigma=sigma, cv=cvs, bandwidth=bandwidth, band.min=band.min))
}
