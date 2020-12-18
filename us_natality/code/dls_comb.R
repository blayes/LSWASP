## meanList = list of vector of multivariate mean
## covList = list of cov matrices
## wts = rep(1 / K, K), where K is the # of subset posteriors
obtainLocationScatterBarycenter <- function(meanList, covList, wts = NULL) {
    library(matrixStats)
    library(expm)
    library(MASS)
    ncomp <- length(covList)
    ndim <- nrow(covList[[1]])

    baryMean <- rowMeans(do.call(cbind, meanList))
    if (is.null(wts)) wts <- rep(1, ncomp) / ncomp

    baryCov <- diag(1.0, ndim)
    barySd <- sqrtm(baryCov)
    err <- baryCov
    cnt <- 1
    while ((norm(err, type = "F") > 1e-6) & cnt < 500) {
        if (cnt %% 10 == 0)  cat("iter: ", cnt, "\n")
        sumMat <- matrix(0.0, nrow = ndim, ncol = ndim)
        for (ii in 1:ncomp) {
            tmp <- chol2inv(chol(baryCov))
            tmp2 <- sqrtm(covList[[ii]] %*% baryCov) %*% tmp
            sumMat <- sumMat + wts[ii] * (tmp2 + t(tmp2)) / 2
        }
        baryCovNew <- barySd %*% (sumMat %*% sumMat) %*% barySd
        err <- baryCov - baryCovNew
        baryCov <- baryCovNew
        barySd <- sqrtm(baryCov)
        cnt <- cnt + 1
    }

    list(mean = baryMean, cov = baryCov, sqrt = barySd)
}

sampleLmat <- function (lmatList) {
    library(matrixStats)
    library(expm)

    npart <- length(lmatList)
    meanLmat <- list()
    covLmat <- list()

    samps <- lapply(lmatList,
                    function(x) do.call(rbind, lapply(x, function(y) {diag(y) <- abs(diag(y)); y[lower.tri(y, diag = TRUE)]})))
    meanLmat <- lapply(samps, function(x) colMeans(x))
    covLmat <- lapply(samps, function(x) cov(x))

    resBary <- obtainLocationScatterBarycenter(meanLmat, covLmat)
    muLmat <- resBary$mean
    sigLmat <- resBary$cov
    sqrtSigLmat <- resBary$sqrt

    wtime1 <- proc.time()
    baryList <- list()
    for (ii in 1:npart) {
        tmp <- as(chol2inv(chol(covLmat[[ii]])), "symmetricMatrix")
        tmp1 <- matrix(meanLmat[[ii]], nrow = nrow(samps[[ii]]), ncol = ncol(samps[[ii]]), byrow = TRUE)
        centScaledSamps <- sqrtm(tmp) %*% (t(samps[[ii]] - tmp1))
        baryList[[ii]] <- t(muLmat + sqrtSigLmat %*% centScaledSamps)
    }
    wtime2 <- proc.time()

    xtime1 <- proc.time()
    xueLiang <- list()
    for (ii in 1:npart) {
        tmp1 <- matrix(meanLmat[[ii]], nrow = nrow(samps[[ii]]), ncol = ncol(samps[[ii]]), byrow = TRUE)
        tmp2 <- matrix(muLmat, nrow = nrow(samps[[ii]]), ncol = ncol(samps[[ii]]), byrow = TRUE)
        xueLiang[[ii]] <- samps[[ii]] + tmp2 - tmp1
    }
    xtime2 <- proc.time()

    list(mean = muLmat,
         cov = sigLmat,
         covAvg = Reduce("+", covLmat) / npart,
         wasp = do.call(rbind, baryList),
         xueLiang = do.call(rbind, xueLiang),
         xlTime = xtime2[3] - xtime1[3],
         wTime = wtime2[3] - wtime1[3]
         )
}

sampleBetas <- function (betaList) {
    library(matrixStats)
    library(expm)

    npart <- length(betaList)
    meanBetas <- list()
    covBetas <- list()

    meanBetas <- lapply(betaList, function(x) colMeans(x))
    covBetas <- lapply(betaList, function(x) cov(x))

    resBary <- obtainLocationScatterBarycenter(meanBetas, covBetas)
    muBetas <- resBary$mean
    sigBetas <- resBary$cov
    sqrtSigBetas <- resBary$sqrt

    wtime1 <- proc.time()
    baryList <- list()
    for (ii in 1:npart) {
        tmp <- as(chol2inv(chol(covBetas[[ii]])), "symmetricMatrix")
        tmp1 <- matrix(meanBetas[[ii]], nrow = nrow(betaList[[ii]]), ncol = ncol(betaList[[ii]]), byrow = TRUE)
        centScaledSamps <- sqrtm(tmp) %*% (t(betaList[[ii]] - tmp1))
        baryList[[ii]] <- t(muBetas + sqrtSigBetas %*% centScaledSamps)
    }
    wtime2 <- proc.time()

    xtime1 <- proc.time()
    xueLiang <- list()
    for (ii in 1:npart) {
        tmp1 <- matrix(meanBetas[[ii]], nrow = nrow(betaList[[ii]]), ncol = ncol(betaList[[ii]]), byrow = TRUE)
        tmp2 <- matrix(muBetas, nrow = nrow(betaList[[ii]]), ncol = ncol(betaList[[ii]]), byrow = TRUE)
        xueLiang[[ii]] <- betaList[[ii]] + tmp2 - tmp1
    }
    xtime2 <- proc.time()

    list(mean = muBetas,
         cov = sigBetas,
         covAvg = Reduce("+", covBetas) / npart,
         wasp = do.call(rbind, baryList),
         xueLiang = do.call(rbind, xueLiang),
         xlTime = xtime2[3] - xtime1[3],
         wTime = wtime2[3] - wtime1[3]
         )
}

sampleDmatCorrMat <- function (lmatList, ndim = 3) {
    ## note that ndim is the dimension of Dmat
    lmatRes <- sampleLmat(lmatList)
    lmats <- lmatRes$wasp
    nsample <- nrow(lmats)

    wtime1 <- proc.time()
    dmats <- vector("list", nsample)
    corrs <- vector("list", nsample)
    tmp <- matrix(0.0, ncol = ndim, nrow = ndim)
    for (ii in 1:nsample) {
        tmp[lower.tri(tmp, diag = TRUE)] <- lmats[ii, ]
        dmats[[ii]] <- tcrossprod(tmp, tmp)
        corrs[[ii]] <- cov2cor(dmats[[ii]])
    }

    dvecs <- do.call(rbind, lapply(dmats, function(x) c(diag(x), sqrt(2) * x[lower.tri(x)])))
    rvecs <- do.call(rbind, lapply(corrs, function(x) c(diag(x), sqrt(2) * x[lower.tri(x)])))
    wtime2 <- proc.time()

    xtime1 <- proc.time()
    xllmat <- lmatRes$xueLiang
    xldmats <- vector("list", nsample)
    xlcorrs <- vector("list", nsample)
    tmp <- matrix(0.0, ncol = ndim, nrow = ndim)
    for (ii in 1:nsample) {
        tmp[lower.tri(tmp, diag = TRUE)] <- xllmat[ii, ]
        xldmats[[ii]] <- tcrossprod(tmp, tmp)
        xlcorrs[[ii]] <- cov2cor(xldmats[[ii]])
    }

    xldvecs <- do.call(rbind, lapply(xldmats, function(x) c(diag(x), sqrt(2) * x[lower.tri(x)])))
    xlrvecs <- do.call(rbind, lapply(xlcorrs, function(x) c(diag(x), sqrt(2) * x[lower.tri(x)])))
    xtime2 <- proc.time()

    list(dmat = dmats,
         corr = corrs,
         vecDmat = dvecs,
         vecRmat = rvecs,
         xldmat = xldmats,
         xlcorr = xlcorrs,
         xlvecDmat = xldvecs,
         xlvecRmat = xlrvecs,
         xlTime = xtime2[3] - xtime1[3],
         wTime = wtime2[3] - wtime1[3]
         )
}
