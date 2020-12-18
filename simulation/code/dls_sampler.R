library(invgamma)

sampleRanEff <- function (ylist, fixefList, ranefList, rmatList, lmat, errVar, fixSamp) {
    nsample <- length(fixefList)
    nranef <- ncol(ranefList[[1]])

    ranSampList <- vector("list", nsample)
    for (ii in 1:nsample) {
        tmp1 <- ranefList[[ii]] %*% lmat
        tmp2 <- chol2inv(chol(tcrossprod(tmp1, tmp1) + errVar * rmatList[[ii]]))
        postRanVar <- diag(1, nranef) - crossprod(tmp1, tmp2 %*% tmp1)
        postRanMean <- drop(crossprod(tmp1, tmp2 %*% (ylist[[ii]] - fixefList[[ii]] %*% fixSamp)))
        ranSampList[[ii]] <- postRanMean + drop(crossprod(chol(postRanVar), rnorm(length(postRanMean))))
    }

    ranSampList
}

sampleFixEff <- function (ylist, fixefList, ranefList, rmatList, npart, dmat, errVar, muBeta0, sigBetaInv0) {
    nsample <- length(fixefList)

    umatList <- vector("list", nsample)
    umatInvList <- vector("list", nsample)
    xtuxList <- vector("list", nsample)
    xtuyList <- vector("list", nsample)
    for (ii in 1:nsample) {
        umatInvList[[ii]] <- errVar * rmatList[[ii]] + tcrossprod(ranefList[[ii]] %*% dmat, ranefList[[ii]])
        umatList[[ii]] <- chol2inv(chol(umatInvList[[ii]]))
        xtuxList[[ii]] <- crossprod(fixefList[[ii]], umatList[[ii]] %*% fixefList[[ii]])
        xtuyList[[ii]] <- crossprod(fixefList[[ii]], umatList[[ii]] %*% ylist[[ii]])
    }
    sumXtux <- Reduce("+", xtuxList)
    sumXtuy <- Reduce("+", xtuyList)
    ## the following code adds the power to variance
    postFixVar <- chol2inv(chol(sigBetaInv0 + npart * sumXtux))
    ## the following code adds the power to mean
    postFixMean <- as.numeric(postFixVar %*% (crossprod(sigBetaInv0, muBeta0) + npart * sumXtuy))

    fixSampList <- postFixMean + drop(crossprod(chol(postFixVar), rnorm(length(postFixMean))))

    fixSampList
}

sampleLmat <- function (ylist, fixefList, ranefList, rmatList, npart, errVar, fixSamp, ranSampList, muL0, sigLInv0) {
    nsample <- length(fixefList)
    nranef <- ncol(ranefList[[1]])

    mnList <- list()
    cvList <- list()
    for (ii in 1:nsample) {
        idz <- unlist(lapply(1:nranef, function(x) seq(x, nranef)))
        zz <- ranefList[[ii]][ , idz, drop = FALSE]
        idc <- unlist(lapply(1:nranef, function(x) rep(x, nranef - x + 1)))
        cc <- matrix(ranSampList[[ii]][idc], ncol = length(idc), nrow = nrow(zz), byrow = TRUE)
        ztilde <- zz * cc
        rr <- ylist[[ii]] - fixefList[[ii]] %*% fixSamp
        tmp1 <- solve(rmatList[[ii]])
        cvList[[ii]] <- crossprod(ztilde, tmp1 %*% ztilde)
        mnList[[ii]] <- crossprod(ztilde, tmp1 %*% rr)
    }

    ## the following code adds the power to mean
    tmp2 <- chol2inv(chol(npart * Reduce("+", cvList) + errVar * sigLInv0))
    postlmatVar <- errVar*tmp2
    postlmatMean <- drop(tmp2 %*% (npart * Reduce("+", mnList) + errVar * sigLInv0 %*% muL0))

    lmat <- diag(0.0, nranef)
    lmat[lower.tri(lmat, diag = TRUE)] <- postlmatMean + drop(crossprod(chol(postlmatVar), rnorm(length(postlmatMean))))

    lmat
}

sampleErrVar <- function (ylist, fixefList, ranefList, rmatList, npart, lmat, fixSamp, ranSampList, sig0, nu0) {
    nsample <- length(fixefList)
    nobs <- sum(sapply(ranefList, nrow))

    ## the following code adds the power to aa
    aa <- 0.5 * (npart * nobs + nu0)

    trm1 <- 0.0
    for (ii in 1:nsample) {
        resids <- ylist[[ii]] - fixefList[[ii]] %*% fixSamp - ranefList[[ii]] %*% lmat %*% ranSampList[[ii]]
        trm1 <- trm1 + drop(crossprod(resids, solve(rmatList[[ii]]) %*% resids))
    }
    trm2 <- sig0 * nu0

    ## the following code adds the power to bb
    bb <- 0.5 * (npart * trm1 + trm2)

    rinvgamma(1, aa, bb)
}

waspLmeSampler <- function (ylist, fixefList, ranefList, rmatList,npart,
                            niter = 1000, nburn = 500, nthin = 2,
                            dmat0, errVar0, muBeta0, sigBetaInv0, muL0, sigLInv0, sig0, nu0) {
    dmat <- dmat0
    lmat <- chol(dmat0)
    errVar <- errVar0
    fixs <- muBeta0
    its <- cnt <- 0
    nsample <- length(ranefList)
    nranef <- ncol(ranefList[[1]])

    sampLmat <- list()
    sampDmat <- list()
    sampErrVar <- list()
    sampBeta <- list()
    sampRans <- list()

    startTime <- proc.time()
    while (its < niter) {
        its <- its + 1

        if (its %% 10 == 0) cat("mixef iter: ", its, "\n")

        rans <- sampleRanEff(ylist, fixefList, ranefList, rmatList, lmat, errVar, fixs)
        fixs <- sampleFixEff(ylist, fixefList, ranefList, rmatList, npart, dmat, errVar, muBeta0, sigBetaInv0)
        lmat <- sampleLmat(ylist, fixefList, ranefList, rmatList, npart, errVar, fixs, rans, muL0, sigLInv0)
        dmat <- tcrossprod(lmat, lmat)
        errVar <- sampleErrVar(ylist, fixefList, ranefList, rmatList, npart, lmat, fixs, rans, sig0, nu0)

        if ((its > nburn) & (its %% nthin == 0)) {
            cnt <- cnt + 1
            sampLmat[[cnt]] <- lmat
            sampDmat[[cnt]] <- dmat
            sampErrVar[[cnt]] <- errVar
            sampBeta[[cnt]] <- fixs
        }
    }
    endTime <- proc.time()

    list(
        errVar = unlist(sampErrVar),
        dmat = sampDmat,
        lmat = sampLmat,
        beta = do.call(rbind, sampBeta),
        iter = its,
        time = endTime[3] - startTime[3]
    )
}
