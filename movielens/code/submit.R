## mtd 0 is full-data sampling
## mtd 1 is subset sampling

cmdArgs <- commandArgs(trailingOnly = TRUE)

mtd <- as.numeric(cmdArgs[1])
id <- as.numeric(cmdArgs[2])

if (mtd == 0) {
    library(Rcpp)
    library(RcppArmadillo)
    sourceCpp("dls.cpp")

    cvtrain <- readRDS("/Shared/ssrivastva/lswasp/data/train_ml.rds")
    train <- cvtrain[[id]]
    rm(cvtrain)

    group <- train$group
    grpLbl <- sort(unique(group))
    ngroup <- length(grpLbl)
    ranefList0 <- list()
    fixefList0 <- list()
    groupList0 <- list()
    ylist0 <- list()
    grpIdx0 <- list()
    for (gg in 1:ngroup) {
        grpIdx0[[gg]] <- which(group == grpLbl[gg])
        ranefList0[[gg]] <- train$z[grpIdx0[[gg]], , drop = FALSE]
        fixefList0[[gg]] <- train$x[grpIdx0[[gg]], , drop = FALSE]
        ylist0[[gg]] <- train$y[grpIdx0[[gg]]]
        groupList0[[gg]] <- train$group[grpIdx0[[gg]]]
    }

    nranef <- 6; nfixef <- 6;
    z <- do.call(rbind, ranefList0)
    x <- do.call(rbind, fixefList0)
    y <- as.matrix(unlist(ylist0))
    groupsTrain <- as.matrix(unlist(groupList0))
    muBeta0 <- as.matrix(rep(0.0, nfixef))
    sigBetaInv0 <- diag(0.0001, nfixef)
    muL0 <- as.matrix(rep(0.0, nranef * (nranef + 1) / 2))
    sigLInv0 <- diag(0.0001, nranef * (nranef + 1) / 2)
    sig0 <- 0.0; nu0 <- (2 + nfixef)
    dmat0 <- diag(1, nranef)
    errVar0 <- 10

    res <- wasp_lme_sampler (y, x, z,
                             groupsTrain,
                             1.0,
                             muBeta0,
                             sigBetaInv0,
                             dmat0,
                             errVar0,
                             muL0, sigLInv0,
                             sig0, nu0,
                             3000, 1000, 3)

    fname <- paste0("/Shared/ssrivastva/lswasp/ml/result/ml_rep_", id, "_full.rds")
    saveRDS(res, fname)
} else if (mtd == 1) {
    library(Rcpp)
    library(RcppArmadillo)
    sourceCpp("dls.cpp")

    npart = 100
    cvs <- rep(1:10, each = npart)
    nsubs <- rep(1:npart, times = 10)

    cid <- cvs[id]
    sid <- nsubs[id]

    fname <- paste0("/Shared/ssrivastva/lswasp/ml/data/ml_data_rep_", cid, "_k_", npart, ".rds")
    dat <- readRDS(fname)
    train <- dat[[sid]]

    group <- train$group
    grpLbl <- sort(unique(group))
    ngroup <- length(grpLbl)
    ranefList0 <- list()
    fixefList0 <- list()
    groupList0 <- list()
    ylist0 <- list()
    grpIdx0 <- list()
    for (gg in 1:ngroup) {
        grpIdx0[[gg]] <- which(group == grpLbl[gg])
        ranefList0[[gg]] <- train$z[grpIdx0[[gg]], , drop = FALSE]
        fixefList0[[gg]] <- train$x[grpIdx0[[gg]], , drop = FALSE]
        ylist0[[gg]] <- train$y[grpIdx0[[gg]]]
        groupList0[[gg]] <- train$group[grpIdx0[[gg]]]
    }

    nranef <- 6; nfixef <- 6;
    z <- do.call(rbind, ranefList0)
    x <- do.call(rbind, fixefList0)
    y <- as.matrix(unlist(ylist0))
    groupsTrain <- as.matrix(unlist(groupList0))
    muBeta0 <- as.matrix(rep(0.0, nfixef))
    sigBetaInv0 <- diag(0.0001, nfixef)
    muL0 <- as.matrix(rep(0.0, nranef * (nranef + 1) / 2))
    sigLInv0 <- diag(0.0001, nranef * (nranef + 1) / 2)
    sig0 <- 0.0; nu0 <- (2 + nfixef)
    dmat0 <- diag(1, nranef)
    errVar0 <- 10

    res <- wasp_lme_sampler (y, x, z,
                             groupsTrain,
                             train$nrep,
                             muBeta0,
                             sigBetaInv0,
                             dmat0,
                             errVar0,
                             muL0, sigLInv0,
                             sig0, nu0,
                             10000, 5000, 5)

    fname <- paste0("/Shared/ssrivastva/lswasp/ml/result/ml_rep_", cid, "_nsub_", sid, "_k_", npart, ".rds")
    saveRDS(res, fname)
} else if (mtd == 2) {
    source("dls_comb.R")

    cvs <- 1:10
    cid <- cvs[id]

    res <- list()
    for (sid in 1:100) {
        fname <- paste0("/Shared/ssrivastva/lswasp/ml/result/ml_rep_", cid, "_nsub_", sid, "_k_", 100, ".rds")
        res[[sid]] <- readRDS(fname)
    }

    lList0 <- lapply(res, function(x) x$lmat)
    lList <- lapply(lList0,
                    function(x) {
                        lsamp <- matrix(NA, 1000, 21)
                        for (ii in 1:1000) {
                            xx <- diag(0.0, 6)
                            diag(xx) <- x[ii, 1:6]
                            xx[lower.tri(xx)] <- x[ii, -(1:6)]
                            xx <- xx * matrix(sign(diag(xx)), nrow = 6, ncol = 6, byrow = TRUE)
                            lsamp[ii, ] <- xx[lower.tri(xx, diag = TRUE)]
                        }
                        lsamp
                    }
                    )
    bList <- lapply(res, function(x) x$beta)

    stime <- proc.time()
    lres <- sampleLmat(lList)
    bres <- sampleBetas(bList)
    dcres <- sampleDmatCorrMat(lList)
    etime <- proc.time()

    rname <- paste0("/Shared/ssrivastva/lswasp/ml/result/res/ml_samps_rep_", cid, ".rds")
    fres <- list(lmat = lres, betas = bres, dmat = dcres,
                 time =  (mean(unlist(lapply(res, function(x) x$time))) + etime[3] - stime[3]) / 3600
                 )
    saveRDS(fres, rname)
} else {
    print("peace")
}

