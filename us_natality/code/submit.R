# mtd 1-7 are sampling, mtd 8 is combination 
cmdArgs <- commandArgs(trailingOnly = TRUE)

mtd <- as.numeric(cmdArgs[1])
id <- as.numeric(cmdArgs[2])

if (mtd == 1) {
    source("dls_sampler.R")

    cvtrain <- readRDS("../data/train_bwt.rds")
    train <- cvtrain[[id]]
    rm(cvtrain)

    group <- train$group
    grpLbl <- sort(unique(group))
    ngroup <- length(grpLbl)
    ranefList0 <- list()
    fixefList0 <- list()
    rmatList0 <- list()
    ylist0 <- list()
    grpIdx0 <- list()
    for (gg in 1:ngroup) {
        grpIdx0[[gg]] <- which(group == grpLbl[gg])
        ranefList0[[gg]] <- train$z[grpIdx0[[gg]], , drop = FALSE]
        fixefList0[[gg]] <- train$x[grpIdx0[[gg]], , drop = FALSE]
        ylist0[[gg]] <- train$y[grpIdx0[[gg]]]
        rmatList0[[gg]] <- diag(1, length(grpIdx0[[gg]]))
    }

    nranef0 <- 3; nfixef0 <- 14;
    muBeta0 <- rep(0, nfixef0); sigBetaInv0 <- diag(0, nfixef0); nu0 <- (2 + nfixef0); sig0 <- 0;
    eta0 <- -(nranef0 + 1); tmat0 <- diag(0, nranef0);
    dmat0 <- diag(1, nranef0);
    errVar0 <- 10;
    sigLInv0 <- diag(0.01, nranef0 * (nranef0 + 1) / 2);
    muL0 <- rep(0.0, nranef0 * (nranef0 + 1) / 2)

    res <- waspLmeSampler(ylist0, fixefList0, ranefList0, rmatList0, 1,
                          niter = 10000, nburn = 5000, nthin = 5,
                          dmat0, errVar0, muBeta0, sigBetaInv0, muL0, sigLInv0, sig0, nu0)

    fname <- paste0("~/us_natality/result/real_full_rep_", id, ".rds")
    saveRDS(res, fname)

} else if (mtd == 2) {
    source("dls_sampler.R")

    cvs <- rep(1:10, each = 20)
    nsubs <- rep(1:20, times = 10)
    
    cid <- cvs[id]
    sid <- nsubs[id]

    fname <- paste0("../data/part_real_data_rep_", cid, "_k_20.rds")
    dat <- readRDS(fname)
    train <- dat[[sid]]

    group <- train$group
    grpLbl <- sort(unique(group))
    ngroup <- length(grpLbl)
    ranefList0 <- list()
    fixefList0 <- list()
    rmatList0 <- list()
    ylist0 <- list()
    grpIdx0 <- list()
    for (gg in 1:ngroup) {
        grpIdx0[[gg]] <- which(group == grpLbl[gg])
        ranefList0[[gg]] <- train$z[grpIdx0[[gg]], , drop = FALSE]
        fixefList0[[gg]] <- train$x[grpIdx0[[gg]], , drop = FALSE]
        ylist0[[gg]] <- train$y[grpIdx0[[gg]]]
        rmatList0[[gg]] <- diag(1, length(grpIdx0[[gg]]))
    }

    nranef0 <- 3; nfixef0 <- 14;
    muBeta0 <- rep(0, nfixef0); sigBetaInv0 <- diag(0, nfixef0); nu0 <- (2 + nfixef0); sig0 <- 0;
    eta0 <- -(nranef0 + 1); tmat0 <- diag(0, nranef0);
    dmat0 <- diag(1, nranef0);
    errVar0 <- 10;
    sigLInv0 <- diag(0.01, nranef0 * (nranef0 + 1) / 2);
    muL0 <- rep(0.0, nranef0 * (nranef0 + 1) / 2)

    res <- waspLmeSampler(ylist0, fixefList0, ranefList0, rmatList0, train$nrep,
                          niter = 10000, nburn = 5000, nthin = 5,
                          dmat0, errVar0, muBeta0, sigBetaInv0, muL0, sigLInv0, sig0, nu0)

    fname <- paste0("~/us_natality/result/samp/sub20/real_dls_rep_", cid, "_nsub_", sid, "_k_20.rds")
    saveRDS(res, fname)

} else if (mtd == 3) {
    source("dls_sampler.R")

    cvs <- rep(1:10, each = 30)
    nsubs <- rep(1:30, times = 10)

    cid <- cvs[id]
    sid <- nsubs[id]

    fname <- paste0("../data/part_real_data_rep_", cid, "_k_30.rds")
    dat <- readRDS(fname)
    train <- dat[[sid]]

    group <- train$group
    grpLbl <- sort(unique(group))
    ngroup <- length(grpLbl)
    ranefList0 <- list()
    fixefList0 <- list()
    rmatList0 <- list()
    ylist0 <- list()
    grpIdx0 <- list()
    for (gg in 1:ngroup) {
        grpIdx0[[gg]] <- which(group == grpLbl[gg])
        ranefList0[[gg]] <- train$z[grpIdx0[[gg]], , drop = FALSE]
        fixefList0[[gg]] <- train$x[grpIdx0[[gg]], , drop = FALSE]
        ylist0[[gg]] <- train$y[grpIdx0[[gg]]]
        rmatList0[[gg]] <- diag(1, length(grpIdx0[[gg]]))
    }

    nranef0 <- 3; nfixef0 <- 14;
    muBeta0 <- rep(0, nfixef0); sigBetaInv0 <- diag(0, nfixef0); nu0 <- (2 + nfixef0); sig0 <- 0;
    eta0 <- -(nranef0 + 1); tmat0 <- diag(0, nranef0);
    dmat0 <- diag(1, nranef0);
    errVar0 <- 10;
    sigLInv0 <- diag(0.01, nranef0 * (nranef0 + 1) / 2);
    muL0 <- rep(0.0, nranef0 * (nranef0 + 1) / 2)

    res <- waspLmeSampler(ylist0, fixefList0, ranefList0, rmatList0, train$nrep,
                          niter = 10000, nburn = 5000, nthin = 5,
                          dmat0, errVar0, muBeta0, sigBetaInv0, muL0, sigLInv0, sig0, nu0)

    fname <- paste0("~/us_natality/result/samp/sub30/real_dls_rep_", cid, "_nsub_", sid, "_k_30.rds")
    saveRDS(res, fname)

} else if (mtd == 4) {
  source("dls_sampler.R")
  
  cvs <- rep(1:10, each = 20)
  nsubs <- rep(1:20, times = 10)
  
  cid <- cvs[id]
  sid <- nsubs[id]
  
  fname <- paste0("../data/part_real_data_rep_", cid, "_k_20_samp250.rds")
  dat <- readRDS(fname)
  train <- dat[[sid]]
  
  group <- train$group
  grpLbl <- sort(unique(group))
  ngroup <- length(grpLbl)
  ranefList0 <- list()
  fixefList0 <- list()
  rmatList0 <- list()
  ylist0 <- list()
  grpIdx0 <- list()
  for (gg in 1:ngroup) {
    grpIdx0[[gg]] <- which(group == grpLbl[gg])
    ranefList0[[gg]] <- train$z[grpIdx0[[gg]], , drop = FALSE]
    fixefList0[[gg]] <- train$x[grpIdx0[[gg]], , drop = FALSE]
    ylist0[[gg]] <- train$y[grpIdx0[[gg]]]
    rmatList0[[gg]] <- diag(1, length(grpIdx0[[gg]]))
  }
  
  nranef0 <- 3; nfixef0 <- 14;
  muBeta0 <- rep(0, nfixef0); sigBetaInv0 <- diag(0, nfixef0); nu0 <- (2 + nfixef0); sig0 <- 0;
  eta0 <- -(nranef0 + 1); tmat0 <- diag(0, nranef0);
  dmat0 <- diag(1, nranef0);
  errVar0 <- 10;
  sigLInv0 <- diag(0.01, nranef0 * (nranef0 + 1) / 2);
  muL0 <- rep(0.0, nranef0 * (nranef0 + 1) / 2)
  
  res <- waspLmeSampler(ylist0, fixefList0, ranefList0, rmatList0, train$nrep,
                        niter = 10000, nburn = 5000, nthin = 5,
                        dmat0, errVar0, muBeta0, sigBetaInv0, muL0, sigLInv0, sig0, nu0)
  
  fname <- paste0("~/us_natality/result/samp250/sub20/real_dls_rep_", cid, "_nsub_", sid, "_k_20_samp250.rds")
  saveRDS(res, fname)

} else if (mtd == 5) {
  source("dls_sampler.R")
  
  cvs <- rep(1:10, each = 30)
  nsubs <- rep(1:30, times = 10)
  
  cid <- cvs[id]
  sid <- nsubs[id]
  
  fname <- paste0("../data/part_real_data_rep_", cid, "_k_30_samp250.rds")
  dat <- readRDS(fname)
  train <- dat[[sid]]
  
  group <- train$group
  grpLbl <- sort(unique(group))
  ngroup <- length(grpLbl)
  ranefList0 <- list()
  fixefList0 <- list()
  rmatList0 <- list()
  ylist0 <- list()
  grpIdx0 <- list()
  for (gg in 1:ngroup) {
    grpIdx0[[gg]] <- which(group == grpLbl[gg])
    ranefList0[[gg]] <- train$z[grpIdx0[[gg]], , drop = FALSE]
    fixefList0[[gg]] <- train$x[grpIdx0[[gg]], , drop = FALSE]
    ylist0[[gg]] <- train$y[grpIdx0[[gg]]]
    rmatList0[[gg]] <- diag(1, length(grpIdx0[[gg]]))
  }
  
  nranef0 <- 3; nfixef0 <- 14;
  muBeta0 <- rep(0, nfixef0); sigBetaInv0 <- diag(0, nfixef0); nu0 <- (2 + nfixef0); sig0 <- 0;
  eta0 <- -(nranef0 + 1); tmat0 <- diag(0, nranef0);
  dmat0 <- diag(1, nranef0);
  errVar0 <- 10;
  sigLInv0 <- diag(0.01, nranef0 * (nranef0 + 1) / 2);
  muL0 <- rep(0.0, nranef0 * (nranef0 + 1) / 2)
  
  res <- waspLmeSampler(ylist0, fixefList0, ranefList0, rmatList0, train$nrep,
                        niter = 10000, nburn = 5000, nthin = 5,
                        dmat0, errVar0, muBeta0, sigBetaInv0, muL0, sigLInv0, sig0, nu0)
  
  fname <- paste0("~/us_natality/result/samp250/sub30/real_dls_rep_", cid, "_nsub_", sid, "_k_30_samp250.rds")
  saveRDS(res, fname)

}else if (mtd == 6) {
  source("dls_sampler.R")
  
  cvs <- rep(1:10, each = 50)
  nsubs <- rep(1:50, times = 10)
  
  cid <- cvs[id]
  sid <- nsubs[id]
  
  fname <- paste0("../data/part_real_data_rep_", cid, "_k_50_samp250.rds")
  dat <- readRDS(fname)
  train <- dat[[sid]]
  
  group <- train$group
  grpLbl <- sort(unique(group))
  ngroup <- length(grpLbl)
  ranefList0 <- list()
  fixefList0 <- list()
  rmatList0 <- list()
  ylist0 <- list()
  grpIdx0 <- list()
  for (gg in 1:ngroup) {
    grpIdx0[[gg]] <- which(group == grpLbl[gg])
    ranefList0[[gg]] <- train$z[grpIdx0[[gg]], , drop = FALSE]
    fixefList0[[gg]] <- train$x[grpIdx0[[gg]], , drop = FALSE]
    ylist0[[gg]] <- train$y[grpIdx0[[gg]]]
    rmatList0[[gg]] <- diag(1, length(grpIdx0[[gg]]))
  }
  
  nranef0 <- 3; nfixef0 <- 14;
  muBeta0 <- rep(0, nfixef0); sigBetaInv0 <- diag(0, nfixef0); nu0 <- (2 + nfixef0); sig0 <- 0;
  eta0 <- -(nranef0 + 1); tmat0 <- diag(0, nranef0);
  dmat0 <- diag(1, nranef0);
  errVar0 <- 10;
  sigLInv0 <- diag(0.01, nranef0 * (nranef0 + 1) / 2);
  muL0 <- rep(0.0, nranef0 * (nranef0 + 1) / 2)
  
  res <- waspLmeSampler(ylist0, fixefList0, ranefList0, rmatList0, train$nrep,
                        niter = 10000, nburn = 5000, nthin = 5,
                        dmat0, errVar0, muBeta0, sigBetaInv0, muL0, sigLInv0, sig0, nu0)
  
  fname <- paste0("~/us_natalit/result/samp250/sub50/real_dls_rep_", cid, "_nsub_", sid, "_k_50_samp250.rds")
  saveRDS(res, fname)

}else if (mtd == 7) {
  source("dls_sampler.R")
  
  cvs <- rep(1:10, each = 100)
  nsubs <- rep(1:100, times = 10)
  
  cid <- cvs[id]
  sid <- nsubs[id]
  
  fname <- paste0("../data/part_real_data_rep_", cid, "_k_100_samp250.rds")
  dat <- readRDS(fname)
  train <- dat[[sid]]
  
  group <- train$group
  grpLbl <- sort(unique(group))
  ngroup <- length(grpLbl)
  ranefList0 <- list()
  fixefList0 <- list()
  rmatList0 <- list()
  ylist0 <- list()
  grpIdx0 <- list()
  for (gg in 1:ngroup) {
    grpIdx0[[gg]] <- which(group == grpLbl[gg])
    ranefList0[[gg]] <- train$z[grpIdx0[[gg]], , drop = FALSE]
    fixefList0[[gg]] <- train$x[grpIdx0[[gg]], , drop = FALSE]
    ylist0[[gg]] <- train$y[grpIdx0[[gg]]]
    rmatList0[[gg]] <- diag(1, length(grpIdx0[[gg]]))
  }
  
  nranef0 <- 3; nfixef0 <- 14;
  muBeta0 <- rep(0, nfixef0); sigBetaInv0 <- diag(0, nfixef0); nu0 <- (2 + nfixef0); sig0 <- 0;
  eta0 <- -(nranef0 + 1); tmat0 <- diag(0, nranef0);
  dmat0 <- diag(1, nranef0);
  errVar0 <- 10;
  sigLInv0 <- diag(0.01, nranef0 * (nranef0 + 1) / 2);
  muL0 <- rep(0.0, nranef0 * (nranef0 + 1) / 2)
  
  res <- waspLmeSampler(ylist0, fixefList0, ranefList0, rmatList0, train$nrep,
                        niter = 10000, nburn = 5000, nthin = 5,
                        dmat0, errVar0, muBeta0, sigBetaInv0, muL0, sigLInv0, sig0, nu0)
  
  fname <- paste0("~/us_natality/result/samp250/sub100/real_dls_rep_", cid, "_nsub_", sid, "_k_100_samp250.rds")
  saveRDS(res, fname)

}else if (mtd == 8) {
  source("dls_sampler.R")
  
  cvs <- rep(1:10, each = 20)
  nsubs <- rep(1:20, times = 10)
  
  cid <- cvs[id]
  sid <- nsubs[id]
  
  fname <- paste0("../data/part_real_data_rep_", cid, "_k_20_samp500.rds")
  dat <- readRDS(fname)
  train <- dat[[sid]]
  
  group <- train$group
  grpLbl <- sort(unique(group))
  ngroup <- length(grpLbl)
  ranefList0 <- list()
  fixefList0 <- list()
  rmatList0 <- list()
  ylist0 <- list()
  grpIdx0 <- list()
  for (gg in 1:ngroup) {
    grpIdx0[[gg]] <- which(group == grpLbl[gg])
    ranefList0[[gg]] <- train$z[grpIdx0[[gg]], , drop = FALSE]
    fixefList0[[gg]] <- train$x[grpIdx0[[gg]], , drop = FALSE]
    ylist0[[gg]] <- train$y[grpIdx0[[gg]]]
    rmatList0[[gg]] <- diag(1, length(grpIdx0[[gg]]))
  }
  
  nranef0 <- 3; nfixef0 <- 14;
  muBeta0 <- rep(0, nfixef0); sigBetaInv0 <- diag(0, nfixef0); nu0 <- (2 + nfixef0); sig0 <- 0;
  eta0 <- -(nranef0 + 1); tmat0 <- diag(0, nranef0);
  dmat0 <- diag(1, nranef0);
  errVar0 <- 10;
  sigLInv0 <- diag(0.01, nranef0 * (nranef0 + 1) / 2);
  muL0 <- rep(0.0, nranef0 * (nranef0 + 1) / 2)
  
  res <- waspLmeSampler(ylist0, fixefList0, ranefList0, rmatList0, train$nrep,
                        niter = 10000, nburn = 5000, nthin = 5,
                        dmat0, errVar0, muBeta0, sigBetaInv0, muL0, sigLInv0, sig0, nu0)
  
  fname <- paste0("~/us_natality/result/samp500/sub20/real_dls_rep_", cid, "_nsub_", sid, "_k_20_samp500.rds")
  saveRDS(res, fname)

} else if (mtd == 9) {
  source("dls_sampler.R")
  
  cvs <- rep(1:10, each = 30)
  nsubs <- rep(1:30, times = 10)
  
  cid <- cvs[id]
  sid <- nsubs[id]
  
  fname <- paste0("../data/part_real_data_rep_", cid, "_k_30_samp500.rds")
  dat <- readRDS(fname)
  train <- dat[[sid]]
  
  group <- train$group
  grpLbl <- sort(unique(group))
  ngroup <- length(grpLbl)
  ranefList0 <- list()
  fixefList0 <- list()
  rmatList0 <- list()
  ylist0 <- list()
  grpIdx0 <- list()
  for (gg in 1:ngroup) {
    grpIdx0[[gg]] <- which(group == grpLbl[gg])
    ranefList0[[gg]] <- train$z[grpIdx0[[gg]], , drop = FALSE]
    fixefList0[[gg]] <- train$x[grpIdx0[[gg]], , drop = FALSE]
    ylist0[[gg]] <- train$y[grpIdx0[[gg]]]
    rmatList0[[gg]] <- diag(1, length(grpIdx0[[gg]]))
  }
  
  nranef0 <- 3; nfixef0 <- 14;
  muBeta0 <- rep(0, nfixef0); sigBetaInv0 <- diag(0, nfixef0); nu0 <- (2 + nfixef0); sig0 <- 0;
  eta0 <- -(nranef0 + 1); tmat0 <- diag(0, nranef0);
  dmat0 <- diag(1, nranef0);
  errVar0 <- 10;
  sigLInv0 <- diag(0.01, nranef0 * (nranef0 + 1) / 2);
  muL0 <- rep(0.0, nranef0 * (nranef0 + 1) / 2)
  
  res <- waspLmeSampler(ylist0, fixefList0, ranefList0, rmatList0, train$nrep,
                        niter = 10000, nburn = 5000, nthin = 5,
                        dmat0, errVar0, muBeta0, sigBetaInv0, muL0, sigLInv0, sig0, nu0)
  
  fname <- paste0("~/us_natality/result/samp500/sub30/real_dls_rep_", cid, "_nsub_", sid, "_k_30_samp500.rds")
  saveRDS(res, fname)

}else if (mtd == 10) {
  source("dls_sampler.R")
  
  cvs <- rep(1:10, each = 50)
  nsubs <- rep(1:50, times = 10)
  
  cid <- cvs[id]
  sid <- nsubs[id]
  
  fname <- paste0("../data/part_real_data_rep_", cid, "_k_50_samp500.rds")
  dat <- readRDS(fname)
  train <- dat[[sid]]
  
  group <- train$group
  grpLbl <- sort(unique(group))
  ngroup <- length(grpLbl)
  ranefList0 <- list()
  fixefList0 <- list()
  rmatList0 <- list()
  ylist0 <- list()
  grpIdx0 <- list()
  for (gg in 1:ngroup) {
    grpIdx0[[gg]] <- which(group == grpLbl[gg])
    ranefList0[[gg]] <- train$z[grpIdx0[[gg]], , drop = FALSE]
    fixefList0[[gg]] <- train$x[grpIdx0[[gg]], , drop = FALSE]
    ylist0[[gg]] <- train$y[grpIdx0[[gg]]]
    rmatList0[[gg]] <- diag(1, length(grpIdx0[[gg]]))
  }
  
  nranef0 <- 3; nfixef0 <- 14;
  muBeta0 <- rep(0, nfixef0); sigBetaInv0 <- diag(0, nfixef0); nu0 <- (2 + nfixef0); sig0 <- 0;
  eta0 <- -(nranef0 + 1); tmat0 <- diag(0, nranef0);
  dmat0 <- diag(1, nranef0);
  errVar0 <- 10;
  sigLInv0 <- diag(0.01, nranef0 * (nranef0 + 1) / 2);
  muL0 <- rep(0.0, nranef0 * (nranef0 + 1) / 2)
  
  res <- waspLmeSampler(ylist0, fixefList0, ranefList0, rmatList0, train$nrep,
                        niter = 10000, nburn = 5000, nthin = 5,
                        dmat0, errVar0, muBeta0, sigBetaInv0, muL0, sigLInv0, sig0, nu0)
  
  fname <- paste0("~/us_natality/result/samp500/sub50/real_dls_rep_", cid, "_nsub_", sid, "_k_50_samp500.rds")
  saveRDS(res, fname)

}else if (mtd == 11) {
  source("dls_sampler.R")
  
  cvs <- rep(1:10, each = 100)
  nsubs <- rep(1:100, times = 10)
  
  cid <- cvs[id]
  sid <- nsubs[id]
  
  fname <- paste0("../data/part_real_data_rep_", cid, "_k_100_samp500.rds")
  dat <- readRDS(fname)
  train <- dat[[sid]]
  
  group <- train$group
  grpLbl <- sort(unique(group))
  ngroup <- length(grpLbl)
  ranefList0 <- list()
  fixefList0 <- list()
  rmatList0 <- list()
  ylist0 <- list()
  grpIdx0 <- list()
  for (gg in 1:ngroup) {
    grpIdx0[[gg]] <- which(group == grpLbl[gg])
    ranefList0[[gg]] <- train$z[grpIdx0[[gg]], , drop = FALSE]
    fixefList0[[gg]] <- train$x[grpIdx0[[gg]], , drop = FALSE]
    ylist0[[gg]] <- train$y[grpIdx0[[gg]]]
    rmatList0[[gg]] <- diag(1, length(grpIdx0[[gg]]))
  }
  
  nranef0 <- 3; nfixef0 <- 14;
  muBeta0 <- rep(0, nfixef0); sigBetaInv0 <- diag(0, nfixef0); nu0 <- (2 + nfixef0); sig0 <- 0;
  eta0 <- -(nranef0 + 1); tmat0 <- diag(0, nranef0);
  dmat0 <- diag(1, nranef0);
  errVar0 <- 10;
  sigLInv0 <- diag(0.01, nranef0 * (nranef0 + 1) / 2);
  muL0 <- rep(0.0, nranef0 * (nranef0 + 1) / 2)
  
  res <- waspLmeSampler(ylist0, fixefList0, ranefList0, rmatList0, train$nrep,
                        niter = 10000, nburn = 5000, nthin = 5,
                        dmat0, errVar0, muBeta0, sigBetaInv0, muL0, sigLInv0, sig0, nu0)
  
  fname <- paste0("~/us_natality/result/samp500/sub100/real_dls_rep_", cid, "_nsub_", sid, "_k_100_samp500.rds")
  saveRDS(res, fname)

}else if (mtd == 12) {
    source("dls_comb.R")
    
    cvs <- 1:10
    cid <- cvs[id]
    
    for (nsub in c(20, 30)) {
      res <- list()
      for (sid in 1:nsub) {
        fname <- paste0("~/us_natality/result/samp/sub",nsub,"/real_dls_rep_", cid, "_nsub_", sid, "_k_", nsub, ".rds")
        res[[sid]] <- readRDS(fname)
      }
      
      lList <- lapply(res, function(x) x$lmat)
      bList <- lapply(res, function(x) x$beta)
      
      stime <- proc.time()
      lres <- sampleLmat(lList)
      bres <- sampleBetas(bList)
      dcres <- sampleDmatCorrMat(lList, ndim = 3)
      etime <- proc.time()
      
      rname <- paste0("~/us_natality/result/comb/real_dls_samps_rep_", cid, "_k_", nsub, ".rds")
      fres <- list(lmat = lres, betas = bres, dmat = dcres,
                   time =  mean(unlist(lapply(res, function(x) x$time))) + etime[3] - stime[3]
      )
      saveRDS(fres, rname)
    }
      for (nsub in c(20,30,50,100)) {
        res <- list()
        for (sid in 1:nsub) {
          fname <- paste0("~/us_natality/result/samp/sub",nsub,"/real_dls_rep_", cid, "_nsub_", sid, "_k_", nsub, "_samp250.rds")
          res[[sid]] <- readRDS(fname)
        }
        
        lList <- lapply(res, function(x) x$lmat)
        bList <- lapply(res, function(x) x$beta)
        
        stime <- proc.time()
        lres <- sampleLmat(lList)
        bres <- sampleBetas(bList)
        dcres <- sampleDmatCorrMat(lList, ndim = 3)
        etime <- proc.time()

        rname <- paste0("~/us_natality/result/comb/real_dls_samps_rep_", cid, "_k_", nsub, "_samp250.rds")
        fres <- list(lmat = lres, betas = bres, dmat = dcres,
                    time =  mean(unlist(lapply(res, function(x) x$time))) + etime[3] - stime[3]
        )
        saveRDS(fres, rname)
      }   
        
        for (nsub in c(20,30,50,100)) {
          res <- list()
          for (sid in 1:nsub) {
            fname <- paste0("~/us_natality/result/samp/sub",nsub,"/real_dls_rep_", cid, "_nsub_", sid, "_k_", nsub, "_samp500.rds")
            res[[sid]] <- readRDS(fname)
          }
          
          lList <- lapply(res, function(x) x$lmat)
          bList <- lapply(res, function(x) x$beta)
          
          stime <- proc.time()
          lres <- sampleLmat(lList)
          bres <- sampleBetas(bList)
          dcres <- sampleDmatCorrMat(lList, ndim = 3)
          etime <- proc.time()

          
          rname <- paste0("~/us_natality/result/comb/real_dls_samps_rep_", cid, "_k_", nsub, "_samp500.rds")
          fres <- list(lmat = lres, betas = bres, dmat = dcres,
                       time =  mean(unlist(lapply(res, function(x) x$time))) + etime[3] - stime[3]
          )
          saveRDS(fres, rname)
        }
  }else if (mtd == 13) {
    # DLS alogrithm 
    source("ganbo_accuracy.R")
    
    nsubset <- c(20, 30)
    subset <- nsubset[id]
    nrep <- 1:10
    
    #########  to get accuracy measure for disjoint partition #############
    resfull <- vector("list",10)
    
    for(ii in nrep){
      data_name <- paste0("../result/full/real_full_rep_",ii,".rds")
      resfull[[ii]] <- readRDS(data_name)
    }
    
    res <- vector("list",10)
    
    ## we import the same data because xl and ours are stored in the same dataframe 
    for(ii in nrep){
      data_name <- paste0("../result/comb/real_dls_samps_rep_",ii,"_k_", subset, ".rds")
      res[[ii]] <- readRDS(data_name)
    }
    
    nobs <- 20000 
    ## each for one replication in total 10
    accuracy_list <- vector("list",10)
    ## get the accuracy for the whole matrix of betas, lmat, Dmat, Rmat
    for(jj in 1:10){
      accur_beta_temp <- Accuracy(res[[jj]]$betas$wasp,resfull[[jj]]$beta, nobs)
      accur_lmat_temp <- Accuracy(res[[jj]]$lmat$wasp,resfull[[jj]]$vechlmat, nobs)
      accur_dmat_temp <- Accuracy(res[[jj]]$dmat$vecDmat,resfull[[jj]]$vecdmat, nobs)
      accur_Rmat_temp <- Accuracy(res[[jj]]$dmat$vecRmat[,4:6],resfull[[jj]]$vecrmat[,4:6], nobs)
      
      accuracy_list[[jj]]$accur_beta <-accur_beta_temp
      accuracy_list[[jj]]$accur_lmat <-accur_lmat_temp
      accuracy_list[[jj]]$accur_dmat <-accur_dmat_temp
      accuracy_list[[jj]]$accur_Rmat <-accur_Rmat_temp
      
      print(jj)
    }  
    
    ## get the accuracy for the selected matrix of betas, lmat, Dmat, Rmat
    for(jj in 1:10){
      ## select 2-dimensions, or 3-dimensions of 4-dimension beta matrix 
      for(jjj in 2:3){
        nD_beta <- combn(1:4, jjj)
        accur_beta_temp_nD <-0
        
        for(jjjj in ncol(nD_beta)){
          accur_beta_temp_nD <- Accuracy(res[[jj]]$betasNew$wasp[,nD_beta[,jjj]],resfull[[jj]]$beta[,nD_beta[,jjj]], nobs) +accur_beta_temp_nD
        }
        if(jjj ==2){
          accuracy_list[[jj]]$accur_beta_2D <-accur_beta_temp_nD/ncol(nD_beta)
        }else if(jjj==3){
          accuracy_list[[jj]]$accur_beta_3D <-accur_beta_temp_nD/ncol(nD_beta)
        }
      }
      ## select 2,3,4,5-dimensions of 6-dimension L matrix 
      for(jjj in 2:5){
        nD_lmat <- combn(1:6, jjj)
        accur_lmat_temp_nD <-0
        
        for(jjjj in ncol(nD_lmat)){
          accur_lmat_temp_nD <- Accuracy(res[[jj]]$lmatNew$wasp[,nD_lmat[,jjj]],resfull[[jj]]$vechlmat[,nD_lmat[,jjj]], nobs) +accur_lmat_temp_nD
        }
        if(jjj ==2){
          accuracy_list[[jj]]$accur_lmat_2D <-accur_lmat_temp_nD/ncol(nD_lmat)
        }else if(jjj==3){
          accuracy_list[[jj]]$accur_lmat_3D <-accur_lmat_temp_nD/ncol(nD_lmat)
        }else if(jjj==4){
          accuracy_list[[jj]]$accur_lmat_4D <-accur_lmat_temp_nD/ncol(nD_lmat)
        }else if(jjj==5){
          accuracy_list[[jj]]$accur_lmat_5D <-accur_lmat_temp_nD/ncol(nD_lmat)
        }
      }
      ## select 2,3,4,5-dimensions of 6-dimension D matrix 
      for(jjj in 2:5){
        nD_dmat <- combn(1:6, jjj)
        accur_dmat_temp_nD <-0
        
        for(jjjj in ncol(nD_dmat)){
          accur_dmat_temp_nD <- Accuracy(res[[jj]]$dmatNew$vecDmat[,nD_dmat[,jjj]],resfull[[jj]]$vecdmat[,nD_dmat[,jjj]], nobs) +accur_dmat_temp_nD
        }
        if(jjj ==2){
          accuracy_list[[jj]]$accur_dmat_2D <-accur_dmat_temp_nD/ncol(nD_dmat)
        }else if(jjj==3){
          accuracy_list[[jj]]$accur_dmat_3D <-accur_dmat_temp_nD/ncol(nD_dmat)
        }else if(jjj==4){
          accuracy_list[[jj]]$accur_dmat_4D <-accur_dmat_temp_nD/ncol(nD_dmat)
        }else if(jjj==5){
          accuracy_list[[jj]]$accur_dmat_5D <-accur_dmat_temp_nD/ncol(nD_dmat)
        }
      }
      print(jj)
    }     
    
    ## this will work because each var of the sub-list in accuracy_list only contains one item
    ## change into matrix by row will go through each var from 1st sub-list to the 10th sub-list
    ## forming ten rows 
    accuracy_matrix <- matrix(unlist(accuracy_list), ncol = 14,byrow=TRUE)
    colnames(accuracy_matrix) <- c("accur_beta", "accur_lmat", "accur_dmat", "accur_Rmat", "accur_beta_2D",
                                   "accur_beta_3D","accur_lmat_2D","accur_lmat_3D","accur_lmat_4D","accur_lmat_5D",
                                   "accur_dmat_2D","accur_dmat_3D","accur_dmat_4D","accur_dmat_5D")
    accuracy_mean <- apply(accuracy_matrix, 2, mean)
    accuracy_sd <- apply(accuracy_matrix, 2,sd)
    result <- vector("list",1)
    result[[1]]$accuracy_matrix <-accuracy_matrix
    result[[1]]$accuracy_mean <-accuracy_mean
    result[[1]]$accuracy_sd <-accuracy_sd
    
    fname <- paste0("../result/real_dls_k_", subset, "_accuracy_ind.rds")
    saveRDS(result, fname)

    
    #########  to get accuracy measure for overlap 250 #############
    
    resfull <- vector("list",10)
    
    for(ii in nrep){
      data_name <- paste0("../result/full/real_full_rep_",ii,".rds")
      resfull[[ii]] <- readRDS(data_name)
    }
    
    res <- vector("list",10)
    
    for(ii in nrep){
      data_name <- paste0("../result/comb/real_dls_samps_rep_",ii,"_k_", subset, "_samp250.rds")
      res[[ii]] <- readRDS(data_name)
    }
    
    nobs <- 20000 
    
    accuracy_list <- vector("list",10)
    
    for(jj in 1:10){
      accur_beta_temp <- Accuracy(res[[jj]]$betas$wasp,resfull[[jj]]$beta, nobs)
      accur_lmat_temp <- Accuracy(res[[jj]]$lmat$wasp,resfull[[jj]]$vechlmat, nobs)
      accur_dmat_temp <- Accuracy(res[[jj]]$dmat$vecDmat,resfull[[jj]]$vecdmat, nobs)
      accur_Rmat_temp <- Accuracy(res[[jj]]$dmat$vecRmat[,4:6],resfull[[jj]]$vecrmat[,4:6], nobs)
      
      accuracy_list[[jj]]$accur_beta <-accur_beta_temp
      accuracy_list[[jj]]$accur_lmat <-accur_lmat_temp
      accuracy_list[[jj]]$accur_dmat <-accur_dmat_temp
      accuracy_list[[jj]]$accur_Rmat <-accur_Rmat_temp
      
      print(jj)
    }  
    
    for(jj in 1:10){
      ## nD beta 
      for(jjj in 2:3){
        nD_beta <- combn(1:4, jjj)
        accur_beta_temp_nD <-0
        
        for(jjjj in ncol(nD_beta)){
          accur_beta_temp_nD <- Accuracy(res[[jj]]$betasNew$wasp[,nD_beta[,jjj]],resfull[[jj]]$beta[,nD_beta[,jjj]], nobs) +accur_beta_temp_nD
        }
        if(jjj ==2){
          accuracy_list[[jj]]$accur_beta_2D <-accur_beta_temp_nD/ncol(nD_beta)
        }else if(jjj==3){
          accuracy_list[[jj]]$accur_beta_3D <-accur_beta_temp_nD/ncol(nD_beta)
        }
      }
      ## nD lmat
      for(jjj in 2:5){
        nD_lmat <- combn(1:6, jjj)
        accur_lmat_temp_nD <-0
        
        for(jjjj in ncol(nD_lmat)){
          accur_lmat_temp_nD <- Accuracy(res[[jj]]$lmatNew$wasp[,nD_lmat[,jjj]],resfull[[jj]]$vechlmat[,nD_lmat[,jjj]], nobs) +accur_lmat_temp_nD
        }
        if(jjj ==2){
          accuracy_list[[jj]]$accur_lmat_2D <-accur_lmat_temp_nD/ncol(nD_lmat)
        }else if(jjj==3){
          accuracy_list[[jj]]$accur_lmat_3D <-accur_lmat_temp_nD/ncol(nD_lmat)
        }else if(jjj==4){
          accuracy_list[[jj]]$accur_lmat_4D <-accur_lmat_temp_nD/ncol(nD_lmat)
        }else if(jjj==5){
          accuracy_list[[jj]]$accur_lmat_5D <-accur_lmat_temp_nD/ncol(nD_lmat)
        }
      }
      ## nD dmat
      for(jjj in 2:5){
        nD_dmat <- combn(1:6, jjj)
        accur_dmat_temp_nD <-0
        
        for(jjjj in ncol(nD_dmat)){
          accur_dmat_temp_nD <- Accuracy(res[[jj]]$dmatNew$vecDmat[,nD_dmat[,jjj]],resfull[[jj]]$vecdmat[,nD_dmat[,jjj]], nobs) +accur_dmat_temp_nD
        }
        if(jjj ==2){
          accuracy_list[[jj]]$accur_dmat_2D <-accur_dmat_temp_nD/ncol(nD_dmat)
        }else if(jjj==3){
          accuracy_list[[jj]]$accur_dmat_3D <-accur_dmat_temp_nD/ncol(nD_dmat)
        }else if(jjj==4){
          accuracy_list[[jj]]$accur_dmat_4D <-accur_dmat_temp_nD/ncol(nD_dmat)
        }else if(jjj==5){
          accuracy_list[[jj]]$accur_dmat_5D <-accur_dmat_temp_nD/ncol(nD_dmat)
        }
      }
      print(jj)
    }      
    accuracy_matrix <- matrix(unlist(accuracy_list), ncol = 14,byrow=TRUE)
    colnames(accuracy_matrix) <- c("accur_beta", "accur_lmat", "accur_dmat", "accur_Rmat", "accur_beta_2D",
                                   "accur_beta_3D","accur_lmat_2D","accur_lmat_3D","accur_lmat_4D","accur_lmat_5D",
                                   "accur_dmat_2D","accur_dmat_3D","accur_dmat_4D","accur_dmat_5D")
    accuracy_mean <- apply(accuracy_matrix, 2, mean)
    accuracy_sd <- apply(accuracy_matrix, 2,sd)
    result <- vector("list",1)
    result[[1]]$accuracy_matrix <-accuracy_matrix
    result[[1]]$accuracy_mean <-accuracy_mean
    result[[1]]$accuracy_sd <-accuracy_sd
    
    fname <- paste0("../result/real_dls_k_", subset, "_accuracy_samp250.rds")
    saveRDS(result, fname)
    
    #########  to get accuracy measure for overlap 500 #############
    
    resfull <- vector("list",10)
    
    for(ii in nrep){
      data_name <- paste0("../result/full/real_full_rep_",ii,".rds")
      resfull[[ii]] <- readRDS(data_name)
    }
    
    res <- vector("list",10)
    
    for(ii in nrep){
      data_name <- paste0("../result/comb/real_dls_samps_rep_",ii,"_k_", subset, "_samp500.rds")
      res[[ii]] <- readRDS(data_name)
    }
    
    nobs <- 20000 
    
    accuracy_list <- vector("list",10)
    
    for(jj in 1:10){
      accur_beta_temp <- Accuracy(res[[jj]]$betas$wasp,resfull[[jj]]$beta, nobs)
      accur_lmat_temp <- Accuracy(res[[jj]]$lmat$wasp,resfull[[jj]]$vechlmat, nobs)
      accur_dmat_temp <- Accuracy(res[[jj]]$dmat$vecDmat,resfull[[jj]]$vecdmat, nobs)
      accur_Rmat_temp <- Accuracy(res[[jj]]$dmat$vecRmat[,4:6],resfull[[jj]]$vecrmat[,4:6], nobs)
      
      accuracy_list[[jj]]$accur_beta <-accur_beta_temp
      accuracy_list[[jj]]$accur_lmat <-accur_lmat_temp
      accuracy_list[[jj]]$accur_dmat <-accur_dmat_temp
      accuracy_list[[jj]]$accur_Rmat <-accur_Rmat_temp
      
      print(jj)
    }  
    
    for(jj in 1:10){
      ## nD beta 
      for(jjj in 2:3){
        nD_beta <- combn(1:4, jjj)
        accur_beta_temp_nD <-0
        
        for(jjjj in ncol(nD_beta)){
          accur_beta_temp_nD <- Accuracy(res[[jj]]$betasNew$wasp[,nD_beta[,jjj]],resfull[[jj]]$beta[,nD_beta[,jjj]], nobs) +accur_beta_temp_nD
        }
        if(jjj ==2){
          accuracy_list[[jj]]$accur_beta_2D <-accur_beta_temp_nD/ncol(nD_beta)
        }else if(jjj==3){
          accuracy_list[[jj]]$accur_beta_3D <-accur_beta_temp_nD/ncol(nD_beta)
        }
      }
      ## nD lmat
      for(jjj in 2:5){
        nD_lmat <- combn(1:6, jjj)
        accur_lmat_temp_nD <-0
        
        for(jjjj in ncol(nD_lmat)){
          accur_lmat_temp_nD <- Accuracy(res[[jj]]$lmatNew$wasp[,nD_lmat[,jjj]],resfull[[jj]]$vechlmat[,nD_lmat[,jjj]], nobs) +accur_lmat_temp_nD
        }
        if(jjj ==2){
          accuracy_list[[jj]]$accur_lmat_2D <-accur_lmat_temp_nD/ncol(nD_lmat)
        }else if(jjj==3){
          accuracy_list[[jj]]$accur_lmat_3D <-accur_lmat_temp_nD/ncol(nD_lmat)
        }else if(jjj==4){
          accuracy_list[[jj]]$accur_lmat_4D <-accur_lmat_temp_nD/ncol(nD_lmat)
        }else if(jjj==5){
          accuracy_list[[jj]]$accur_lmat_5D <-accur_lmat_temp_nD/ncol(nD_lmat)
        }
      }
      ## nD dmat
      for(jjj in 2:5){
        nD_dmat <- combn(1:6, jjj)
        accur_dmat_temp_nD <-0
        
        for(jjjj in ncol(nD_dmat)){
          accur_dmat_temp_nD <- Accuracy(res[[jj]]$dmatNew$vecDmat[,nD_dmat[,jjj]],resfull[[jj]]$vecdmat[,nD_dmat[,jjj]], nobs) +accur_dmat_temp_nD
        }
        if(jjj ==2){
          accuracy_list[[jj]]$accur_dmat_2D <-accur_dmat_temp_nD/ncol(nD_dmat)
        }else if(jjj==3){
          accuracy_list[[jj]]$accur_dmat_3D <-accur_dmat_temp_nD/ncol(nD_dmat)
        }else if(jjj==4){
          accuracy_list[[jj]]$accur_dmat_4D <-accur_dmat_temp_nD/ncol(nD_dmat)
        }else if(jjj==5){
          accuracy_list[[jj]]$accur_dmat_5D <-accur_dmat_temp_nD/ncol(nD_dmat)
        }
      }
      print(jj)
    }      
    accuracy_matrix <- matrix(unlist(accuracy_list), ncol = 14,byrow=TRUE)
    colnames(accuracy_matrix) <- c("accur_beta", "accur_lmat", "accur_dmat", "accur_Rmat", "accur_beta_2D",
                                   "accur_beta_3D","accur_lmat_2D","accur_lmat_3D","accur_lmat_4D","accur_lmat_5D",
                                   "accur_dmat_2D","accur_dmat_3D","accur_dmat_4D","accur_dmat_5D")
    accuracy_mean <- apply(accuracy_matrix, 2, mean)
    accuracy_sd <- apply(accuracy_matrix, 2,sd)
    result <- vector("list",1)
    result[[1]]$accuracy_matrix <-accuracy_matrix
    result[[1]]$accuracy_mean <-accuracy_mean
    result[[1]]$accuracy_sd <-accuracy_sd
    
    fname <- paste0("../result/real_dls_k_", subset, "_accuracy_samp500.rds")
    saveRDS(result, fname)
    
    
  }else if (mtd == 14) {
    # XL

    source("ganbo_accuracy.R")
    
    nsubset <- c(20, 30)
    subset <- nsubset[id]
    nrep <- 1:10
 
    #########  to get accuracy measure for disjoint partition #############
       
    resfull <- vector("list",10)
    
    for(ii in nrep){
      data_name <- paste0("../result/full/real_full_rep_",ii,".rds")
      resfull[[ii]] <- readRDS(data_name)
    }
    
    res <- vector("list",10)
    
    ## we import the same data because xl and ours are stored in the same dataframe 
    for(ii in nrep){
      data_name <- paste0("../result/comb/real_dls_samps_rep_",ii,"_k_", subset, ".rds")
      res[[ii]] <- readRDS(data_name)
    }
    
    nobs <- 20000 
    ## each for one replication in total 10
    accuracy_list <- vector("list",10)
    ## get the accuracy for the whole matrix of betas, lmat, Dmat, Rmat
    for(jj in 1:10){
      accur_beta_temp <- Accuracy(res[[jj]]$betas$xueLiang,resfull[[jj]]$beta, nobs)
      accur_lmat_temp <- Accuracy(res[[jj]]$lmat$xueLiang,resfull[[jj]]$vechlmat, nobs)
      accur_dmat_temp <- Accuracy(res[[jj]]$dmat$xlvecDmat,resfull[[jj]]$vecdmat, nobs)
      accur_Rmat_temp <- Accuracy(res[[jj]]$dmat$xlvecRmat[,4:6],resfull[[jj]]$vecrmat[,4:6], nobs)
      
      accuracy_list[[jj]]$accur_beta <-accur_beta_temp
      accuracy_list[[jj]]$accur_lmat <-accur_lmat_temp
      accuracy_list[[jj]]$accur_dmat <-accur_dmat_temp
      accuracy_list[[jj]]$accur_Rmat <-accur_Rmat_temp
      
      print(jj)
    }  
    
    ## get the accuracy for the selected matrix of betas, lmat, Dmat, Rmat
    for(jj in 1:10){
      ## select 2-dimensions, or 3-dimensions of 4-dimension beta matrix  
      for(jjj in 2:3){
        nD_beta <- combn(1:4, jjj)
        accur_beta_temp_nD <-0
        
        for(jjjj in ncol(nD_beta)){
          accur_beta_temp_nD <- Accuracy(res[[jj]]$betas$xueLiang[,nD_beta[,jjj]],resfull[[jj]]$beta[,nD_beta[,jjj]], nobs) +accur_beta_temp_nD
        }
        if(jjj ==2){
          accuracy_list[[jj]]$accur_beta_2D <-accur_beta_temp_nD/ncol(nD_beta)
        }else if(jjj==3){
          accuracy_list[[jj]]$accur_beta_3D <-accur_beta_temp_nD/ncol(nD_beta)
        }
      }
      ## select 2,3,4,5-dimensions of 6-dimension L matrix 
      for(jjj in 2:5){
        nD_lmat <- combn(1:6, jjj)
        accur_lmat_temp_nD <-0
        
        for(jjjj in ncol(nD_lmat)){
          accur_lmat_temp_nD <- Accuracy(res[[jj]]$lmat$xueLiang[,nD_lmat[,jjj]],resfull[[jj]]$vechlmat[,nD_lmat[,jjj]], nobs) +accur_lmat_temp_nD
        }
        if(jjj ==2){
          accuracy_list[[jj]]$accur_lmat_2D <-accur_lmat_temp_nD/ncol(nD_lmat)
        }else if(jjj==3){
          accuracy_list[[jj]]$accur_lmat_3D <-accur_lmat_temp_nD/ncol(nD_lmat)
        }else if(jjj==4){
          accuracy_list[[jj]]$accur_lmat_4D <-accur_lmat_temp_nD/ncol(nD_lmat)
        }else if(jjj==5){
          accuracy_list[[jj]]$accur_lmat_5D <-accur_lmat_temp_nD/ncol(nD_lmat)
        }
      }
      ## select 2,3,4,5-dimensions of 6-dimension D matrix 
      for(jjj in 2:5){
        nD_dmat <- combn(1:6, jjj)
        accur_dmat_temp_nD <-0
        
        for(jjjj in ncol(nD_dmat)){
          accur_dmat_temp_nD <- Accuracy(res[[jj]]$dmat$xlvecDmat[,nD_dmat[,jjj]],resfull[[jj]]$vecdmat[,nD_dmat[,jjj]], nobs) +accur_dmat_temp_nD
        }
        if(jjj ==2){
          accuracy_list[[jj]]$accur_dmat_2D <-accur_dmat_temp_nD/ncol(nD_dmat)
        }else if(jjj==3){
          accuracy_list[[jj]]$accur_dmat_3D <-accur_dmat_temp_nD/ncol(nD_dmat)
        }else if(jjj==4){
          accuracy_list[[jj]]$accur_dmat_4D <-accur_dmat_temp_nD/ncol(nD_dmat)
        }else if(jjj==5){
          accuracy_list[[jj]]$accur_dmat_5D <-accur_dmat_temp_nD/ncol(nD_dmat)
        }
      }
      print(jj)
    }      
    
    ## this will work because each var of the sub-list in accuracy_list only contains one item
    ## change into matrix by row will go through each var from 1st sub-list to the 10th sub-list
    ## forming ten rows 
    accuracy_matrix <- matrix(unlist(accuracy_list), ncol = 14,byrow=TRUE)
    colnames(accuracy_matrix) <- c("accur_beta", "accur_lmat", "accur_dmat", "accur_Rmat", "accur_beta_2D",
                                   "accur_beta_3D","accur_lmat_2D","accur_lmat_3D","accur_lmat_4D","accur_lmat_5D",
                                   "accur_dmat_2D","accur_dmat_3D","accur_dmat_4D","accur_dmat_5D")
    accuracy_mean <- apply(accuracy_matrix, 2, mean)
    accuracy_sd <- apply(accuracy_matrix, 2,sd)
    result <- vector("list",1)
    result[[1]]$accuracy_matrix <-accuracy_matrix
    result[[1]]$accuracy_mean <-accuracy_mean
    result[[1]]$accuracy_sd <-accuracy_sd
    
    fname <- paste0("../result//real_XL_dls_k_", subset, "_accuracy_ind.rds")
    saveRDS(result, fname)
    

    #########  to get accuracy measure for overlap 250 #############
    
    
    resfull <- vector("list",10)
    
    for(ii in nrep){
      data_name <- paste0("../result/full/real_full_rep_",ii,".rds")
      resfull[[ii]] <- readRDS(data_name)
    }
    
    res <- vector("list",10)
    
    for(ii in nrep){
      data_name <- paste0("../result/comb/real_dls_samps_rep_",ii,"_k_", subset, "_samp250.rds")
      res[[ii]] <- readRDS(data_name)
    }
    
    nobs <- 20000 
    
    accuracy_list <- vector("list",10)
    
    for(jj in 1:10){
      accur_beta_temp <- Accuracy(res[[jj]]$betas$xueLiang,resfull[[jj]]$beta, nobs)
      accur_lmat_temp <- Accuracy(res[[jj]]$lmat$xueLiang,resfull[[jj]]$vechlmat, nobs)
      accur_dmat_temp <- Accuracy(res[[jj]]$dmat$xlvecDmat,resfull[[jj]]$vecdmat, nobs)
      accur_Rmat_temp <- Accuracy(res[[jj]]$dmat$xlvecRmat[,4:6],resfull[[jj]]$vecrmat[,4:6], nobs)
      
      accuracy_list[[jj]]$accur_beta <-accur_beta_temp
      accuracy_list[[jj]]$accur_lmat <-accur_lmat_temp
      accuracy_list[[jj]]$accur_dmat <-accur_dmat_temp
      accuracy_list[[jj]]$accur_Rmat <-accur_Rmat_temp
      
      print(jj)
    }  
    
    for(jj in 1:10){
      ## nD beta 
      for(jjj in 2:3){
        nD_beta <- combn(1:4, jjj)
        accur_beta_temp_nD <-0
        
        for(jjjj in ncol(nD_beta)){
          accur_beta_temp_nD <- Accuracy(res[[jj]]$betas$xueLiang[,nD_beta[,jjj]],resfull[[jj]]$beta[,nD_beta[,jjj]], nobs) +accur_beta_temp_nD
        }
        if(jjj ==2){
          accuracy_list[[jj]]$accur_beta_2D <-accur_beta_temp_nD/ncol(nD_beta)
        }else if(jjj==3){
          accuracy_list[[jj]]$accur_beta_3D <-accur_beta_temp_nD/ncol(nD_beta)
        }
      }
      ## nD lmat
      for(jjj in 2:5){
        nD_lmat <- combn(1:6, jjj)
        accur_lmat_temp_nD <-0
        
        for(jjjj in ncol(nD_lmat)){
          accur_lmat_temp_nD <- Accuracy(res[[jj]]$lmat$xueLiang[,nD_lmat[,jjj]],resfull[[jj]]$vechlmat[,nD_lmat[,jjj]], nobs) +accur_lmat_temp_nD
        }
        if(jjj ==2){
          accuracy_list[[jj]]$accur_lmat_2D <-accur_lmat_temp_nD/ncol(nD_lmat)
        }else if(jjj==3){
          accuracy_list[[jj]]$accur_lmat_3D <-accur_lmat_temp_nD/ncol(nD_lmat)
        }else if(jjj==4){
          accuracy_list[[jj]]$accur_lmat_4D <-accur_lmat_temp_nD/ncol(nD_lmat)
        }else if(jjj==5){
          accuracy_list[[jj]]$accur_lmat_5D <-accur_lmat_temp_nD/ncol(nD_lmat)
        }
      }
      ## nD dmat
      for(jjj in 2:5){
        nD_dmat <- combn(1:6, jjj)
        accur_dmat_temp_nD <-0
        
        for(jjjj in ncol(nD_dmat)){
          accur_dmat_temp_nD <- Accuracy(res[[jj]]$dmat$xlvecDmat[,nD_dmat[,jjj]],resfull[[jj]]$vecdmat[,nD_dmat[,jjj]], nobs) +accur_dmat_temp_nD
        }
        if(jjj ==2){
          accuracy_list[[jj]]$accur_dmat_2D <-accur_dmat_temp_nD/ncol(nD_dmat)
        }else if(jjj==3){
          accuracy_list[[jj]]$accur_dmat_3D <-accur_dmat_temp_nD/ncol(nD_dmat)
        }else if(jjj==4){
          accuracy_list[[jj]]$accur_dmat_4D <-accur_dmat_temp_nD/ncol(nD_dmat)
        }else if(jjj==5){
          accuracy_list[[jj]]$accur_dmat_5D <-accur_dmat_temp_nD/ncol(nD_dmat)
        }
      }
      print(jj)
    }      
    accuracy_matrix <- matrix(unlist(accuracy_list), ncol = 14,byrow=TRUE)
    colnames(accuracy_matrix) <- c("accur_beta", "accur_lmat", "accur_dmat", "accur_Rmat", "accur_beta_2D",
                                   "accur_beta_3D","accur_lmat_2D","accur_lmat_3D","accur_lmat_4D","accur_lmat_5D",
                                   "accur_dmat_2D","accur_dmat_3D","accur_dmat_4D","accur_dmat_5D")
    accuracy_mean <- apply(accuracy_matrix, 2, mean)
    accuracy_sd <- apply(accuracy_matrix, 2,sd)
    result <- vector("list",1)
    result[[1]]$accuracy_matrix <-accuracy_matrix
    result[[1]]$accuracy_mean <-accuracy_mean
    result[[1]]$accuracy_sd <-accuracy_sd
    
    fname <- paste0("../result/real_XL_dls_k_", subset, "_accuracy_samp250.rds")
    saveRDS(result, fname)
    
    #########  to get accuracy measure for overlap 500 #############
    
    resfull <- vector("list",10)
    
    for(ii in nrep){
      data_name <- paste0("../result/full/real_full_rep_",ii,".rds")
      resfull[[ii]] <- readRDS(data_name)
    }
    
    res <- vector("list",10)
    
    for(ii in nrep){
      data_name <- paste0("../result/comb/real_dls_samps_rep_",ii,"_k_", subset, "_samp500.rds")
      res[[ii]] <- readRDS(data_name)
    }
    
    nobs <- 20000 
    
    accuracy_list <- vector("list",10)
    
    for(jj in 1:10){
      accur_beta_temp <- Accuracy(res[[jj]]$betas$xueLiang,resfull[[jj]]$beta, nobs)
      accur_lmat_temp <- Accuracy(res[[jj]]$lmat$xueLiang,resfull[[jj]]$vechlmat, nobs)
      accur_dmat_temp <- Accuracy(res[[jj]]$dmat$xlvecDmat,resfull[[jj]]$vecdmat, nobs)
      accur_Rmat_temp <- Accuracy(res[[jj]]$dmat$xlvecRmat[,4:6],resfull[[jj]]$vecrmat[,4:6], nobs)
      
      accuracy_list[[jj]]$accur_beta <-accur_beta_temp
      accuracy_list[[jj]]$accur_lmat <-accur_lmat_temp
      accuracy_list[[jj]]$accur_dmat <-accur_dmat_temp
      accuracy_list[[jj]]$accur_Rmat <-accur_Rmat_temp
      
      print(jj)
    }  
    
    for(jj in 1:10){
      ## nD beta 
      for(jjj in 2:3){
        nD_beta <- combn(1:4, jjj)
        accur_beta_temp_nD <-0
        
        for(jjjj in ncol(nD_beta)){
          accur_beta_temp_nD <- Accuracy(res[[jj]]$betas$xueLiang[,nD_beta[,jjj]],resfull[[jj]]$beta[,nD_beta[,jjj]], nobs) +accur_beta_temp_nD
        }
        if(jjj ==2){
          accuracy_list[[jj]]$accur_beta_2D <-accur_beta_temp_nD/ncol(nD_beta)
        }else if(jjj==3){
          accuracy_list[[jj]]$accur_beta_3D <-accur_beta_temp_nD/ncol(nD_beta)
        }
      }
      ## nD lmat
      for(jjj in 2:5){
        nD_lmat <- combn(1:6, jjj)
        accur_lmat_temp_nD <-0
        
        for(jjjj in ncol(nD_lmat)){
          accur_lmat_temp_nD <- Accuracy(res[[jj]]$lmat$xueLiang[,nD_lmat[,jjj]],resfull[[jj]]$vechlmat[,nD_lmat[,jjj]], nobs) +accur_lmat_temp_nD
        }
        if(jjj ==2){
          accuracy_list[[jj]]$accur_lmat_2D <-accur_lmat_temp_nD/ncol(nD_lmat)
        }else if(jjj==3){
          accuracy_list[[jj]]$accur_lmat_3D <-accur_lmat_temp_nD/ncol(nD_lmat)
        }else if(jjj==4){
          accuracy_list[[jj]]$accur_lmat_4D <-accur_lmat_temp_nD/ncol(nD_lmat)
        }else if(jjj==5){
          accuracy_list[[jj]]$accur_lmat_5D <-accur_lmat_temp_nD/ncol(nD_lmat)
        }
      }
      ## nD dmat
      for(jjj in 2:5){
        nD_dmat <- combn(1:6, jjj)
        accur_dmat_temp_nD <-0
        
        for(jjjj in ncol(nD_dmat)){
          accur_dmat_temp_nD <- Accuracy(res[[jj]]$dmat$xlvecDmat[,nD_dmat[,jjj]],resfull[[jj]]$vecdmat[,nD_dmat[,jjj]], nobs) +accur_dmat_temp_nD
        }
        if(jjj ==2){
          accuracy_list[[jj]]$accur_dmat_2D <-accur_dmat_temp_nD/ncol(nD_dmat)
        }else if(jjj==3){
          accuracy_list[[jj]]$accur_dmat_3D <-accur_dmat_temp_nD/ncol(nD_dmat)
        }else if(jjj==4){
          accuracy_list[[jj]]$accur_dmat_4D <-accur_dmat_temp_nD/ncol(nD_dmat)
        }else if(jjj==5){
          accuracy_list[[jj]]$accur_dmat_5D <-accur_dmat_temp_nD/ncol(nD_dmat)
        }
      }
      print(jj)
    }      
    accuracy_matrix <- matrix(unlist(accuracy_list), ncol = 14,byrow=TRUE)
    colnames(accuracy_matrix) <- c("accur_beta", "accur_lmat", "accur_dmat", "accur_Rmat", "accur_beta_2D",
                                   "accur_beta_3D","accur_lmat_2D","accur_lmat_3D","accur_lmat_4D","accur_lmat_5D",
                                   "accur_dmat_2D","accur_dmat_3D","accur_dmat_4D","accur_dmat_5D")
    accuracy_mean <- apply(accuracy_matrix, 2, mean)
    accuracy_sd <- apply(accuracy_matrix, 2,sd)
    result <- vector("list",1)
    result[[1]]$accuracy_matrix <-accuracy_matrix
    result[[1]]$accuracy_mean <-accuracy_mean
    result[[1]]$accuracy_sd <-accuracy_sd
    
    fname <- paste0("../result/real_XL_dls_k_", subset, "_accuracy_samp500.rds")
    saveRDS(result, fname)
    
    
  }else if (mtd == 15) {
    
    source("ganbo_accuracy.R")
    
    nsubset <- c(20, 30)
    subset <- nsubset[id]
    nrep <- 1:10
    
    ####################### to get the accuracy of individual entries for Dmat and Rmat under our method ########
    
    resfull <- vector("list",10)
    
    for(ii in nrep){
      data_name <- paste0("../result/full/real_full_rep_",ii,".rds")
      resfull[[ii]] <- readRDS(data_name)
    }
    
    res <- vector("list",10)
    
    for(ii in nrep){
      data_name <- paste0("../result/comb/real_dls_samps_rep_",ii,"_k_", subset, ".rds")
      res[[ii]] <- readRDS(data_name)
    }
    
    nobs <- 20000 
    
    accuracy_list <- vector("list",10)
    
    for(j in 1:10){
      
      accuracy_list[[j]]$accur_dmat1<- Accuracy(as.matrix(res[[j]]$dmat$vecDmat[,1]),as.matrix(resfull[[j]]$vecdmat[,1]), nobs)
      accuracy_list[[j]]$accur_dmat2<- Accuracy(as.matrix(res[[j]]$dmat$vecDmat[,2]),as.matrix(resfull[[j]]$vecdmat[,2]), nobs)
      accuracy_list[[j]]$accur_dmat3<- Accuracy(as.matrix(res[[j]]$dmat$vecDmat[,3]),as.matrix(resfull[[j]]$vecdmat[,3]), nobs)
      accuracy_list[[j]]$accur_dmat4<- Accuracy(as.matrix(res[[j]]$dmat$vecDmat[,4]),as.matrix(resfull[[j]]$vecdmat[,4]), nobs)
      accuracy_list[[j]]$accur_dmat5<- Accuracy(as.matrix(res[[j]]$dmat$vecDmat[,5]),as.matrix(resfull[[j]]$vecdmat[,5]), nobs)
      accuracy_list[[j]]$accur_dmat6<- Accuracy(as.matrix(res[[j]]$dmat$vecDmat[,6]),as.matrix(resfull[[j]]$vecdmat[,6]), nobs)
      
      accuracy_list[[j]]$accur_Rmat1<- Accuracy(as.matrix(res[[j]]$dmat$vecRmat[,1]),as.matrix(resfull[[j]]$vecrmat[,1]), nobs)
      accuracy_list[[j]]$accur_Rmat2<- Accuracy(as.matrix(res[[j]]$dmat$vecRmat[,2]),as.matrix(resfull[[j]]$vecrmat[,2]), nobs)
      accuracy_list[[j]]$accur_Rmat3<- Accuracy(as.matrix(res[[j]]$dmat$vecRmat[,3]),as.matrix(resfull[[j]]$vecrmat[,3]), nobs)
      accuracy_list[[j]]$accur_Rmat4<- Accuracy(as.matrix(res[[j]]$dmat$vecRmat[,4]),as.matrix(resfull[[j]]$vecrmat[,4]), nobs)
      accuracy_list[[j]]$accur_Rmat5<- Accuracy(as.matrix(res[[j]]$dmat$vecRmat[,5]),as.matrix(resfull[[j]]$vecrmat[,5]), nobs)
      accuracy_list[[j]]$accur_Rmat6<- Accuracy(as.matrix(res[[j]]$dmat$vecRmat[,6]),as.matrix(resfull[[j]]$vecrmat[,6]), nobs)
      

      print(j)
    }   
    
    ## this will work because each var of the sub-list in accuracy_list only contains one item
    ## change into matrix by row will go through each var from 1st sub-list to the 10th sub-list
    ## forming ten rows 
    accuracy_matrix <- matrix(unlist(accuracy_list), ncol = 12,byrow=TRUE)
    colnames(accuracy_matrix) <- c("accur_dmat1","accur_dmat2","accur_dmat3","accur_dmat4","accur_dmat5",
                                   "accur_dmat6","accur_Rmat1","accur_Rmat2","accur_Rmat3","accur_Rmat4",
                                   "accur_Rmat5","accur_Rmat6")
    accuracy_mean <- apply(accuracy_matrix, 2, mean)
    accuracy_sd <- apply(accuracy_matrix, 2,sd)
    result <- vector("list",1)
    result[[1]]$accuracy_matrix <-accuracy_matrix
    result[[1]]$accuracy_mean <-accuracy_mean
    result[[1]]$accuracy_sd <-accuracy_sd
    
    fname <- paste0("../result/real_dls_k_", subset, "_accuracy_ind_covcorr.rds")
    saveRDS(result, fname)
    
    ####################### to get the accuracy of individual entries for Dmat and Rmat under xl method ########
    resfull <- vector("list",10)
    
    for(ii in nrep){
      data_name <- paste0("../result/full/real_full_rep_",ii,".rds")
      resfull[[ii]] <- readRDS(data_name)
    }
    
    res <- vector("list",10)
    
    for(ii in nrep){
      data_name <- paste0("../result/comb/real_dls_samps_rep_",ii,"_k_", subset, ".rds")
      res[[ii]] <- readRDS(data_name)
    }
    
    nobs <- 20000 
    
    accuracy_list <- vector("list",10)
    
    for(j in 1:10){
      
      accuracy_list[[j]]$accur_xldmat1<- Accuracy(as.matrix(res[[j]]$dmat$xlvecDmat[,1]),as.matrix(resfull[[j]]$vecdmat[,1]), nobs)
      accuracy_list[[j]]$accur_xldmat2<- Accuracy(as.matrix(res[[j]]$dmat$xlvecDmat[,2]),as.matrix(resfull[[j]]$vecdmat[,2]), nobs)
      accuracy_list[[j]]$accur_xldmat3<- Accuracy(as.matrix(res[[j]]$dmat$xlvecDmat[,3]),as.matrix(resfull[[j]]$vecdmat[,3]), nobs)
      accuracy_list[[j]]$accur_xldmat4<- Accuracy(as.matrix(res[[j]]$dmat$xlvecDmat[,4]),as.matrix(resfull[[j]]$vecdmat[,4]), nobs)
      accuracy_list[[j]]$accur_xldmat5<- Accuracy(as.matrix(res[[j]]$dmat$xlvecDmat[,5]),as.matrix(resfull[[j]]$vecdmat[,5]), nobs)
      accuracy_list[[j]]$accur_xldmat6<- Accuracy(as.matrix(res[[j]]$dmat$xlvecDmat[,6]),as.matrix(resfull[[j]]$vecdmat[,6]), nobs)
      
      accuracy_list[[j]]$accur_xlRmat1<- Accuracy(as.matrix(res[[j]]$dmat$xlvecRmat[,1]),as.matrix(resfull[[j]]$vecrmat[,1]), nobs)
      accuracy_list[[j]]$accur_xlRmat2<- Accuracy(as.matrix(res[[j]]$dmat$xlvecRmat[,2]),as.matrix(resfull[[j]]$vecrmat[,2]), nobs)
      accuracy_list[[j]]$accur_xlRmat3<- Accuracy(as.matrix(res[[j]]$dmat$xlvecRmat[,3]),as.matrix(resfull[[j]]$vecrmat[,3]), nobs)
      accuracy_list[[j]]$accur_xlRmat4<- Accuracy(as.matrix(res[[j]]$dmat$xlvecRmat[,4]),as.matrix(resfull[[j]]$vecrmat[,4]), nobs)
      accuracy_list[[j]]$accur_xlRmat5<- Accuracy(as.matrix(res[[j]]$dmat$xlvecRmat[,5]),as.matrix(resfull[[j]]$vecrmat[,5]), nobs)
      accuracy_list[[j]]$accur_xlRmat6<- Accuracy(as.matrix(res[[j]]$dmat$xlvecRmat[,6]),as.matrix(resfull[[j]]$vecrmat[,6]), nobs)
      
      
      print(j)
    }      
    accuracy_matrix <- matrix(unlist(accuracy_list), ncol = 12,byrow=TRUE)
    colnames(accuracy_matrix) <- c("accur_dmat1","accur_dmat2","accur_dmat3","accur_dmat4","accur_dmat5",
                                   "accur_dmat6","accur_Rmat1","accur_Rmat2","accur_Rmat3","accur_Rmat4",
                                   "accur_Rmat5","accur_Rmat6")
    accuracy_mean <- apply(accuracy_matrix, 2, mean)
    accuracy_sd <- apply(accuracy_matrix, 2,sd)
    result <- vector("list",1)
    result[[1]]$accuracy_matrix <-accuracy_matrix
    result[[1]]$accuracy_mean <-accuracy_mean
    result[[1]]$accuracy_sd <-accuracy_sd
    
    fname <- paste0("../result/real_XL_dls_k_", subset, "_accuracy_ind_covcorr.rds")
    saveRDS(result, fname)
    
     

  }else {
    print("peace")
}
