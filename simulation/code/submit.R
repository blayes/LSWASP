
cmdArgs <- commandArgs(trailingOnly = TRUE)

mtd <- as.numeric(cmdArgs[1])
id <- as.numeric(cmdArgs[2])

if (mtd == 1) {
    source("dls_sampler.R")

    cvtrain <- readRDS("../data/full_data.rds")
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

    nranef0 <- 3; nfixef0 <- 4;
    muBeta0 <- rep(0, nfixef0); sigBetaInv0 <- diag(0, nfixef0); nu0 <- (2 + nfixef0); sig0 <- 0;
    eta0 <- -(nranef0 + 1); tmat0 <- diag(0, nranef0);
    dmat0 <- diag(1, nranef0);
    errVar0 <- 10;
    sigLInv0 <- diag(0.01, nranef0 * (nranef0 + 1) / 2);
    muL0 <- rep(0.0, nranef0 * (nranef0 + 1) / 2)

    res <- waspLmeSampler(ylist0, fixefList0, ranefList0, rmatList0, 1,
                          niter = 10000, nburn = 5000, nthin = 5,
                          dmat0, errVar0, muBeta0, sigBetaInv0, muL0, sigLInv0, sig0, nu0)

    fname <- paste0("~/simulation/result/sim_full_rep_", id, ".rds")
    saveRDS(res, fname)
} else if (mtd == 2) {
    source("dls_sampler.R")

    cvs <- rep(1:10, each = 20)
    nsubs <- rep(1:20, times = 10)

    cid <- cvs[id]
    sid <- nsubs[id]

    fname <- paste0("../data/part_data_rep_", cid, "_k_20.rds")
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

    nranef0 <- 3; nfixef0 <- 4;
    muBeta0 <- rep(0, nfixef0); sigBetaInv0 <- diag(0, nfixef0); nu0 <- (2 + nfixef0); sig0 <- 0;
    eta0 <- -(nranef0 + 1); tmat0 <- diag(0, nranef0);
    dmat0 <- diag(1, nranef0);
    errVar0 <- 10;
    sigLInv0 <- diag(0.01, nranef0 * (nranef0 + 1) / 2);
    muL0 <- rep(0.0, nranef0 * (nranef0 + 1) / 2)

    res <- waspLmeSampler(ylist0, fixefList0, ranefList0, rmatList0, train$nrep,
                          niter = 10000, nburn = 5000, nthin = 5,
                          dmat0, errVar0, muBeta0, sigBetaInv0, muL0, sigLInv0, sig0, nu0)

    fname <- paste0("~/simulation/result/samp/sub20/ind_dls_rep_", cid, "_nsub_", sid, "_k_20.rds")
    saveRDS(res, fname)
} else if (mtd == 3) {
    source("dls_sampler.R")

    cvs <- rep(1:10, each = 50)
    nsubs <- rep(1:50, times = 10)

    cid <- cvs[id]
    sid <- nsubs[id]

    fname <- paste0("../data/part_data_rep_", cid, "_k_50.rds")
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

    nranef0 <- 3; nfixef0 <- 4;
    muBeta0 <- rep(0, nfixef0); sigBetaInv0 <- diag(0, nfixef0); nu0 <- (2 + nfixef0); sig0 <- 0;
    eta0 <- -(nranef0 + 1); tmat0 <- diag(0, nranef0);
    dmat0 <- diag(1, nranef0);
    errVar0 <- 10;
    sigLInv0 <- diag(0.01, nranef0 * (nranef0 + 1) / 2);
    muL0 <- rep(0.0, nranef0 * (nranef0 + 1) / 2)

    res <- waspLmeSampler(ylist0, fixefList0, ranefList0, rmatList0, train$nrep,
                          niter = 10000, nburn = 5000, nthin = 5,
                          dmat0, errVar0, muBeta0, sigBetaInv0, muL0, sigLInv0, sig0, nu0)

    fname <- paste0("~/simulation/result/samp/sub50/ind_dls_rep_", cid, "_nsub_", sid, "_k_50.rds")
    saveRDS(res, fname)
} else if (mtd == 4) {
    source("dls_sampler.R")

    cvs <- rep(1:10, each = 100)
    nsubs <- rep(1:100, times = 10)

    cid <- cvs[id]
    sid <- nsubs[id]

    fname <- paste0("../data/part_data_rep_", cid, "_k_100.rds")
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

    nranef0 <- 3; nfixef0 <- 4;
    muBeta0 <- rep(0, nfixef0); sigBetaInv0 <- diag(0, nfixef0); nu0 <- (2 + nfixef0); sig0 <- 0;
    eta0 <- -(nranef0 + 1); tmat0 <- diag(0, nranef0);
    dmat0 <- diag(1, nranef0);
    errVar0 <- 10;
    sigLInv0 <- diag(0.01, nranef0 * (nranef0 + 1) / 2);
    muL0 <- rep(0.0, nranef0 * (nranef0 + 1) / 2)

    res <- waspLmeSampler(ylist0, fixefList0, ranefList0, rmatList0, train$nrep,
                          niter = 10000, nburn = 5000, nthin = 5,
                          dmat0, errVar0, muBeta0, sigBetaInv0, muL0, sigLInv0, sig0, nu0)

    fname <- paste0("~/simulation/result/samp/sub100/ind_dls_rep_", cid, "_nsub_", sid, "_k_100.rds")
    saveRDS(res, fname)
} else if (mtd == 5) {
    source("dls_sampler.R")

    cvs <- rep(1:10, each = 150)
    nsubs <- rep(1:150, times = 10)

    cid <- cvs[id]
    sid <- nsubs[id]

    fname <- paste0("../data/part_data_rep_", cid, "_k_150.rds")
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

    nranef0 <- 3; nfixef0 <- 4;
    muBeta0 <- rep(0, nfixef0); sigBetaInv0 <- diag(0, nfixef0); nu0 <- (2 + nfixef0); sig0 <- 0;
    eta0 <- -(nranef0 + 1); tmat0 <- diag(0, nranef0);
    dmat0 <- diag(1, nranef0);
    errVar0 <- 10;
    sigLInv0 <- diag(0.01, nranef0 * (nranef0 + 1) / 2);
    muL0 <- rep(0.0, nranef0 * (nranef0 + 1) / 2)

    res <- waspLmeSampler(ylist0, fixefList0, ranefList0, rmatList0, train$nrep,
                          niter = 10000, nburn = 5000, nthin = 5,
                          dmat0, errVar0, muBeta0, sigBetaInv0, muL0, sigLInv0, sig0, nu0)

    fname <- paste0("~/simulation/result//samp/sub150/ind_dls_rep_", cid, "_nsub_", sid, "_k_150.rds")
    saveRDS(res, fname)
} else if (mtd == 6) {
    source("dls_sampler.R")

    cvs <- rep(1:10, each = 200)
    nsubs <- rep(1:200, times = 10)

    cid <- cvs[id]
    sid <- nsubs[id]

    fname <- paste0("../data/part_data_rep_", cid, "_k_200.rds")
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

    nranef0 <- 3; nfixef0 <- 4;
    muBeta0 <- rep(0, nfixef0); sigBetaInv0 <- diag(0, nfixef0); nu0 <- (2 + nfixef0); sig0 <- 0;
    eta0 <- -(nranef0 + 1); tmat0 <- diag(0, nranef0);
    dmat0 <- diag(1, nranef0);
    errVar0 <- 10;
    sigLInv0 <- diag(0.01, nranef0 * (nranef0 + 1) / 2);
    muL0 <- rep(0.0, nranef0 * (nranef0 + 1) / 2)

    res <- waspLmeSampler(ylist0, fixefList0, ranefList0, rmatList0, train$nrep,
                          niter = 10000, nburn = 5000, nthin = 5,
                          dmat0, errVar0, muBeta0, sigBetaInv0, muL0, sigLInv0, sig0, nu0)

    fname <- paste0("~/simulation/result/samp/sub200/ind_dls_rep_", cid, "_nsub_", sid, "_k_200.rds")
    saveRDS(res, fname)
} else if (mtd == 7) {
    source("dls_comb.R")

    cvs <- 1:10
    cid <- cvs[id]

    for (nsub in c(20, 50, 100, 150, 200)) {
        res <- list()
        for (sid in 1:nsub) {
            foldname <- paste0("sub",nsub)
            fname <- paste0("~/simulation/result/samp/",foldname,"/ind_dls_rep_", cid, "_nsub_", sid, "_k_", nsub, ".rds")
            res[[sid]] <- readRDS(fname)
        }

        lList <- lapply(res, function(x) x$lmat)
        bList <- lapply(res, function(x) x$beta)

        stime <- proc.time()
        lres <- sampleLmat(lList)
        bres <- sampleBetas(bList)
        dcres <- sampleDmatCorrMat(lList, ndim = 3)
        etime <- proc.time()

        rname <- paste0("~/simulation/result/comb/ind_dls_samps_rep_", cid, "_k_", nsub, ".rds")
        fres <- list(lmat = lres, betas = bres, dmat = dcres,
                     time =  mean(unlist(lapply(res, function(x) x$time))) + etime[3] - stime[3]
                     )
        saveRDS(fres, rname)
    }else if(mtd == 8){
  
  source("ganbo_accuracy.R")
  
  nsubset <- c(20, 50, 100, 150, 200)
  subset <- nsubset[id]
  nrep <- 1:10
  
  #########  to get accuracy measure for disjoint #############
  #########  DLS
  
  resfull <- vector("list",10)
  
  for(ii in nrep){
    data_name <- paste0("../result/full/sim_full_rep_",ii,".rds")
    resfull[[ii]] <- readRDS(data_name)
  }
  
  res <- vector("list",10)
  
  for(ii in nrep){
    data_name <- paste0("../result/comb/ind_dls_samps_rep_",ii,"_k_", subset, ".rds")
    res[[ii]] <- readRDS(data_name)
  }
  # nobs comes from length(unique(repData[j]$group)) # repData <- readRDS("train_bwt.rds")
  nobs <- 50000
  ## each for one replication in total 10
  accuracy_list <- vector("list",10)

  ## get the accuracy for the selected matrix of betas, lmat, Dmat, Rmat
  for(jj in 1:10){
    
    ## select 2-6 dimension beta matrix
    for(jjj in 2:4){
      nD_beta <- combn(1:4, jjj)
      accur_beta_temp_nD <-0
      for(jjjj in ncol(nD_beta)){
        accur_beta_temp_nD <- Accuracy(res[[jj]]$betasNew$samps[,nD_beta[,jjjj]],resfull[[jj]]$beta[,nD_beta[,jjjj]], nobs) +accur_beta_temp_nD
      }
      if(jjj ==2){
        accuracy_list[[jj]]$accur_beta_2D <-accur_beta_temp_nD/ncol(nD_beta)
      }else if(jjj==3){
        accuracy_list[[jj]]$accur_beta_3D <-accur_beta_temp_nD/ncol(nD_beta)
      }else if(jjj==4){
        accuracy_list[[jj]]$accur_beta_4D <-accur_beta_temp_nD/ncol(nD_beta)
      }
    }
    
    ## select 2-5 dimension D matrix
    for(jjj in 2:6){
      nD_dmat <- combn(1:6, jjj)
      accur_dmat_temp_nD <-0
      for(jjjj in ncol(nD_dmat)){
        accur_dmat_temp_nD <- Accuracy(res[[jj]]$dmatNew$vecDmat[,nD_dmat[,jjjj]],resfull[[jj]]$vecdmat[,nD_dmat[,jjjj]], nobs) +accur_dmat_temp_nD
      }
      if(jjj ==2){
        accuracy_list[[jj]]$accur_dmat_2D <-accur_dmat_temp_nD/ncol(nD_dmat)
      }else if(jjj==3){
        accuracy_list[[jj]]$accur_dmat_3D <-accur_dmat_temp_nD/ncol(nD_dmat)
      }else if(jjj==4){
        accuracy_list[[jj]]$accur_dmat_4D <-accur_dmat_temp_nD/ncol(nD_dmat)
      }else if(jjj==5){
        accuracy_list[[jj]]$accur_dmat_5D <-accur_dmat_temp_nD/ncol(nD_dmat)
      }else if(jjj==6){
        accuracy_list[[jj]]$accur_dmat_6D <-accur_dmat_temp_nD/ncol(nD_dmat)
      }
    }
    
    ## select 2 dimension R matrix
    for(jjj in 2:3){
      nD_Rmat <- combn(4:6, jjj)
      accur_Rmat_temp_nD <-0
      for(jjjj in ncol(nD_Rmat)){
        accur_Rmat_temp_nD <- Accuracy(res[[jj]]$dmatNew$vecRmat[,nD_Rmat[,jjjj]],resfull[[jj]]$vecrmat[,nD_Rmat[,jjjj]], nobs) +accur_Rmat_temp_nD
      }
      if(jjj ==2){
        accuracy_list[[jj]]$accur_Rmat_2D <-accur_Rmat_temp_nD/ncol(nD_Rmat)
      }else if(jjj==3){
        accuracy_list[[jj]]$accur_Rmat_3D <-accur_Rmat_temp_nD/ncol(nD_Rmat)
      }
    }
  }
  ## this will work because each var of the sub-list in accuracy_list only contains one item
  ## change into matrix by row will go through each var from 1st sub-list to the 10th sub-list
  ## forming ten rows
  accuracy_matrix <- matrix(unlist(accuracy_list), ncol = 10,byrow=TRUE)
  colnames(accuracy_matrix) <- c(  "accur_beta_2D","accur_beta_3D","accur_beta_4D",
                                   "accur_dmat_2D","accur_dmat_3D","accur_dmat_4D","accur_dmat_5D","accur_dmat_6D",
                                   "accur_Rmat_2D","accur_Rmat_3D")
  accuracy_mean <- apply(accuracy_matrix, 2, mean)
  accuracy_sd <- apply(accuracy_matrix, 2,sd)
  result <- vector("list",1)
  result[[1]]$accuracy_matrix <-accuracy_matrix
  result[[1]]$accuracy_mean <-accuracy_mean
  result[[1]]$accuracy_sd <-accuracy_sd
  
  fname <- paste0("~/simulation/result/sim_dls_k_", subset, "_accuracy_ind.rds")
  saveRDS(result, fname)

  #########  XL 
  resfull <- vector("list",10)
  
  for(ii in nrep){
    data_name <- paste0("../result/full/sim_full_rep_",ii,".rds")
    resfull[[ii]] <- readRDS(data_name)
  }
  
  res <- vector("list",10)
  
  for(ii in nrep){
    data_name <- paste0("../result/comb/ind_dls_samps_rep_",ii,"_k_", subset, ".rds")
    res[[ii]] <- readRDS(data_name)
  }

  nobs <- 50000
  ## each for one replication in total 10
  accuracy_list <- vector("list",10)

  ## get the accuracy for the selected matrix of betas, lmat, Dmat, Rmat
  for(jj in 1:10){
    
    ## select 2-4 dimension beta matrix
    for(jjj in 2:4){
      nD_beta <- combn(1:4, jjj)
      accur_beta_temp_nD <-0
      for(jjjj in ncol(nD_beta)){
        accur_beta_temp_nD <- Accuracy(res[[jj]]$betas$xueLiang[,nD_beta[,jjjj]],resfull[[jj]]$beta[,nD_beta[,jjjj]], nobs) +accur_beta_temp_nD
      }
      if(jjj ==2){
        accuracy_list[[jj]]$accur_beta_2D <-accur_beta_temp_nD/ncol(nD_beta)
      }else if(jjj==3){
        accuracy_list[[jj]]$accur_beta_3D <-accur_beta_temp_nD/ncol(nD_beta)
      }else if(jjj==4){
        accuracy_list[[jj]]$accur_beta_4D <-accur_beta_temp_nD/ncol(nD_beta)
      }
    }
    ## select 2-6 dimension D matrix
    for(jjj in 2:6){
      nD_dmat <- combn(1:6, jjj)
      accur_dmat_temp_nD <-0
      for(jjjj in ncol(nD_dmat)){
        accur_dmat_temp_nD <- Accuracy(res[[jj]]$dmat$xlvecDmat[,nD_dmat[,jjjj]],resfull[[jj]]$vecdmat[,nD_dmat[,jjjj]], nobs) +accur_dmat_temp_nD
      }
      if(jjj ==2){
        accuracy_list[[jj]]$accur_dmat_2D <-accur_dmat_temp_nD/ncol(nD_dmat)
      }else if(jjj==3){
        accuracy_list[[jj]]$accur_dmat_3D <-accur_dmat_temp_nD/ncol(nD_dmat)
      }else if(jjj==4){
        accuracy_list[[jj]]$accur_dmat_4D <-accur_dmat_temp_nD/ncol(nD_dmat)
      }else if(jjj==5){
        accuracy_list[[jj]]$accur_dmat_5D <-accur_dmat_temp_nD/ncol(nD_dmat)
      }else if(jjj==6){
        accuracy_list[[jj]]$accur_dmat_6D <-accur_dmat_temp_nD/ncol(nD_dmat)
      }
    }
    ## select 2 dimension R matrix
    for(jjj in 2:3){
      nD_Rmat <- combn(4:6, jjj)
      accur_Rmat_temp_nD <-0
      for(jjjj in ncol(nD_Rmat)){
        accur_Rmat_temp_nD <- Accuracy(res[[jj]]$dmat$xlvecRmat[,nD_Rmat[,jjjj]],resfull[[jj]]$vecrmat[,nD_Rmat[,jjjj]], nobs) +accur_Rmat_temp_nD
      }
      if(jjj ==2){
        accuracy_list[[jj]]$accur_Rmat_2D <-accur_Rmat_temp_nD/ncol(nD_Rmat)
      }else if(jjj==3){
        accuracy_list[[jj]]$accur_Rmat_3D <-accur_Rmat_temp_nD/ncol(nD_Rmat)
      }
    }
  }

  accuracy_matrix <- matrix(unlist(accuracy_list), ncol = 10,byrow=TRUE)
  colnames(accuracy_matrix) <- c(  "accur_beta_2D","accur_beta_3D","accur_beta_4D",
                                   "accur_dmat_2D","accur_dmat_3D","accur_dmat_4D","accur_dmat_5D","accur_dmat_6D",
                                   "accur_Rmat_2D","accur_Rmat_3D")
  accuracy_mean <- apply(accuracy_matrix, 2, mean)
  accuracy_sd <- apply(accuracy_matrix, 2,sd)
  result <- vector("list",1)
  result[[1]]$accuracy_matrix <-accuracy_matrix
  result[[1]]$accuracy_mean <-accuracy_mean
  result[[1]]$accuracy_sd <-accuracy_sd
  
  fname <- paste0("~/simulation/result/sim_XL_dls_k_", subset, "_accuracy_ind.rds")
  saveRDS(result, fname)
  

  
} else {
    print("peace")
}
