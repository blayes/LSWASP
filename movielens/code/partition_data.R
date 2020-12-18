rm(list=ls())
setwd("~/lswasp/ml/code")
mlData <- readRDS("../data/ml_full.rds")


set.seed(12345)
group <- as.numeric(mlData$group)
grpLbl <- sort(unique(group))
ngroup <- length(grpLbl)

usrs <- split(1:nrow(mlData$x), group)

trainData <- list()
for (cc in 1:10) {
    cat("cc: ", cc, "\n")
    usrIdx <- sort(sample(1:length(usrs), 50000))
    ridx <- unlist(usrs[usrIdx])
    trainData[[cc]] <- list(x = mlData$x[ridx, ],
                            z = mlData$z[ridx, ],
                            y = mlData$y[ridx],
                            users = usrIdx,
                            group = group[ridx]
                            )
}

saveRDS(trainData, "/Shared/ssrivastva/lswasp/ml/data/train_ml.rds")

npart <- 100

set.seed(12345)
for (cc in 1:10) {
    cat("cc: ", cc, "\n")
    repData = trainData[[cc]]
    ngroup <- length(unique(repData$group))
    partData <- vector("list", npart)
    names(partData) <- paste0("k", 1:npart)
    grpSplit <- split(1:nrow(repData$x), repData$group)
    for (ll in 1:npart) {
        grpIdx <- sample(1:length(grpSplit), 2000)
        idx <- unlist(grpSplit[grpIdx])
        partData[[ll]]$nobs <- length(unique(grpIdx))
        partData[[ll]]$x <- repData$x[idx, ]
        partData[[ll]]$y <- repData$y[idx]
        partData[[ll]]$z <- repData$z[idx, ]
        partData[[ll]]$group <- repData$group[idx]
        partData[[ll]]$idx <- idx
        partData[[ll]]$nrep <- ngroup / partData[[ll]]$nobs
    }
    fname <- paste0("/Shared/ssrivastva/lswasp/ml/data/ml_data_rep_", cc, "_k_", npart, ".rds")
    saveRDS(partData, fname)
}


rm(list=ls())
setwd("~/lswasp/ml/code")
trainData <- readRDS("/Shared/ssrivastva/lswasp/ml/data/train_ml.rds")

npart <- 100

set.seed(12345)
for (cc in 1:10) {
    cat("cc: ", cc, "\n")
    repData = trainData[[cc]]
    ngroup <- length(unique(repData$group))
    partData <- vector("list", npart)
    names(partData) <- paste0("k", 1:npart)
    grpSplit <- split(1:nrow(repData$x), repData$group)
    grpPart <- sample(1:npart, length(grpSplit), replace = TRUE)
    for (ll in 1:npart) {
        grpIdx <- which(grpPart == ll)
        idx <- unlist(grpSplit[grpIdx])
        partData[[ll]]$nobs <- length(unique(grpIdx))
        partData[[ll]]$x <- repData$x[idx, ]
        partData[[ll]]$y <- repData$y[idx]
        partData[[ll]]$z <- repData$z[idx, ]
        partData[[ll]]$group <- repData$group[idx]
        partData[[ll]]$idx <- idx
        partData[[ll]]$nrep <- ngroup / partData[[ll]]$nobs
    }
    fname <- paste0("/Shared/ssrivastva/lswasp/ml/data/disjoint_ml_data_rep_", cc, "_k_", npart, ".rds")
    saveRDS(partData, fname)
}

npart <- 200

set.seed(12345)
for (cc in 1:10) {
    cat("cc: ", cc, "\n")
    repData = trainData[[cc]]
    ngroup <- length(unique(repData$group))
    partData <- vector("list", npart)
    names(partData) <- paste0("k", 1:npart)
    grpSplit <- split(1:nrow(repData$x), repData$group)
    grpPart <- sample(1:npart, length(grpSplit), replace = TRUE)
    for (ll in 1:npart) {
        grpIdx <- which(grpPart == ll)
        idx <- unlist(grpSplit[grpIdx])
        partData[[ll]]$nobs <- length(unique(grpIdx))
        partData[[ll]]$x <- repData$x[idx, ]
        partData[[ll]]$y <- repData$y[idx]
        partData[[ll]]$z <- repData$z[idx, ]
        partData[[ll]]$group <- repData$group[idx]
        partData[[ll]]$idx <- idx
        partData[[ll]]$nrep <- ngroup / partData[[ll]]$nobs
    }
    fname <- paste0("/Shared/ssrivastva/lswasp/ml/data/disjoint_ml_data_rep_", cc, "_k_", npart, ".rds")
    saveRDS(partData, fname)
}
