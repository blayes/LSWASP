setwd("~/us_natality/data")


################################ Disjoint Scheme ##################################
## npart = 20

rm(list=ls())

repData <- readRDS("train_bwt.rds")

set.seed(12345)

nfixef <- 14
nranef <- 3

npart <- 20
partData <- list()

for (cc in 1:10) {
    partData <- vector("list", npart)
    names(partData) <- paste0("k", 1:npart)
    lst <- repData[[cc]]
    grpSplit <- split(1:nrow(lst$x), lst$group)
    partsIdx <- sample(1:npart, length(grpSplit), replace = TRUE)
    for (ll in 1:npart) {
        grpIdx <- which(partsIdx == ll)
        idx <- unlist(grpSplit[grpIdx])
        partData[[ll]]$nobs <- length(unique(grpIdx))
        partData[[ll]]$x <- lst$x[idx, ]
        partData[[ll]]$y <- lst$y[idx]
        partData[[ll]]$z <- lst$z[idx, ]
        partData[[ll]]$group <- lst$group[idx]
        partData[[ll]]$idx <- idx
        partData[[ll]]$nrep <- ngroup / partData[[ll]]$nobs
    }
    fname <- paste0("part_real_data_rep_", cc, "_k_", npart, ".rds")
    saveRDS(partData, paste0("~/us_natality/data/", fname))
}


## npart = 30

rm(list=ls())

repData <- readRDS("train_bwt.rds")

set.seed(12345)

nfixef <- 14
nranef <- 3

npart <- 50
partData <- list()

for (cc in 1:10) {
    partData <- vector("list", npart)
    names(partData) <- paste0("k", 1:npart)
    lst <- repData[[cc]]
    grpSplit <- split(1:nrow(lst$x), lst$group)
    partsIdx <- sample(1:npart, length(grpSplit), replace = TRUE)
    for (ll in 1:npart) {
        grpIdx <- which(partsIdx == ll)
        idx <- unlist(grpSplit[grpIdx])
        partData[[ll]]$nobs <- length(unique(grpIdx))
        partData[[ll]]$x <- lst$x[idx, ]
        partData[[ll]]$y <- lst$y[idx]
        partData[[ll]]$z <- lst$z[idx, ]
        partData[[ll]]$group <- lst$group[idx]
        partData[[ll]]$idx <- idx
        partData[[ll]]$nrep <- ngroup / partData[[ll]]$nobs
    }
    fname <- paste0("part_real_data_rep_", cc, "_k_", npart, ".rds")
    saveRDS(partData, paste0("~/us_natality/data/", fname))
}




################################ Overlap Scheme ##################################

             ########## sample = 250 ##########

## npart = 20
rm(list=ls())

repData <- readRDS("train_bwt.rds")

set.seed(12345)

nfixef <- 14
nranef <- 3
nresamp <- 250

npart <- 20
partData <- list()

for (cc in 1:10) {
  nobs <- nrow(repData[[cc]]$x)
  ngroup <- length(unique(repData[[cc]]$group))  
  partData <- vector("list", npart)
  names(partData) <- paste0("k", 1:npart)
  lst <- repData[[cc]]
  grpSplit <- split(1:nrow(lst$x), lst$group)
  for (ll in 1:npart) {
    grpIdx <- sample(1:length(grpSplit), nresamp, replace = TRUE)
    idx <- unlist(grpSplit[grpIdx])
    partData[[ll]]$nobs <- length(unique(grpIdx))
    partData[[ll]]$x <- lst$x[idx, ]
    partData[[ll]]$y <- lst$y[idx]
    partData[[ll]]$z <- lst$z[idx, ]
    partData[[ll]]$group <- lst$group[idx]
    partData[[ll]]$idx <- idx
    partData[[ll]]$nrep <- ngroup / partData[[ll]]$nobs
  }
  fname <- paste0("part_real_data_rep_", cc, "_k_", npart, "_dep250.rds")
  saveRDS(partData, paste0("~/us_natality/data/",fname))
}

## npart = 30
rm(list=ls())

repData <- readRDS("train_bwt.rds")

set.seed(12345)

nfixef <- 14
nranef <- 3
nresamp <- 250

npart <- 30
partData <- list()

for (cc in 1:10) {
  nobs <- nrow(repData[[cc]]$x)
  ngroup <- length(unique(repData[[cc]]$group))  
  partData <- vector("list", npart)
  names(partData) <- paste0("k", 1:npart)
  lst <- repData[[cc]]
  grpSplit <- split(1:nrow(lst$x), lst$group)
  for (ll in 1:npart) {
    grpIdx <- sample(1:length(grpSplit), nresamp, replace = TRUE)
    idx <- unlist(grpSplit[grpIdx])
    partData[[ll]]$nobs <- length(unique(grpIdx))
    partData[[ll]]$x <- lst$x[idx, ]
    partData[[ll]]$y <- lst$y[idx]
    partData[[ll]]$z <- lst$z[idx, ]
    partData[[ll]]$group <- lst$group[idx]
    partData[[ll]]$idx <- idx
    partData[[ll]]$nrep <- ngroup / partData[[ll]]$nobs
  }
  fname <- paste0("part_real_data_rep_", cc, "_k_", npart, "_dep250.rds")
  saveRDS(partData, paste0("~/us_natality/data/",fname))
}


## npart = 50
rm(list=ls())

repData <- readRDS("train_bwt.rds")

set.seed(12345)

nfixef <- 14
nranef <- 3
nresamp <- 250

npart <- 50
partData <- list()

for (cc in 1:10) {
  nobs <- nrow(repData[[cc]]$x)
  ngroup <- length(unique(repData[[cc]]$group))  
  partData <- vector("list", npart)
  names(partData) <- paste0("k", 1:npart)
  lst <- repData[[cc]]
  grpSplit <- split(1:nrow(lst$x), lst$group)
  for (ll in 1:npart) {
    grpIdx <- sample(1:length(grpSplit), nresamp, replace = TRUE)
    idx <- unlist(grpSplit[grpIdx])
    partData[[ll]]$nobs <- length(unique(grpIdx))
    partData[[ll]]$x <- lst$x[idx, ]
    partData[[ll]]$y <- lst$y[idx]
    partData[[ll]]$z <- lst$z[idx, ]
    partData[[ll]]$group <- lst$group[idx]
    partData[[ll]]$idx <- idx
    partData[[ll]]$nrep <- ngroup / partData[[ll]]$nobs
  }
  fname <- paste0("part_real_data_rep_", cc, "_k_", npart, "_dep250.rds")
  saveRDS(partData, paste0("~/us_natality/data/",fname))
}



## npart = 100
rm(list=ls())

repData <- readRDS("train_bwt.rds")

set.seed(12345)

nfixef <- 14
nranef <- 3
nresamp <- 250

npart <- 100
partData <- list()

for (cc in 1:10) {
  nobs <- nrow(repData[[cc]]$x)
  ngroup <- length(unique(repData[[cc]]$group))  
  partData <- vector("list", npart)
  names(partData) <- paste0("k", 1:npart)
  lst <- repData[[cc]]
  grpSplit <- split(1:nrow(lst$x), lst$group)
  for (ll in 1:npart) {
    grpIdx <- sample(1:length(grpSplit), nresamp, replace = TRUE)
    idx <- unlist(grpSplit[grpIdx])
    partData[[ll]]$nobs <- length(unique(grpIdx))
    partData[[ll]]$x <- lst$x[idx, ]
    partData[[ll]]$y <- lst$y[idx]
    partData[[ll]]$z <- lst$z[idx, ]
    partData[[ll]]$group <- lst$group[idx]
    partData[[ll]]$idx <- idx
    partData[[ll]]$nrep <- ngroup / partData[[ll]]$nobs
  }
  fname <- paste0("part_real_data_rep_", cc, "_k_", npart, "_dep250.rds")
  saveRDS(partData, paste0("~/us_natality/data/",fname))
}




        ########## sample = 500 ##########

## npart = 20
rm(list=ls())

repData <- readRDS("train_bwt.rds")

set.seed(12345)

nfixef <- 14
nranef <- 3
nresamp <- 500

npart <- 20
partData <- list()

for (cc in 1:10) {
  nobs <- nrow(repData[[cc]]$x)
  ngroup <- length(unique(repData[[cc]]$group))  
  partData <- vector("list", npart)
  names(partData) <- paste0("k", 1:npart)
  lst <- repData[[cc]]
  grpSplit <- split(1:nrow(lst$x), lst$group)
  for (ll in 1:npart) {
    grpIdx <- sample(1:length(grpSplit), nresamp, replace = TRUE)
    idx <- unlist(grpSplit[grpIdx])
    partData[[ll]]$nobs <- length(unique(grpIdx))
    partData[[ll]]$x <- lst$x[idx, ]
    partData[[ll]]$y <- lst$y[idx]
    partData[[ll]]$z <- lst$z[idx, ]
    partData[[ll]]$group <- lst$group[idx]
    partData[[ll]]$idx <- idx
    partData[[ll]]$nrep <- ngroup / partData[[ll]]$nobs
  }
  fname <- paste0("part_real_data_rep_", cc, "_k_", npart, "_dep500.rds")
  saveRDS(partData, paste0("~/us_natality/data/",fname))
}


## npart = 30
rm(list=ls())

repData <- readRDS("train_bwt.rds")

set.seed(12345)

nfixef <- 14
nranef <- 3
nresamp <- 500

npart <- 30
partData <- list()

for (cc in 1:10) {
  nobs <- nrow(repData[[cc]]$x)
  ngroup <- length(unique(repData[[cc]]$group))  
  partData <- vector("list", npart)
  names(partData) <- paste0("k", 1:npart)
  lst <- repData[[cc]]
  grpSplit <- split(1:nrow(lst$x), lst$group)
  for (ll in 1:npart) {
    grpIdx <- sample(1:length(grpSplit), nresamp, replace = TRUE)
    idx <- unlist(grpSplit[grpIdx])
    partData[[ll]]$nobs <- length(unique(grpIdx))
    partData[[ll]]$x <- lst$x[idx, ]
    partData[[ll]]$y <- lst$y[idx]
    partData[[ll]]$z <- lst$z[idx, ]
    partData[[ll]]$group <- lst$group[idx]
    partData[[ll]]$idx <- idx
    partData[[ll]]$nrep <- ngroup / partData[[ll]]$nobs
  }
  fname <- paste0("part_real_data_rep_", cc, "_k_", npart, "_dep500.rds")
  saveRDS(partData, paste0("~/us_natality/data/",fname))
}


## npart = 50
rm(list=ls())

repData <- readRDS("train_bwt.rds")

set.seed(12345)

nfixef <- 14
nranef <- 3
nresamp <- 500

npart <- 50
partData <- list()

for (cc in 1:10) {
  nobs <- nrow(repData[[cc]]$x)
  ngroup <- length(unique(repData[[cc]]$group))  
  partData <- vector("list", npart)
  names(partData) <- paste0("k", 1:npart)
  lst <- repData[[cc]]
  grpSplit <- split(1:nrow(lst$x), lst$group)
  for (ll in 1:npart) {
    grpIdx <- sample(1:length(grpSplit), nresamp, replace = TRUE)
    idx <- unlist(grpSplit[grpIdx])
    partData[[ll]]$nobs <- length(unique(grpIdx))
    partData[[ll]]$x <- lst$x[idx, ]
    partData[[ll]]$y <- lst$y[idx]
    partData[[ll]]$z <- lst$z[idx, ]
    partData[[ll]]$group <- lst$group[idx]
    partData[[ll]]$idx <- idx
    partData[[ll]]$nrep <- ngroup / partData[[ll]]$nobs
  }
  fname <- paste0("part_real_data_rep_", cc, "_k_", npart, "_dep500.rds")
  saveRDS(partData, paste0("~/us_natality/data/",fname))
}



## npart = 100
rm(list=ls())

repData <- readRDS("train_bwt.rds")

set.seed(12345)

nfixef <- 14
nranef <- 3
nresamp <- 500

npart <- 100
partData <- list()

for (cc in 1:10) {
  nobs <- nrow(repData[[cc]]$x)
  ngroup <- length(unique(repData[[cc]]$group))  
  partData <- vector("list", npart)
  names(partData) <- paste0("k", 1:npart)
  lst <- repData[[cc]]
  grpSplit <- split(1:nrow(lst$x), lst$group)
  for (ll in 1:npart) {
    grpIdx <- sample(1:length(grpSplit), nresamp, replace = TRUE)
    idx <- unlist(grpSplit[grpIdx])
    partData[[ll]]$nobs <- length(unique(grpIdx))
    partData[[ll]]$x <- lst$x[idx, ]
    partData[[ll]]$y <- lst$y[idx]
    partData[[ll]]$z <- lst$z[idx, ]
    partData[[ll]]$group <- lst$group[idx]
    partData[[ll]]$idx <- idx
    partData[[ll]]$nrep <- ngroup / partData[[ll]]$nobs
  }
  fname <- paste0("part_real_data_rep_", cc, "_k_", npart, "_dep500.rds")
  saveRDS(partData, fname)
}

