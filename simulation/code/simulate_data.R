rm(list=ls())

set.seed(12345)

setwd("~/simulation/code")
library(matrixStats)
library(Matrix)

## See mbest package at http://ptrckprry.com/code/ Perry (2017) in JRSS-B.
genData <- function (ngroup, nobs, nfixef, nranef) {
    ## fixed effects coefficients
    fixef <- rep(c(-2, 2), length = nfixef)
    if (nranef == 3) {
        ranefCorr <- matrix(c(1, -0.4, 0.3,
                              -0.4, 1, 0.001,
                              0.3, 0.001, 1),
                            nranef, nranef)
    } else {
        ranefCorr <- as.matrix(bdiag(rep(list(matrix(c(1, -0.4, 0.3,
                                                       -0.4, 1, 0.001,
                                                       0.3, 0.001, 1),
                                                     3, 3)), 2)))
    }
    ranefCov <- outer(sqrt(1:nranef), sqrt(1:nranef)) * ranefCorr
    ranefCovSqrt <- chol(ranefCov)

    # generate coefficients
    u <- matrix(rnorm(ngroup * nranef), ngroup, nranef)
    ranef <- u %*% ranefCovSqrt

    ## generate group
    suppressWarnings({ # ignore warning about using Walker's alias method
        group <- sample.int(ngroup, nobs, replace=TRUE)
    })

    ## generate feature  matrices with Pr(x[i,j] = +1) = P(x[i,j] = -1) = 1/2,
    x <- matrix(sample(c(-1, +1), nobs * nfixef, replace=TRUE), nobs, nfixef)
    z <- matrix(sample(c(-1, +1), nobs * nranef, replace=TRUE), nobs, nranef)

    ## compute linear predictors and generate observations
    mu <- drop(x %*% fixef) + rowSums(z * ranef[group,])
    y <- rnorm(nobs, mean=mu, sd=1)

    list(ngroup = ngroup, nobs = nobs,
         #fixef = fixef, #ranef = ranef,
         ranefCov = ranefCov,
         ranefCovSqrt = ranefCovSqrt,
         group = group, x = x, z = z, y.mean = mu, y = y)
}

ngroup <- 5e4
nobs <- 5e5
nfixef <- 4
nranef <- 3

train <- list()
for (cc in 1:10) {
    train[[cc]] <- genData(ngroup, nobs, nfixef, nranef)
}

saveRDS(train, "~/simulation/data/full_data.rds")

################################ Disjoint Scheme ##################################
## npart = 20

rm(list=ls())

repData <- readRDS("~/simulation/data/full_data.rds")

set.seed(12345)

ngroup <- 5e4
nobs <- 1e6
nfixef <- 4
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
    fname <- paste0("part_data_rep_", cc, "_k_", npart, ".rds")
    saveRDS(partData, paste0("~/simulation/data/", fname))
}

## npart = 50

rm(list=ls())

repData <- readRDS("~/simulation/data/full_data.rds")

set.seed(12345)

ngroup <- 5e4
nobs <- 1e6
nfixef <- 4
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
    fname <- paste0("part_data_rep_", cc, "_k_", npart, ".rds")
    saveRDS(partData, paste0("~/simulation/data/", fname))
}

## npart = 100

rm(list=ls())

repData <- readRDS("~/simulation/data/full_data.rds")

set.seed(12345)

ngroup <- 5e4
nobs <- 1e6
nfixef <- 4
nranef <- 3

npart <- 100
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
    fname <- paste0("part_data_rep_", cc, "_k_", npart, ".rds")
    saveRDS(partData, paste0("~/simulation/data/", fname))
}

## npart = 150

rm(list=ls())

repData <- readRDS("~/simulation/data/full_data.rds")

set.seed(12345)

ngroup <- 5e4
nobs <- 1e6
nfixef <- 4
nranef <- 3

npart <- 150
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
    fname <- paste0("part_data_rep_", cc, "_k_", npart, ".rds")
    saveRDS(partData, paste0("~/simulation/data/", fname))
}

## npart = 200

rm(list=ls())

repData <- readRDS("~/simulation/data/full_data.rds")

set.seed(12345)

ngroup <- 5e4
nobs <- 1e6
nfixef <- 4
nranef <- 3

npart <- 200
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
    fname <- paste0("part_data_rep_", cc, "_k_", npart, ".rds")
    saveRDS(partData, paste0("~/simulation/data/", fname))
}

################################ Overlapping Scheme ##################################
     
                ############# when sample = 250 ############

## npart = 20
rm(list=ls())

#repData <- readRDS("~/dls/data/full_data.rds")
repData <- readRDS("~/simulation/data/full_data.rds")

set.seed(12345)

ngroup <- 5e4
nobs <- 1e6
nfixef <- 4
nranef <- 3
nresamp <- 250

npart <- 20
partData <- list()

for (cc in 1:10) {
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
  fname <- paste0("part_data_rep_", cc, "_k_", npart, "_dep_samp250.rds")
  #saveRDS(partData, paste0("~/dls/data/", fname))
  saveRDS(partData, paste0("~/simulation/data/",fname))
}

## npart = 50

rm(list=ls())

#repData <- readRDS("~/dls/data/full_data.rds")
repData <- readRDS("~/simulation/data/full_data.rds")

set.seed(12345)

ngroup <- 5e4
nobs <- 1e6
nfixef <- 4
nranef <- 3
nresamp <- 250

npart <- 50
partData <- list()

for (cc in 1:10) {
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
  fname <- paste0("part_data_rep_", cc, "_k_", npart, "_dep_samp250.rds")
  #saveRDS(partData, paste0("~/dls/data/", fname))
  saveRDS(partData, paste0("~/simulation/data/",fname))
}

## npart = 100

rm(list=ls())

#repData <- readRDS("~/dls/data/full_data.rds")
repData <- readRDS("~/simulation/data/full_data.rds")

set.seed(12345)

ngroup <- 5e4
nobs <- 1e6
nfixef <- 4
nranef <- 3
nresamp <- 250

npart <- 100
partData <- list()

for (cc in 1:10) {
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
  fname <- paste0("part_data_rep_", cc, "_k_", npart, "_dep_samp250.rds")
  #saveRDS(partData, paste0("~/dls/data/", fname))
  saveRDS(partData, paste0("~/simulation/data/",fname))
}

## npart = 150

rm(list=ls())

#repData <- readRDS("~/dls/data/full_data.rds")
repData <- readRDS("~/simulation/data/full_data.rds")

set.seed(12345)

ngroup <- 5e4
nobs <- 1e6
nfixef <- 4
nranef <- 3
nresamp <- 250

npart <- 150
partData <- list()

for (cc in 1:10) {
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
  fname <- paste0("part_data_rep_", cc, "_k_", npart, "_dep_samp250.rds")
  #saveRDS(partData, paste0("~/dls/data/", fname))
  saveRDS(partData, paste0("~/simulation/data/",fname))
}

## npart = 200

rm(list=ls())

#repData <- readRDS("~/dls/data/full_data.rds")
repData <- readRDS("~/simulation/data/full_data.rds")

set.seed(12345)

ngroup <- 5e4
nobs <- 1e6
nfixef <- 4
nranef <- 3
nresamp <- 250

npart <- 200
partData <- list()

for (cc in 1:10) {
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
  fname <- paste0("part_data_rep_", cc, "_k_", npart, "_dep_samp250.rds")
  #saveRDS(partData, paste0("~/dls/data/", fname))
  saveRDS(partData, paste0("~/simulation/data/",fname))
}




                ############# when sample = 500 ############

## npart = 20
setwd("/Users/AaronXu/Research_Project/Mixed effect Model/Meeting 10 - real data /simulation/data")
rm(list=ls())

#repData <- readRDS("~/dls/data/full_data.rds")
repData <- readRDS("~/simulation/data/full_data.rds")

set.seed(12345)

ngroup <- 5e4
nobs <- 1e6
nfixef <- 4
nranef <- 3
nresamp <- 500

npart <- 20
partData <- list()

for (cc in 1:10) {
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
  fname <- paste0("part_data_rep_", cc, "_k_", npart, "_dep_samp500.rds")
  #saveRDS(partData, paste0("~/dls/data/", fname))
  saveRDS(partData, paste0("~/simulation/data/",fname))
}

## npart = 50

rm(list=ls())

#repData <- readRDS("~/dls/data/full_data.rds")
repData <- readRDS("~/simulation/data/full_data.rds")

set.seed(12345)

ngroup <- 5e4
nobs <- 1e6
nfixef <- 4
nranef <- 3
nresamp <- 500

npart <- 50
partData <- list()

for (cc in 1:10) {
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
  fname <- paste0("part_data_rep_", cc, "_k_", npart, "_dep_samp500.rds")
  #saveRDS(partData, paste0("~/dls/data/", fname))
  saveRDS(partData, paste0("~/simulation/data/",fname))
}

## npart = 100

rm(list=ls())

#repData <- readRDS("~/dls/data/full_data.rds")
repData <- readRDS("~/simulation/data/full_data.rds")

set.seed(12345)

ngroup <- 5e4
nobs <- 1e6
nfixef <- 4
nranef <- 3
nresamp <- 500

npart <- 100
partData <- list()

for (cc in 1:10) {
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
  fname <- paste0("part_data_rep_", cc, "_k_", npart, "_dep_samp500.rds")
  #saveRDS(partData, paste0("~/dls/data/", fname))
  saveRDS(partData, paste0("~/simulation/data/",fname))
}

## npart = 150

rm(list=ls())

#repData <- readRDS("~/dls/data/full_data.rds")
repData <- readRDS("~/simulation/data/full_data.rds")

set.seed(12345)

ngroup <- 5e4
nobs <- 1e6
nfixef <- 4
nranef <- 3
nresamp <- 500

npart <- 150
partData <- list()

for (cc in 1:10) {
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
  fname <- paste0("part_data_rep_", cc, "_k_", npart, "_dep_samp500.rds")
  #saveRDS(partData, paste0("~/dls/data/", fname))
  saveRDS(partData, paste0("~/simulation/data/",fname))
}

## npart = 200

rm(list=ls())

#repData <- readRDS("~/dls/data/full_data.rds")
repData <- readRDS("~/simulation/data/full_data.rds")

set.seed(12345)

ngroup <- 5e4
nobs <- 1e6
nfixef <- 4
nranef <- 3
nresamp <- 500

npart <- 200
partData <- list()

for (cc in 1:10) {
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
  fname <- paste0("part_data_rep_", cc, "_k_", npart, "_dep_samp500.rds")
  #saveRDS(partData, paste0("~/dls/data/", fname))
  saveRDS(partData, paste0("~/simulation/data/",fname))
}
