#!/usr/bin/env Rscript
require(data.table)
require(rhdf5)

basedir <- getwd()

dir0 <- grep("^Ex", dir(), value = TRUE)
dir0

d1 <- list(list())
d0 <- list()

for(i in 1:length(dir0)){
  subd <- dir(dir0[[i]], full.name = TRUE)

  Ex <- dir0[[i]]

  vol <-  gsub("\\.csv", "", 
            gsub("\\/", "_", dir(dir0[[i]])))

  for(j in 1:length(subd)){
    d0[[length(d0)+1]] <- 
      as.data.table(cbind(Ex = Ex, Q = vol[[j]], read.csv(subd[[j]])))
  }
}
DFull <- as.data.table(Reduce("rbind", d0))

ExLev <- levels(DFull$Ex)
Qlev <- levels(DFull$Q)

fname <- "AnishLocationsFull.h5"
clean <- h5createFile(fname)
clean <- h5createGroup(fname, "Locations")
clean <- h5createGroup(fname, "Levels")
h5write(DFull, fname, "Locations/locs")
h5write(ExLev, fname, "Levels/Ex")
h5write(Qlev, fname, "Levels/Q")
h5ls(fname)
H5close()

