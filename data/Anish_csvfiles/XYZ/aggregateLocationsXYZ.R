#!/usr/bin/env Rscript
require(data.table)
require(rhdf5)

basedir <- getwd()

dir0 <- grep("^Ex", dir(), value = TRUE)
dir0

d0 <- list()

for(i in 1:length(dir0)){
  subd <- dir(dir0[[i]], full.name = TRUE)

  Ex <- dir0[[i]]

  vol <-  gsub("\\.csv", "", 
            gsub("\\/", "_", dir(dir0[[i]])))

  for(j in 1:length(subd)){
    d0[[length(d0)+1]] <- 
      data.table(cbind(read.csv(subd[[j]]), Ex = Ex, Q = vol[[j]], stringsAsFactors=FALSE), stringsAsFactors=FALSE)
  }
}

tmp <- Reduce("rbind", d0)

DFull <- tmp
ExLev <- DFull$Ex
Qlev <- DFull$Q

fname <- "AnishLocationsFull.h5"
clean <- h5createFile(fname)
clean <- h5createGroup(fname, "Locations")
clean <- h5createGroup(fname, "Levels")
h5write(DFull, fname, "Locations/locs")
h5write(ExLev, fname, "Levels/Ex")
h5write(Qlev, fname, "Levels/Q")
h5ls(fname)
H5close()


ex10r55 <- DFull[Ex == "Ex10R55",]
ex10r55[, Ex := NULL]

ex10r55[, `:=`(col = round(col), # rounding to nearest integer
               row = round(row),
               z   = round(z))]

fname <- "Ex10R55_BufferedPM5.h5"
clean <- h5createFile(fname)
#clean <- h5createGroup(fname, "Locations")
h5write(ex10r55, fname, "Locations")
h5ls(fname)
H5close()

set.seed(317)
s <- sample(nrow(ex10r55))
write.csv(ex10r55[s,1:3], file = "ex10r55_buffered_seed317.csv", row.names = FALSE)


### Area extent of experiment Ex2R18C1 given by 
### http://openconnecto.me/ocp/ca/Ex2R18C1/info/
### is x: [0,2106), y: [0,3236), z: [0,42)
ex2r18c1 <- DFull[Ex == "Ex2R18C1", 
                  x > 5 && x < (2106-(5+1)) &&
                  y > 5 && y < (3236-(5+1)) &&
                  z > 5 && z < (42-(5+1))]

ex2r18c1[, Ex := NULL]

ex2r18c1[, `:=`(col = round(col), # rounding to nearest integer
               row = round(row),
               z   = round(z))]

fname <- "lEx2R18C1_BufferedPM5.h5"
clean <- h5createFile(fname)
#clean <- h5createGroup(fname, "Locations")
h5write(ex2r18c1, fname, "Locations")
h5ls(fname)
H5close()

set.seed(317)
s <- sample(nrow(ex2r18c1))
write.csv(ex2r18c1[s,1:3], file = "ex2r18c1_buffered_seed317.csv", row.names = FALSE)
