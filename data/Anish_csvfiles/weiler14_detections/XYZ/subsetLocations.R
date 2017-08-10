#!/usr/bin/env Rscript

require(rhdf5)
require(data.table)

f <- 'AnishLocationsFull.h5'

h5ls(f)
dat <- data.table(h5read(f, name = "Locations/locs"))

Exs <- unique(dat$Ex)

ext <- list()
## Gotten from api.boss.neurodata.io
ext$Ex10R55 <- list(xstop = 3409, ystop = 3337, zstop = 70)
ext$Ex12R75 <- list(xstop = 5491, ystop = 4749, zstop = 35)
ext$Ex12R76 <- list(xstop = 5979, ystop = 4872, zstop = 37)
ext$Ex13R51 <- list(xstop = 5184, ystop = 5840, zstop = 30)
ext$Ex2R18C1 <- list(xstop = 2106, ystop = 3236, zstop = 42)
ext$Ex2R18C2 <- list(xstop = 1970, ystop = 3175, zstop = 42)
ext$Ex3R43C1 <- list(xstop = 2101, ystop = 3223, zstop = 69)
ext$Ex3R43C2 <- list(xstop = 1971, ystop = 3164, zstop = 69)
ext$Ex3R43C3 <- list(xstop = 1989, ystop = 3252, zstop = 69)
ext$Ex6R15C1 <- list(xstop = 3208, ystop = 3581, zstop = 30)
ext$Ex6R15C2 <- list(xstop = 3233, ystop = 3636, zstop = 30)


bf <- 5
buffered <- list()
for(i in Exs){
  buffered[[i]] <- dat[Ex == i & 
                  col > bf & col < (ext[[i]]$xstop - (bf + 1)) &
                  row > bf & row < (ext[[i]]$ystop - (bf + 1)) &
                  z > bf & z < (ext[[i]]$zstop - (bf + 1)), ]

  buffered[[i]][, `:=` (Ex = NULL, 
                 Q = NULL)]

  buffered[[i]][, `:=`(col = round(col), # rounding to nearest integer
               row = round(row),
               z   = round(z))]
}


for( j in 1:length(buffered)){
  write.csv(buffered[[j]], file = paste0(names(buffered)[j],"_buff", bf, ".csv"),, row.names = FALSE)
  set.seed(317)
  s <- sample(nrow(buffered[[j]]), 1e4)
  write.csv(buffered[[j]][s,], file = paste0(names(buffered)[j],"_buff", bf, "seed317_samp1e4.csv"),, row.names = FALSE)
}


