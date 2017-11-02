#!/usr/bin/env Rscript
###
### Rscript to plot synaptograms from
### hdf5 file 
###
###
### Jesse Leigh Patsolic 
### 2017 <jpatsolic@jhu.edu>
### S.D.G 
#

require(ggplot2)
require(foreach)
require(doMC)
require(rhdf5)
require(reshape2)
require(gplots)
require(spatialfil)
require(raster)
require(mvtnorm)
require(compiler)

registerDoMC(8)

mu <- c(0,0,0)
sigma <-  diag(3) * c(2000,2000,15)
sigma


xseq <- seq(-108,108, by = 1)
yseq <- seq(-108,108, by = 1)
zseq <- seq(-5,5, by = 1)

grL <- expand.grid(xseq, yseq, zseq)
grU <- grL + 1
head(grL)
head(grU)

M <- array(NA, dim = c(length(xseq),length(yseq), length(zseq)))


pmv <- cmpfun(pmvnorm)

gL <- as.list(as.data.frame(t(grL)))
gU <- as.list(as.data.frame(t(grU)))

itL <- iter(gL)
itU <- iter(gU)

#la <- foreach(i = 1:nrow(grL), .combine = c) %dopar% {
#  pmvnorm(mean = mu, sigma=sigma, lower=as.numeric(grL[i,]), upper=as.numeric(grU[i,]))
#}

system.time({
  la <- foreach(L = itL, .combine = c) %dopar% {
    pmvnorm(mean = mu, sigma=sigma, lower=as.numeric(L), upper=as.numeric(L + 1))
    pmvnorm(mean = mu, sigma=sigma, lower=as.numeric(L), upper=as.numeric(L + 1))
  }
}
)


itL <- iter(gL[1:1000])
foreach(L = itL) %do% {
  L 
}

system.time({
  la <- foreach(L = itL, .combine = c) %dopar% {
    pmv(mean = mu, sigma=sigma, lower=as.numeric(L), upper=as.numeric(L + 1))
  }
}
)



A <- array(la/max(la), dim = c(length(xseq),length(yseq), length(zseq)))
#saveRDS(A, file = "maskGaussian_mu0_sigma2000_2000_15.rds")

th <- theme(axis.text = element_blank(), axis.ticks = element_blank(),
            axis.title.y = element_blank(), axis.title.x = element_blank(),
            legend.position="bottom", legend.key.size = unit(2,'lines'),
            panel.spacing = unit(0, "lines"))

pdat <- melt(A)
ggplot(pdat, 
aes(Var1,Var2, group = Var3, fill = value)) +
geom_raster() + 
scale_y_reverse() + 
facet_grid(. ~ Var3) +
#scale_fill_gradient(low = "black", high = "green") + th
scale_fill_gradient(low = "black", high = "white") + th


f <- 'testCubes.csv.h5'
f <- 'collman14v2_annotationCubes.csv.h5'
h5ls(f)


cubes <- h5read(f, name = "collman14v2_cubes")
cubes <- cubes[,,,1:5,]
tmp <- cubes

for(i in 1:dim(tmp)[4]){
  for(j in 1:dim(tmp)[5]){
    a <- tmp[,,,i,j] 
    b <- a * A
    tmp[,,,i,j] <- b
  }}


plot(raster(a[,,6]))
plot(raster(b[,,6]))

rawrr <- foreach(l = 1:dim(cubes)[4]) %do% {
  r <- cubes[,,,l,]
  mr <- melt(r)
  colnames(mr) <- c("x", "y", "z", "ch", "value")
  mr$z <- mr$z - (max(mr$z) +1)/2
  mr
}

rr <- foreach(l = 1:dim(tmp)[4]) %do% {
  r <- tmp[,,,l,]
  mr <- melt(r)
  colnames(mr) <- c("x", "y", "z", "ch", "value")
  mr$z <- mr$z - (max(mr$z) +1)/2
  mr
}


th <- theme(axis.text = element_blank(), axis.ticks = element_blank(),
            axis.title.y = element_blank(), axis.title.x = element_blank(),
            legend.position="bottom", legend.key.size = unit(2,'lines'),
            panel.spacing = unit(0, "lines"))


pg <- list()
for(i in 1:length(rr)){
  pdat <- rr[[i]]
  pg[[i]] <- 
    ggplot(pdat, 
    aes(x,y, group = ch, fill = value)) +
    geom_raster() + 
    scale_y_reverse() + 
    facet_grid(ch ~ z, labeller = label_both) +
    #scale_fill_gradient(low = "black", high = "green") + th
    scale_fill_gradient(low = "black", high = "white") + th
}

rawpg <- list()
for(i in 1:length(rawrr)){
  asd <- rawrr[[i]]
  rawpg[[i]] <- 
    ggplot(asd, 
    aes(x,y, group = ch, fill = value)) +
    geom_raster() + 
    scale_y_reverse() + 
    facet_grid(ch ~ z, labeller = label_both) +
    #scale_fill_gradient(low = "black", high = "green") + th
    scale_fill_gradient(low = "black", high = "white") + th
} 

show(rawpg[[1]])
dev.new()
show(pg[[1]])

#
