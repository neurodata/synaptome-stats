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
args  <- commandArgs(trailingOnly = TRUE)
input <- args[1]
token <- args[2]
output<- as.character(args[3])
bf    <- args[4]
scaleMeans <- read.csv(args[5], header = TRUE)
scaleSDs <- read.csv(args[6], header = TRUE)

if(FALSE){
  input <- "~/tmps/k15Run20170518/data/outputs/zscored/outputs/cubesC2111.h5"
  token <- "kristina15"
  output<- "~/tmps/k15Run20170518/data/outpus/zscored/outputs/rsynC2111_"
  bf    <- 5
  scaleMeans <- h5read("~/tmps/k15Run20170518/data/inputs/K15F0zscored_seed1234_samp1e4.h5", 
                       name= "data/ColMeans")
  scaleSDs <- h5read("~/tmps/k15Run20170518/data/inputs/K15F0zscored_seed1234_samp1e4.h5", 
                       name= "data/ColSDs")
  H5close()

  hmcO <- readRDS("~/tmps/k15Run20170518/data/outputs/zscored/hmc.rds")
}

#TOKEN1= (sys.argv[1]) ## e.g. "Ex10R55"  "kristina15"
#INPUT = (sys.argv[2]) ## a csv file of locations in (x, y, z) order.
#NAME  = (sys.argv[3]) ## A base name for output
#bf = (sys.argv[4]) ## Buffer size around each center point

suppressPackageStartupMessages(require(rhdf5))
suppressWarnings(suppressPackageStartupMessages(require(meda)))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(gridExtra))
suppressPackageStartupMessages(require(foreach))

cnames <- h5read(input, name = "Channels")
abc <- abbreviate(gsub("_", "", cnames), minlength = 3)
abcnames <- factor(abc, levels = abc)
dat <- h5read(input, name = paste0(token, "_cubes"))

dimnames(dat) <- list(1:dim(dat)[1], cnames, NULL, NULL, NULL)
datO <- dat
loc <- h5read(input, name = "Locations")
H5close()

## Type for kristina data full only!!!
#type <- c("ot", "ex", "in", "ot", 
#          "in", "ex", "ot", "ot", 
#          "ot", "ex", "ex", "ex", 
#          "ex", "in", "in", "ex", 
#          "ex", "ex", "ot", "in", 
#          "ex", "in", "in")
type <- c(rep("ex", 11), rep("in", 8), rep("ot", 4))

tab <- table(type)
ctype <- data.frame(cnames, type)
ord <- order(type)

inpal <- colorpanel(255,"black", "red")
expal <- colorpanel(255, "black", "green")
otpal <- colorpanel(255, "black", "blue")


## Calculate feature F0
F0o <- 
  foreach(j = 1:dim(dat)[1], .combine = rbind) %:%
  foreach(i = 1:dim(dat)[2], .combine=cbind) %do% {
      f0 <- as.numeric(sum(dat[j,i,,,]))
}


colnames(F0o) <- cnames

if(!all(colnames(F0o) == colnames(scaleMeans)) |
  !all(colnames(F0o) == colnames(scaleSDs))){
  stop("Column names do not match!") 
}


F0 <- scale(F0o, center = as.numeric(scaleMeans), scale = as.numeric(scaleSDs))

F0 <- lapply(1:nrow(F0), function(r) F0[r,])
### Scaling for color bar according to z-scores
F0min <-min(Reduce(c, Reduce(c, F0)))
F0max <-max(Reduce(c, Reduce(c, F0)))

## Normalize across channels for visualization purposes
#for(j in 1:dim(dat)[2]){
#    a1 <- (dat[,j,,,])
#    if(!all(a1 == 0)){
#      #a2 <- scale(a1, center = TRUE, scale = TRUE)
#      #a2 <- (a1 - min(a1))/(max(a1) - min(a1))
#      #a3 <- log10(a1 + 1)
#      sm <- as.numeric(scaleMeans[,j])
#      ss <- as.numeric(scaleSDs[,j])
#      a2 <- (a1 - sm)/ss
#      dat[,j,,,] <- a2
#    }
#  }

rr <- foreach(l = 1:dim(dat)[1]) %do% {
  r <- dat[l,,,,]
  mr <- melt(r)
  names(mr) <- c("ch", "x", "y", "z", "value")
  mr$z <- mr$z - (max(mr$z) +1)/2

  F0z <- Reduce('rbind', 
           lapply(names(F0[[l]]), 
             function(n){
             Fz <- mr[mr$ch == n,]
             #Fz$z <- as.numeric(max(mr$z)+1)
             Fz$F0 <- round(as.numeric(F0[[l]][n]),4)
             Fz
             }))

  #mr <- rbind(mr, F0z)
  #ch <- mr$ch
  #mr$type <- ctype$type[ch]
  #mr$ch <- abcnames[ch]
  #mr
  ch <- F0z$ch
  F0z$type <- ctype$type[ch]
  F0z$ch <- abcnames[ch]
  F0z
}


names(rr) <- 1:length(rr)


th <- theme(axis.text = element_blank(), axis.ticks = element_blank(),
            axis.title.y = element_blank(), axis.title.x = element_blank(),
            legend.position="bottom", legend.key.size = unit(2,'lines'),
            panel.spacing = unit(0, "lines"))

mmex <- min(sapply(rr, function(am){ min(am[am$type == "ex",]$value) }))
MMex <- max(sapply(rr, function(am){ max(am[am$type == "ex",]$value) }))
mmin <- min(sapply(rr, function(am){ min(am[am$type == "in",]$value) }))
MMin <- max(sapply(rr, function(am){ max(am[am$type == "in",]$value) }))
mmot <- min(sapply(rr, function(am){ min(am[am$type == "ot",]$value) }))
MMot <- max(sapply(rr, function(am){ max(am[am$type == "ot",]$value) }))

ggpsep <- list()
for(k in 1:length(rr)){
  mr <- rr[[k]]
  pex <- 
  ggplot(mr[mr$type == "ex",], 
    aes(x,y, group = type, fill = value)) +
    geom_raster() + 
    scale_y_reverse() + 
    #facet_grid(ch + F0 ~ z, labeller = label_both) +
    facet_grid(ch ~ z, labeller = label_both) +
    scale_fill_gradient(low = "black", high = "green", limits = c(mmex,MMex)) + th

  pexF0 <- 
    ggplot(mr[mr$type == "ex" ,], 
      aes(x,y, group = type, fill = F0)) +
      geom_raster() + 
      scale_y_reverse() + 
      facet_grid(ch + F0 ~ type, labeller = label_both) +
      scale_fill_gradient2(low = "darkorchid4", 
                              mid = "gray99", 
                              high = "darkorange3",
                              midpoint = 0, limits=c(F0min, F0max)) + th
      #scale_fill_gradient(low = "#", high = "orange") + th


  pin <- 
  ggplot(mr[mr$type == "in",], 
    aes(x,y, group = type, fill = value)) +
    geom_raster() + 
    scale_y_reverse() + 
    #facet_grid(ch + F0 ~ z, labeller = label_both) +
    facet_grid(ch ~ z, labeller = label_both) +
    scale_fill_gradient(low = "black", high = "red", limits = c(mmin,MMin)) + th

  pinF0 <- 
    ggplot(mr[mr$type == "in" ,], 
      aes(x,y, group = type, fill = F0)) +
      geom_raster() + 
      scale_y_reverse() + 
      facet_grid(ch + F0 ~ type, labeller = label_both) +
      scale_fill_gradient2(low = "darkorchid4", 
                              mid = "gray99", 
                              high = "darkorange3",
                              midpoint = 0, limits = c(F0min, F0max)) + th
      #scale_fill_gradient(low = "black", high = "white") + th

  pot <- 
    ggplot(mr[mr$type == "ot",], 
    aes(x,y, group = type, fill = value)) +
    geom_raster() +
    scale_y_reverse() + 
    facet_grid(ch ~ z, labeller = label_both) +
    scale_fill_gradient(low = "black", high = "blue", limits = c(mmot,MMot)) + th

  potF0 <- 
    ggplot(mr[mr$type == "ot" ,], 
      aes(x,y, group = type, fill = F0)) +
      geom_raster() + 
      scale_y_reverse() + 
      facet_grid(ch + F0 ~ type, labeller = label_both) +
      scale_fill_gradient2(low = "darkorchid4", 
                              mid = "gray99", 
                              high = "darkorange3",
                              midpoint = 0, limits = c(F0min, F0max)) + th
      #scale_fill_gradient(low = "black", high = "white") + th

  xyzloc = c(x = loc[k,1], 
          y = loc[k,2], 
          z = loc[k,3])

  lz <- length(range(mr$z)[1]:range(mr$z)[2])

  laysep <- c()
  kj <- 1
  for(i in seq(1,5,2)){
    inner <- c()
    for(j in 1:table(type)[kj]){
      inner <- rbind(inner, c(rep(i, lz),rep(i+1,3)))
    }
  kj <- kj +1
  laysep <- rbind(laysep, inner)
  }

  tmpd <- tempdir()
  png(file.path(tmpd, "tmp.png"))
  psep <- arrangeGrob(pex, pexF0, pin, pinF0, pot, potF0, layout_matrix = laysep)
  dev.off()

  ggpsep[[k]] <- list(xyzloc=xyzloc, p=psep, 
                   pex=pex, pexF0=pexF0,
                   pin=pin, pinF0=pinF0,
                   pot=pot, potF0=potF0)
}

## DEBUG
if(FALSE){
  pdf("~/Desktop/tmp.pdf", 
      height = sum(tab)*1.5, 
      width = (dim(dat)[5]+1)*1.5)
  plot(ggpsep[[2]][[2]])
  dev.off()
}

save(ggpsep, file = paste0(output, "synaptograms.RData"))

invisible(
mapply(function(grob,locs){
  outName <- paste0(output, token,
               paste0(c("x", "y", "z"),locs, collapse="_"),
               ".png")
  png(outName, 
      height = sum(tab)*1.5, width = (dim(dat)[5]+1)*1.5, 
      units = "in", res = 300)
  plot(grob)
  dev.off()},
  lapply(ggpsep, '[[', 2), lapply(ggpsep, '[[', 1))
)
