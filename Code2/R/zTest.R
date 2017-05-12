suppressPackageStartupMessages(require(rhdf5))
suppressPackageStartupMessages(require(meda))
require(gridExtra)
require(foreach)
require(raster)

fname <- "../Python/kristina15rscaled317samp5000_synaptograms11cubes5NN_test7.h5"

h5ls(fname) 
cnames <- h5read(fname, name = "Channels")
abcnames <- abbreviate(gsub("_", "", cnames), minlength = 3)
dat <- h5read(fname, name = "kristina15")
dimnames(dat) <- list(1:dim(dat)[1], cnames, NULL, NULL, NULL)
datO <- dat
loc <- h5read(fname, name = "Locations")
H5close()

type <- c("in", "ex", "in", "ot", 
          "in", "ex", "ot", "in", 
          "ot", "ex", "ex", "ex", 
          "ex", "in", "in", "ex", 
          "ex", "ex", "ot", "in", 
          "ex", "ex", "in")

ctype <- data.frame(cnames, type)
ord <- order(type)

inpal <- colorpanel(255,"black", "red")
expal <- colorpanel(255, "black", "green")
otpal <- colorpanel(255, "black", "blue")


## Normalize across channels
for(j in 1:dim(dat)[2]){
    a1 <- (dat[,j,,,])
    if(!all(a1 == 0)){
      a2 <- scale(a1, center = TRUE, scale = TRUE)
      #a2 <- (a1 - min(a1))/(max(a1) - min(a1))
      #a3 <- log10(a1 + 1)
      dat[,j,,,] <- a2
    }
  }

F0 <- 
  foreach(j = 1:dim(dat)[1]) %:%
  foreach(i = 1:dim(dat)[2]) %do% {
      f0 <- sum(dat[j,i,,,])
}

names(F0) <- 1:length(F0)
for(i in 1:length(F0)){
  names(F0[[i]]) <- cnames
}

j <- dat

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

#registerDoMC(6)
#foreach(k = 1:length(rr))%dopar%{
mmex <- min(sapply(rr, function(am){ min(am[am$type == "ex",]$value) }))
MMex <- max(sapply(rr, function(am){ max(am[am$type == "ex",]$value) }))
mmin <- min(sapply(rr, function(am){ min(am[am$type == "in",]$value) }))
MMin <- max(sapply(rr, function(am){ max(am[am$type == "in",]$value) }))
mmot <- min(sapply(rr, function(am){ min(am[am$type == "ot",]$value) }))
MMot <- max(sapply(rr, function(am){ max(am[am$type == "ot",]$value) }))

#for(k in 1:length(rr)){
for(k in 1:5){
  mr <- rr[[k]]
#mm <- 0
#MM <- 1
pex <- 
ggplot(mr[mr$type == "ex",], 
  aes(x,y, group = factor(type), fill = value)) +
  geom_raster() + 
  scale_y_reverse() + 
  facet_grid(ch+F0 ~ z, labeller = label_both) +
  #scale_fill_gradient(low = "black", high = "green", limits = c(mmex,MMex)) + th
  scale_fill_gradient(low = "black", high = "green") + th

pin <- 
ggplot(mr[mr$type == "in",], 
  aes(x,y, group = factor(type), fill = value)) +
  geom_raster() + 
  scale_y_reverse() + 
  facet_grid(ch+F0 ~ z, labeller = label_both) +
  #scale_fill_gradient(low = "black", high = "red", limits = c(mmin,MMin)) + th
  scale_fill_gradient(low = "black", high = "red") + th

pot <- 
  ggplot(mr[mr$type == "ot",], 
  aes(x,y, group = factor(type), fill = value)) +
  geom_raster() +
  scale_y_reverse() + 
  facet_grid(ch+F0 ~ z, labeller = label_both) +
  #scale_fill_gradient(low = "black", high = "blue", limits = c(mmot,MMot)) + th
  scale_fill_gradient(low = "black", high = "blue") + th
#

lay <- rbind(matrix(1,11,1), 
      matrix(2,8,1), 
      matrix(3,4,1))



outName <- paste("~/Desktop/zTest",
           paste0("C",paste(loc[k,], collapse = "_"), sep = "_", collapse = "_") )
png(paste0(outName, ".png"), 
    height = 23*1.5, width = 11*1.5, units = "in", res = 300)
print(grid.arrange(pex, pin, pot, layout_matrix = lay))
dev.off()
}

