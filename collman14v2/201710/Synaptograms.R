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
require(rhdf5)
require(reshape2)
require(gplots)
require(raster)
require(gridExtra)
require(data.table)

input <- "collman14v2_tightAnnotationCubes20171101T2000_C5synaptograms.h5"
Locinput <- "collman14v2_tightAnnotationCubes20171101T2000_C5synaptograms.h5"
input <- "collman14v2_tightAnnotationCubes20171101T2000.csv.h5"

input <- "meda_plots_gaussian/collman14v2_gaussianSynaptograms.csv.h5"

input <- "collman14v2_fullCubes_20171101T1630.csv.h5"

h5ls(input)
name <- "collman14v2_cubes"
loc <- h5read(input, name = 'Locations')
dat <- h5read(input, name = name)

chan <- h5read(input, name = 'Channels')
H5close()

set.seed(1030)
s1 <- sample(dim(dat)[4], 10)
dat <- dat[,,,s1,]

type <- c('em', "ot", "ot", "ot", 
          'in', 'in', 'in','in',
          'ex', 'ex', 'ex', 'ex')
          
ctype <- data.frame(chan, type=type)

empal <- colorpanel(255,"black", "white")
inpal <- colorpanel(255,"black", "red")
expal <- colorpanel(255, "black", "green")
otpal <- colorpanel(255, "black", "blue")


dimnames(dat) <- list(NULL, NULL, NULL, 1:dim(dat)[4], chan)
if(FALSE){
  G <- readRDS("maskGaussian_mu0_sigma2123_2123_5-455.rds")
  plot(raster(G[,,6]))
  
  for(j in 1:dim(dat)[5]){
    for(i in 1:dim(dat)[4]){
    dat[,,,i,j] <- G * dat[,,,i,j]
    }
  }
}


rr <- foreach(l = 1:dim(dat)[4]) %do% {
  r <- dat[,,,l,]
  mr <- melt(r)
  colnames(mr) <- c("x", "y", "z", "ch", "value")
  mr$z <- mr$z - (max(mr$z) +1)/2
  ch <- mr$ch
  mr$type <- ctype$type[ch]
  mr
}


th <- theme(axis.text = element_blank(), axis.ticks = element_blank(),
            axis.title.y = element_blank(), axis.title.x = element_blank(),
            legend.position="bottom", legend.key.size = unit(2,'lines'),
            panel.spacing = unit(0, "lines"))


#pg <- list()
#for(i in 1:length(rr)){
#  pdat <- rr[[i]]
#  pg[[i]] <- 
#    ggplot(pdat, 
#    aes(x,y, group = ch, fill = value)) +
#    geom_raster() + 
#    scale_y_reverse() + 
#    facet_grid(ch ~ z, labeller = label_both) +
#    #scale_fill_gradient(low = "black", high = "green") + th
#    scale_fill_gradient(low = "black", high = "white") + th
#} 
#
#
##pdf("collman14v2_synaptograms.pdf", width = 12, height = 12)
#for(i in 1:length(pg)){
#  plot(pg[[i]])
#  #Sys.sleep(3)
#}
##dev.off()
#
#
#for(k in 1:length(rr)){
#  mr <- rr[[k]]
#  pex <- 
#  ggplot(mr[mr$type == "ex",], 
#  mr
#}

mmem <- min(sapply(rr, function(am){ min(am[am$type == "em",]$value) }))
MMem <- max(sapply(rr, function(am){ max(am[am$type == "em",]$value) }))
mmex <- min(sapply(rr, function(am){ min(am[am$type == "ex",]$value) }))
MMex <- max(sapply(rr, function(am){ max(am[am$type == "ex",]$value) }))
mmin <- min(sapply(rr, function(am){ min(am[am$type == "in",]$value) }))
MMin <- max(sapply(rr, function(am){ max(am[am$type == "in",]$value) }))
mmot <- min(sapply(rr, function(am){ min(am[am$type == "ot",]$value) }))
MMot <- max(sapply(rr, function(am){ max(am[am$type == "ot",]$value) }))


P <- list()
for(k in 1:length(rr)){
  mr <- rr[[k]]

  pem <- 
  ggplot(mr[mr$type == "em",], 
    aes(x,y, group = factor(type), fill = value)) +
    geom_raster() + 
    scale_y_reverse() + 
    facet_grid(ch ~ z, labeller = label_both) +
    #scale_fill_gradient(low = "black", high = "white", limits = c(mmem, MMem)) + th + 
    scale_fill_gradient(low = "black", high = "white") + th + 
    theme(legend.position = "none")


  pex <- 
  ggplot(mr[mr$type == "ex",], 
    aes(x,y, group = factor(type), fill = value)) +
    geom_raster() + 
    scale_y_reverse() + 
    facet_grid(ch ~ z, labeller = label_both) +
    #scale_fill_gradient(low = "black", high = "green", limits = c(mmex, MMex)) + th
    scale_fill_gradient(low = "black", high = "green") + th
  
  pin <- 
  ggplot(mr[mr$type == "in",], 
    aes(x,y, group = factor(type), fill = value)) +
    geom_raster() + 
    scale_y_reverse() + 
    facet_grid(ch ~ z, labeller = label_both) +
    #scale_fill_gradient(low = "black", high = "red", limits = c(mmin,MMin)) + th
    scale_fill_gradient(low = "black", high = "red") + th
  
  pot <- 
    ggplot(mr[mr$type == "ot",], 
    aes(x,y, group = factor(type), fill = value)) +
    geom_raster() +
    scale_y_reverse() + 
    facet_grid(ch ~ z, labeller = label_both) +
    #scale_fill_gradient(low = "black", high = "blue", limits = c(mmot,MMot)) + th
    scale_fill_gradient(low = "black", high = "blue") + th
  #
  
  #lay <- rbind(
  #      matrix(1,table(type)['em'],1),
  #      matrix(2,table(type)['ex'],1), 
  #      matrix(3,table(type)['in'],1), 
  #      matrix(4,table(type)['ot'],1))

  lay <- rbind(
        matrix(1,table(type)['ot'],1),
        matrix(2,table(type)['in'],1), 
        matrix(3,table(type)['ex'],1), 
        matrix(4,table(type)['em'],1))
  
  
  #P[[k]] <- grid.arrange(pem, pex, pin, pot, layout_matrix = lay)
  P[[k]] <- grid.arrange(pot, pin, pex, pem, layout_matrix = lay)
  #print(grid.arrange(pem, pex, pin, pot, layout_matrix = lay))

  if(FALSE){
    cname <- 
    paste0("meda_plots_gaussian/collman14v2_top5_GaussianSynaptogram", sprintf("_x%d_y%d_z%d.png", loc[k,1], loc[k,2], loc[k,3]))
    
    b = 1080
    w = b
    h = 1.25*b
  
    png(cname, width = w, height=h)
    plot(P[[k]])
    dev.off()
  }
}







