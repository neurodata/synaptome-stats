suppressPackageStartupMessages(require(rhdf5))
suppressPackageStartupMessages(require(meda))
require(gridExtra)
require(foreach)

fname <- "../Python/kristina15rscaled317samp5000_synaptograms11cubes5NN_test7.h5"

hmcO <- readRDS("run10/data/outputs/K15F0_rscale317samp5000/hmc.rds")
clusterN <- hmcO$dat$Get("name", filterFun = isLeaf)

h5ls(fname) 
cnames <- h5read(fname, name = "Channels")
abcnames <- abbreviate(gsub("_", "", cnames), minlength = 3)
dat <- h5read(fname, name = "kristina15")
dimnames(dat) <- list(1:40, cnames, NULL, NULL, NULL)
datO <- dat
loc <- h5read(fname, name = "Locations")
trueLoc <- as.character(clusterN[read.csv("K15F0_rscaled317samp5e3_hmc5NN.csv", header = TRUE)[,4]])
#cseq <- rep(1:10,each =5)
#clusterN <- hmcO$dat$Get("name", filterFun = isLeaf)
#trueLoc <- cbind(trueLoc, cseq)
#trueCl <- trueLoc[trueLoc$z >=5, 'cseq']


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

## Calculate feature F0
F0 <- 
  foreach(j = 1:dim(dat)[1]) %:%
  foreach(i = 1:dim(dat)[2]) %do% {
      f0 <- sum(dat[j,i,,,])
}

names(F0) <- 1:length(F0)
for(i in 1:length(F0)){
  names(F0[[i]]) <- cnames
}


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

for(k in 1:length(rr)){
  mr <- rr[[k]]
  pex <- 
  ggplot(mr[mr$type == "ex",], 
    aes(x,y, group = factor(type), fill = value)) +
    geom_raster() + 
    scale_y_reverse() + 
    facet_grid(ch + F0 ~ z, labeller = label_both) +
    scale_fill_gradient(low = "black", high = "green", limits = c(mmex,MMex)) + th
  
  pin <- 
  ggplot(mr[mr$type == "in",], 
    aes(x,y, group = factor(type), fill = value)) +
    geom_raster() + 
    scale_y_reverse() + 
    facet_grid(ch + F0 ~ z, labeller = label_both) +
    scale_fill_gradient(low = "black", high = "red", limits = c(mmin,MMin)) + th
  
  pot <- 
    ggplot(mr[mr$type == "ot",], 
    aes(x,y, group = factor(type), fill = value)) +
    geom_raster() +
    scale_y_reverse() + 
    facet_grid(ch + F0 ~ z, labeller = label_both) +
    scale_fill_gradient(low = "black", high = "blue", limits = c(mmot,MMot)) + th
  #
  
  lay <- rbind(matrix(1,table(type)['ex'],1), 
        matrix(2,table(type)['in'],1), 
        matrix(3,table(type)['ot'],1))
  
  
  
  outName <- paste("~/Desktop/zScale",
             paste0("C", trueLoc[k]), paste(loc[k,], collapse = "_"), sep = "_", collapse = "_") 
  png(paste0(outName, ".png"), 
      height = 23*1.5, width = 11*1.5, units = "in", res = 300)
  print(grid.arrange(pex, pin, pot, layout_matrix = lay))
  dev.off()
}


cisbaseZ <- "http://cis.jhu.edu/~jesse/fotd/20170509/synaptograms/k15syn_log10scale"
cisbase0 <- "http://cis.jhu.edu/~jesse/fotd/20170509/synaptograms/k15syn_log10scale"
ndviz <- "http://viz.neurodata.io/project/kristina15_SynDiv_JLP/xy/0/"


out <- foreach(k = 1:length(rr), .combine = 'rbind')%do%{
 outNameZ <- paste(cisbaseZ,
            paste0("C", trueLoc[k]), paste(loc[k,], collapse = "_"), sep = "_", collapse = "_") 
 
 outName0 <- paste(cisbase0,
            paste0("C", trueLoc[k]), paste(loc[k,], collapse = "_"), sep = "_", collapse = "_") 

 cl <- paste0("C", trueLoc[k], "_", paste(loc[k,], collapse = "_"))

 cluZ <-  paste0(outNameZ, ".png")
 clu0 <-  paste0(outName0, ".png")

 nd <- paste0(ndviz, paste(loc[k,1], loc[k,2], loc[k,3], sep = "/"))

 sprintf("[rscale%s](%s) \t\t| [01scale%s](%s) \t\t| [viz%s](%s/)  ",cl, cluZ, cl, clu0, cl, nd)
}

write.table(out, file = "~/Desktop/links.txt", quote = FALSE, row.names = FALSE)




