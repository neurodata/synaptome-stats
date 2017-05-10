suppressPackageStartupMessages(require(rhdf5))
suppressPackageStartupMessages(require(meda))
require(gridExtra)
require(foreach)

fname <- "run10/kristina15rscaled317samp5000_synaptograms11cubes5NN_test3.h5"

hmcO <- readRDS("run10/data/outputs/K15F0_rscale317samp5000/hmc.rds")
clusterN <- hmcO$dat$Get("name", filterFun = isLeaf)

h5ls(fname) 
dat <- h5read(fname, name = "kristina15")
datO <- dat
loc <- h5read(fname, name = "Locations")
trueLoc <- as.character(clusterN[read.csv("K15F0_rscaled317samp5e3_hmc5NN.csv", header = TRUE)[,4]])
#cseq <- rep(1:10,each =5)
#clusterN <- hmcO$dat$Get("name", filterFun = isLeaf)
#trueLoc <- cbind(trueLoc, cseq)
#trueCl <- trueLoc[trueLoc$z >=5, 'cseq']


cnames <- h5read(fname, name = "Channels")
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
      a1 <- scale(a1, center = TRUE, scale = TRUE)
      #a2 <- (a1 - min(a1))/(max(a1) - min(a1))
      dat[,j,,,] <- a1
    }
  }


rr <- foreach(l = 1:dim(dat)[1]) %do% {
  r <- dat[l,,,,]
  mr <- melt(r)
  names(mr) <- c("chN", "x", "y", "z", "value")
  mr$z <- mr$z - (max(mr$z) +1)/2
  mr$ch <- cnames[mr$chN]
  mr$type <- ctype$type[mr$chN]
  mr
}


names(rr) <- 1:length(rr)

th <- theme(axis.text = element_blank(), axis.ticks = element_blank(),
            axis.title.y = element_blank(), axis.title.x = element_blank(),
            panel.spacing = unit(0, "lines"))

#registerDoMC(6)
#foreach(k = 1:length(rr))%dopar%{
for(k in 1:length(rr)){
  mr <- rr[[k]]
#mm <- 0
#MM <- 1
mm <- min(mr[mr$type == 'ex',]$value)
MM <- max(mr[mr$type == 'ex',]$value)
pex <- 
ggplot(mr[mr$type == "ex",], 
  aes(x,y, group = factor(type), fill = value)) +
  geom_raster() + 
  scale_y_reverse() + 
  facet_grid(ch ~ z, labeller = label_both) +
  scale_fill_gradient(low = "black", high = "green", limits = c(mm,MM)) + th

mm <- min(mr[mr$type == 'in',]$value)
MM <- max(mr[mr$type == 'in',]$value)
pin <- 
ggplot(mr[mr$type == "in",], 
  aes(x,y, group = factor(type), fill = value)) +
  geom_raster() + 
  scale_y_reverse() + 
  facet_grid(ch ~ z, labeller = label_both) +
  scale_fill_gradient(low = "black", high = "red", limits = c(mm,MM)) + th

mm <- min(mr[mr$type == 'ot',]$value)
MM <- max(mr[mr$type == 'ot',]$value)
pot <- 
  ggplot(mr[mr$type == "ot",], 
  aes(x,y, group = factor(type), fill = value)) +
  geom_raster() +
  scale_y_reverse() + 
  facet_grid(ch ~ z, labeller = label_both) +
  scale_fill_gradient(low = "black", high = "blue", limits = c(mm,MM)) + th
#

lay <- rbind(matrix(1,11,23), 
      matrix(2,8,23), 
      matrix(3,4,23))



outName <- paste("~/Desktop/k15syn_rscale",
           paste0("C", trueLoc[k]), paste(loc[k,], collapse = "_"), sep = "_", collapse = "_") 
png(paste0(outName, ".png"), 
    height = 23*1.5, width = 11*1.5, units = "in", res = 300)
print(grid.arrange(pex, pin, pot, layout_matrix = lay))
dev.off()
}


cisbaseZ <- "http://cis.jhu.edu/~jesse/fotd/20170509/synaptograms/k15syn_rscale"
cisbase0 <- "http://cis.jhu.edu/~jesse/fotd/20170509/synaptograms/k15syn_01scale"
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




