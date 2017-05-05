



suppressPackageStartupMessages(require(rhdf5))
suppressPackageStartupMessages(require(meda))
require(gridExtra)
require(foreach)

fname <- "Ex12R75vol2_synaptograms11cubes3NN.h5"
fname <- "~/Desktop/kristina15rscaled317samp5000_synaptograms11cubes3NN.h5"

h5ls(fname) 
#ex12r75v2<- dat <- h5read(fname, name = "Ex12R75")
dat <- h5read(fname, name = "kristina15")
datO <- dat
loc <- h5read(fname, name = "Locations")
cnames <- h5read(fname, name = "Channels")
H5close()

type <- c("Inhibitory", "Other", "Other", "Other", "Other", "Other", 
"Other", "Excitatory", "Excitatory", "Excitatory", "Inhibitory", 
"Excitatory", "Excitatory", "Other", "Inhibitory", "Excitatory", 
"Excitatory", "Excitatory", "Inhibitory", "Excitatory", "Inhibitory"
)

ctype <- data.frame(cnames, type[c(1,1,1,1:21)])
ord <- order(type)

inpal <- colorpanel(255,"black", "red")
expal <- colorpanel(255, "black", "green")
otpal <- colorpanel(255, "black", "blue")


### Normalize across channels
dim(dat)
### dat[x,y,z,synapse, channel]
#datN <- 
for(i in 1:dim(dat)[4]){
  for(j in 1:dim(dat)[5]){
    a1 <- (dat[,,,i,j])
    if(!all(a1 == 0)){
      a1 <- (a1 - min(a1))/(max(a1) - min(a1))
      dat[,,,i,j] <- a1
    }
  }
}

#dat <- A #REMOVE ME
rr <- foreach(l = 1:dim(dat)[4]) %do% {
  r <- dat[,,,l,]
  mr <- melt(r)
  names(mr) <- c("x", "y", "z","chN", "value")
  mr$z <- mr$z - (max(mr$z) +1)/2
  mr$ch <- cnames[mr$chN]
  mr$type <- ctype$type[mr$chN]
  mr
}


names(rr) <- 1:length(rr)

th <- theme(axis.text = element_blank(), axis.ticks = element_blank(),
            axis.title.y = element_blank(), axis.title.x = element_blank(),
            panel.spacing = unit(0, "lines"))

for(k in 1:length(rr)){
  mr <- rr[[k]]
pex <- 
ggplot(mr[mr$type == "Excitatory",], 
  aes(x,y, group = factor(type), fill = value)) +
  geom_raster() + 
  facet_grid(ch ~ z, labeller = label_both) +
  scale_fill_gradient(low = "black", high = "green") + th

pin <- 
ggplot(mr[mr$type == "Inhibitory",], 
  aes(x,y, group = factor(type), fill = value)) +
  geom_raster() + 
  facet_grid(ch ~ z, labeller = label_both) +
  scale_fill_gradient(low = "black", high = "red") + th

pot <- 
  ggplot(mr[mr$type == "Other",], 
  aes(x,y, group = factor(type), fill = value)) +
  geom_raster() +
  facet_grid(ch ~ z, labeller = label_both) +
  scale_fill_gradient(low = "black", high = "blue") + th

lay <- rbind(matrix(1,9,21), 
      matrix(2,5,21), 
      matrix(3,7,21))

clusterN <- hmcO$dat$Get("name", filterFun = isLeaf)
cseq <- rep(1:8,each =3)
outName <- paste("~/Desktop/Ex12R75vol2",paste0("C", clusterN[cseq[k]]), paste(loc[k,], collapse = "_"),
              sep = "_", collapse = "_") 
pdf(paste0(outName,".pdf"), height = 21*1.5, width = 11*1.5)
#png(paste0(outName, ".png"), 
#    height = 21*1.5, width = 11*1.5, units = "in", res = 300)
print(grid.arrange(pex, pin, pot, layout_matrix = lay))
dev.off()
}

