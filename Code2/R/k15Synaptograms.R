suppressPackageStartupMessages(require(rhdf5))
suppressPackageStartupMessages(require(meda))
require(gridExtra)
require(foreach)

fname <- "kristina15rscaled317samp5000_synaptograms11cubes5NN_test3.h5"

h5ls(fname) 
dat <- h5read(fname, name = "kristina15")
datO <- dat
loc <- h5read(fname, name = "Locations")
trueLoc <- read.csv("K15F0_rscaled317samp5e3_hmc5NN.csv", header = TRUE)
cseq <- rep(1:10,each =5)
clusterN <- hmcO$dat$Get("name", filterFun = isLeaf)
trueLoc <- cbind(trueLoc, cseq)
trueCl <- trueLoc[trueLoc$z >=5, 'cseq']


cnames <- h5read(fname, name = "Channels")
H5close()

hmcO <- readRDS("data/outputs/K15F0_rscale317samp5000/hmc.rds")

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


clusterN <- hmcO$dat$Get("name", filterFun = isLeaf)
cseq <- trueCl
outName <- paste("~/Desktop/k15syn_rscale",paste0("C", clusterN[cseq[k]]), paste(loc[k,], collapse = "_"), sep = "_", collapse = "_") 
#pdf(paste0(outName,".pdf"), height = 23*1.5, width = 11*1.5)
png(paste0(outName, ".png"), 
    height = 23*1.5, width = 11*1.5, units = "in", res = 300)
print(grid.arrange(pex, pin, pot, layout_matrix = lay))
dev.off()
}

x <- loc[,1]
y <- loc[,2]
z <- loc[,3]
fgh <- sprintf("http://viz.neurodata.io/project/kristina15_SynDiv_JLP/xy/0/%d/%d/%d/",
        x,y,z)

write.table(fgh, "~/Desktop/links.csv", row.names = FALSE, quote = FALSE)

cis <-  sprintf("http://viz.neurodata.io/project/kristina15_SynDiv_JLP/xy/0/%d/%d/%d/",

cisbase <- "http://cis.jhu.edu/~jesse/fotd/20170509/synaptograms/k15syn_rscale"

lllfoo <- paste0("C",v,"_",apply(loc, 1, paste, collapse = "_"))
outName <- sapply(1:length(rr), function(k){ 
  paste(cisbase,paste0("C", clusterN[cseq[k]]), paste(loc[k,], collapse = "_"), sep = "_", collapse = "_")})

outmlp <- mapply(function(x,y,z){
  paste0("[", x, "]", "(",y,") | [", x, "](", z,")")
  }, lllfoo, outName, fgh)
