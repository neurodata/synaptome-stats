# Generate ndviz links for synapse detections
Jesse Leigh Patsolic  





```r
suppressPackageStartupMessages(require(rhdf5))
suppressPackageStartupMessages(require(meda))

#h5ls("../../data/Anish_csvfiles/XYZ/AnishLocationsFull.h5")
dat <- data.table(h5read("../../../data/Anish_csvfiles/XYZ/AnishLocationsFull.h5", name = "Locations/locs"))
exlev <- h5read("../../../data/Anish_csvfiles/XYZ/AnishLocationsFull.h5", name = "Levels/Ex")
qlev <-  h5read("../../../data/Anish_csvfiles/XYZ/AnishLocationsFull.h5", name = "Levels/Q")

dat$Ex <- as.factor(dat$Ex); levels(dat$Ex) <- exlev 
dat$Q <- as.factor(dat$Q); levels(dat$Q) <- qlev 

H5close()

subd <- dat[Ex == "Ex2R18C1" & 
            row < 100 & 
            row > 50 &
            z < 13.5 & 
            z > 12.5 & 
            col < 700]

xyz <- subd[, `:=`(col = round(col),
            row = round(row), 
            z = round(z))][, .(col, row, z)]

l <- mapply(function(x,y,z){
  sprintf("http://viz.neurodata.io/project/exqc/xy/0/%d/%d/%d/",x,y,z)
  }, xyz$col, xyz$row, xyz$z)

#lapply(l, function(x){ system(paste("open ", x))})
print(l)
```

```
##  [1] "http://viz.neurodata.io/project/exqc/xy/0/604/65/13/"
##  [2] "http://viz.neurodata.io/project/exqc/xy/0/602/73/13/"
##  [3] "http://viz.neurodata.io/project/exqc/xy/0/560/59/13/"
##  [4] "http://viz.neurodata.io/project/exqc/xy/0/664/59/13/"
##  [5] "http://viz.neurodata.io/project/exqc/xy/0/600/65/13/"
##  [6] "http://viz.neurodata.io/project/exqc/xy/0/580/70/13/"
##  [7] "http://viz.neurodata.io/project/exqc/xy/0/590/77/13/"
##  [8] "http://viz.neurodata.io/project/exqc/xy/0/544/84/13/"
##  [9] "http://viz.neurodata.io/project/exqc/xy/0/586/98/13/"
## [10] "http://viz.neurodata.io/project/exqc/xy/0/590/77/13/"
## [11] "http://viz.neurodata.io/project/exqc/xy/0/563/51/13/"
## [12] "http://viz.neurodata.io/project/exqc/xy/0/664/60/13/"
```




