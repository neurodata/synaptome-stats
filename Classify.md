---
title: "Classification on RORB synapses"
author: "Jesse Leigh Patsolic"
output: 
  html_document:
    keep_md: true
---

<!--
### ### INITIAL COMMENTS HERE ###
###
### Jesse Leigh Patsolic 
### 2018 <Jesse.L.Patsolic@alumni.wfu.edu>
### S.D.G 
#

-->

<style type="text/css">
.table {
    width: 30%;
}
tr:hover {background-color:#f5f5f5;}
</style>






```r
#gs_auth(new_user = TRUE)
gs_ls()

tmp <- gs_title("MNSite3Synaptograms.csv")
fdat <- as.data.table(gs_read_csv(ss = tmp, ws = 1, skip = 1))
rm(tmp)

fdat$cx <- 
  (fdat$maxX - fdat$minX)/2 + fdat$minX

fdat$cy <- 
  (fdat$maxY - fdat$minY)/2 + fdat$minY

fdat$cz <- 
  (fdat$maxZ - fdat$minZ)/2 + fdat$minZ

fdat$atx <- round((fdat$cx * 3) / 96)
fdat$aty <- round((fdat$cy * 3) / 96)
fdat$atz <- round(fdat$cz)

fdat$nmx <- (fdat$cx / 3)
fdat$nmy <- (fdat$cy / 3)
fdat$nmz <- fdat$cz/50

dx <- 438
dy <- 460
dz <- 49

ind <- 
  (dat$atx - 5) >= 0  &
  (dat$aty - 5) >= 0  &
  (dat$atz - 5) >= 0  &
  (dat$atx + 6) <= dx &
  (dat$aty + 6) <= dy &
  (dat$atz + 6) <= dz

fdat$buff <- as.numeric(ind)
```

# Merge data with labels by ID


```r
avg <- read.csv("rorb_avg_at_123ids.csv")
names(avg)[1] <- "DAPI1"
loc <- read.csv("rorb_avg_at_123Locations.csv")
names(loc)[1] <- 'x'
loc$id <- loc$id + 1 #Adjust id

dat <- cbind(loc, avg)

Y <- fdat[, .(GABA, id)]
Y[, ':=' (id = id + 1, Y = GABA, GABA = NULL)]

subdat <- dat[dat$id %in% Y$id, ]

m1 <- merge(Y, dat, by = 'id')

X <- m1[, -c(1,3:5)]
X[, Y := as.factor(Y)]
```

# Head


```r
ind <- c(which(X$Y == 0)[1:3], which(X$Y == 1)[1:3])
X[ind, -c(2,3,4,12)]
```

```
##    Y     GABA      GAD2 Gephyrin    GluN1      MBP    PSD95 synapsin
## 1: 0 1.659654  6.424493 3.872276 8.947408 0.811420 5.511645 5.819684
## 2: 0 0.709992  7.243426 2.845229 6.234410 1.444027 4.184072 1.386927
## 3: 0 1.482344  3.693464 2.167543 3.223140 1.018032 2.694215 0.776860
## 4: 1 5.376409 10.532682 2.771600 7.486852 0.948911 4.314801 6.265965
## 5: 1 8.032307 15.578512 4.747558 8.625094 1.053343 4.529677 5.246431
## 6: 1 3.537190  7.006762 3.507889 7.169797 1.216379 3.503381 1.344853
##      VGlut1
## 1: 5.858002
## 2: 7.084899
## 3: 3.212622
## 4: 7.865515
## 5: 9.308039
## 6: 6.761082
```


# Classification

## caret using LOOCV for various methods


```r
tc <- trainControl(method = "LOOCV")
ldaFit <- train(Y ~ ., data = X, trControl = tc, method = "lda")
bayesglmFit  <- train(Y ~ ., data = X, trControl = tc, method = "bayesglm")

clda <- sum(predict(ldaFit) != X$Y)/length(X$Y)
cbayes <- sum(predict(bayesglmFit) != X$Y)/length(X$Y)
```

## RF


```r
rfMod <- randomForest(Y ~ ., data = X, importance = TRUE)
rfPred <- predict(rfMod)
rfL <- sum(rfPred != X$Y)/length(X$Y)
```

## MclustDA

```r
daMod <- MclustDA(X[, -1], X$Y)
daPred <- predict(daMod)$classification
daL <- sum(daPred != X$Y)/length(X$Y)
```

## LDA


```r
ldaMod <- lda(Y ~ ., data = X, CV = FALSE)
ldaPred <- predict(ldaMod)$class
ldaL <- sum(ldaPred != X$Y)/length(X$Y)
```

## LDA with LOOCV


```r
ldaModcv <- lda(Y ~ ., data = X, CV = TRUE)
ldacvL <- sum(ldaModcv$class != X$Y)/length(X$Y)
```

| method | error |
|------------|-------|
| caret LDA w/LOOCV | 0.0895334|
| caret Bayes GLM w/LOOCV | 0.0819672|
| RF | 0.0970996|
| MclustDA | 0.074401|
| LDA | 0.0895334|
| LDA w/LOOCV | 0.0983607|

# Spatial point pattern

The Ripley's K function I found only takes 2d spatial data ...
Will keep looking.


<!--
#   Time:
##  Working status:
### Comments:
####Soli Deo Gloria
--> 

