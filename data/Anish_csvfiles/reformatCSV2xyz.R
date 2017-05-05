


basedir <- getwd()
dir0 <- grep("^Ex", dir(), value = TRUE)
dir0

for(i in dir0){
 exdir <- i 
 csvf <- dir(exdir, full.names=TRUE)
 for(j in csvf){
   r <- read.csv(j, header = FALSE)
   n <- paste0("XYZ/", strsplit(j, "\\.")[[1]][1], "xyz.csv")
   out <- r[, c(2,1,3)]
   colnames(out) <- c("col", "row", "z")
   write.csv(out, file = n, row.names = FALSE)
 }
}
