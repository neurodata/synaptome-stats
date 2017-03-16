#!/cm/shared/apps/R/3.2.2/bin/Rscript
options(echo = TRUE)
args <- commandArgs(TRUE)

sampN <- 
  if(length(args) > 0){ 
    seq(1e3, 1e4, by = 1e3)[as.integer(args[1])]
  } else {
    stop()
  }

runplots <- ifelse(length(args) > 1, TRUE, FALSE)

fileBase <- as.character(sampN)
outDir <- paste0("figures_samp", sampN)
if(!file.exists(outDir)){
  dir.create(outDir)
}

source("https://raw.githubusercontent.com/neurodata/discriminability/master/Code/FlashRupdated/functions/reliability.R")
require(meda)
require(raster)
base <- Sys.getenv("NEURODATA")

if(!runplots){
  load(paste0(base,'/synaptome-stats/Code/cleanDataWithAttributes.RData'))
  
  set.seed(sampN*2^13)
  dat <- data01[sample(nrow(data01), sampN), 1:24, with = FALSE]
  
  colnames(dat) <- gsub("_F0", "", colnames(dat))
  
  ccol3 <- c("#197300", "#197300", "#197300", "#197300", "#197300", "#197300", 
    "#197300", "#197300", "#197300", "#197300", "#197300", "#cc0000", 
    "#cc0000", "#cc0000", "#cc0000", "#cc0000", "#cc0000", "#0000cd", 
    "#0000cd", "#0000cd", "#0000cd", "#0000cd", "#0000cd", "#0000cd"
    )
  #transformData <- function(x, type = c("1e3"), t98 = FALSE, ...) {
  datRaw <- dat
  datLog <- transformData(dat, type = c("log10"))$log10
  dat01e3 <- transformData(dat, type = c("1e3"))$d01e3
  
  N <- lapply(list(datRaw, datLog, dat01e3), nrow)
  
  Lraw <- p.hmc(datRaw, maxDepth = 5, maxDim = 8, modelNames = c("VVV"))
  Llog <- p.hmc(datLog, maxDepth = 5, maxDim = 8, modelNames = c("VVV"))
  L01 <- p.hmc(dat01e3, maxDepth = 5, maxDim = 8, modelNames = c("VVV"))
  
  labRaw <- Lraw$labels$col
  labLog <- Llog$labels$col
  lab01e3 <- L01$labels$col
  
  rawLab <- labRaw[order(labRaw)] + 1e5
  logLab <- labLog[order(labLog)] + 1e5
  e3Lab <-  lab01e3[order(lab01e3)] + 1e5
  
  p1 <- order(labRaw)
  p2 <- order(labLog)
  p3 <- order(lab01e3)
  
  distRaw <- as.matrix(dist(datRaw))[p1,p1]
  distLog <- as.matrix(dist(datLog))[p2,p2]
  dist01e3 <- as.matrix(dist(dat01e3))[p3,p3]
  
  rank_distRaw <- apply(distRaw, 1, rank, ties.method = "average")-1
  rank_distLog <- apply(distLog, 1, rank, ties.method = "average")-1
  rank_dist01e3 <- apply(dist01e3, 1, rank, ties.method = "average")-1
  
  relRaw <- rdf_large_s(rank_distRaw, as.character(rawLab))
  relLog <- rdf_large_s(rank_distLog, as.character(logLab))
  rel01e3 <- rdf_large_s(rank_dist01e3, as.character(e3Lab))
  
  discRaw <- mnr(relRaw, rawLab, remove_outliers = FALSE, output = TRUE)
  discLog <- mnr(relLog, logLab, remove_outliers = FALSE, output = TRUE)
  disc01e3 <- mnr(rel01e3, e3Lab, remove_outliers = FALSE, output = TRUE)
  
  ci <- mapply(function(x, n){ x + 2*sqrt(5/n) *c(-1,1) },
               list(discRaw, discLog, disc01e3), N, SIMPLIFY = FALSE)
  
  save.image(file=paste0("outData_samp", sampN, ".RData"))
} else {
  load(paste0("outData_samp", sampN, ".RData"))

  ci <- mapply(function(x, n){ x + 2*sqrt(5/n) *c(-1,1) },
               list(discRaw, discLog, disc01e3), N, SIMPLIFY = FALSE)

  pdf(paste0(outDir,"/dend", sampN, ".pdf"), height = 4, width = 12)
  par(mfrow = c(1,3))
  p.dend(Lraw); title("Raw")
  p.dend(Llog); title("Log10")
  p.dend(L01); title("01e3")
  dev.off()
  
  pdf(paste0(outDir, "/dist", sampN, ".pdf"), height = 6, width = 18)
  par(mfrow = c(1,3))
  plot(raster(distRaw), col = viridis(255)[1:215], main = "Distance matrix sorted by label: Raw")
  plot(raster(distLog), col = viridis(255)[1:215], main = "Distance matrix sorted by label: Log10")
  plot(raster(dist01e3), col = viridis(255)[1:215], main = "Distance matrix sorted by label: 01e3")
  dev.off()
  
  pdf(paste0(outDir, "/rankdist", sampN, ".pdf"), height = 6, width = 18)
  par(mfrow = c(1,3))
  plot(raster(log(rank_distRaw+1)), col = rev(gray.colors(255)), main = "log rank distances: Raw")
  plot(raster(log(rank_distLog+1)), col = rev(gray.colors(255)), main = "log rank distances: Log10")
  plot(raster(log(rank_dist01e3+1)), col = rev(gray.colors(255)), main = "log rank distances: 01e3")
  dev.off()
  
  
  pdf(paste0(outDir, "/relia", sampN, ".pdf"), height = 6, width = 18)
  
  par(mfrow = c(1,3))
  t1 <- "Reliability: Raw (blue chance, red MNR, green 95% CI)"
  hist(relRaw, main = t1, breaks = "Scott"); abline(v = c(ci[[1]], discRaw), lwd = 2, col = c(3,3,2))
  abline(v = 0.5, col = 'blue', lwd = 3, lty = 2)
  
  t2 <- "Reliability: Log10 (blue chance, red MNR, green 95% CI)"
  hist(relLog, main = t2, breaks = "Scott"); abline(v = c(ci[[2]], discLog), lwd = 2, col = c(3,3,2))
  abline(v = 0.5, col = 'blue', lwd = 3, lty = 2)
  
  t3 <- "Reliability: 01e3 (blue chance, red MNR, green 95% CI)"
  hist(rel01e3, main =t3, breaks = "Scott"); abline(v = c(ci[[3]], disc01e3), lwd = 2, col = c(3,3,2))
  abline(v = 0.5, col = 'blue', lwd = 3, lty = 2)
  
  dev.off()
}
