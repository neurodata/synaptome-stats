library(RMySQL)
require(data.table)
require(raster)
require(meda)
## Set up the connection to the database.
## REMEMBER to disconnect when finished.
#con1  <- dbConnect(MySQL(), user="root", dbname="ArchMath",
#                   password="2113StPaul");dbSendQuery(con1, "SET NAMES 'UTF8';")

con0 <- file("mysqlPassword.txt", "r")
pswd <- readLines(con0, n = 1)
con1  <- dbConnect(MySQL(), user="root", dbname="synapses",
                   password=pswd);dbSendQuery(con1, "SET NAMES 'UTF8';")

chan <- c('synapsinR_7thA','synapsinGP_5thA','VGluT1_3rdA','VGluT1_8thA',
					'VGluT2_2ndA', 'psd_8thA', 'GluR2_2ndA', 'NMDAR1_6thA',
					'NR2B_9thA','NOS_9thA', 'Synpod_3rdA', 'GAD_6thA','VGAT_5thA', 
					'PV_1stA','gephyrin_1stA', 'GABARa1_4thA', 'GABABR1_3rdA',
					'VGluT3_1stA', 'CR1_2ndA','5HT1A_6th', 'TH_5thA','VAChT_4thA',	
          'tubulin_8thA', 'dapi_1stA')

Syncol3 <- c("#197300","#cc0000","#0000cd")

ccol <- Syncol3[c(rep(1,11), rep(2,6), rep(3,7))]

q1 <- c("SELECT * FROM kristina15")

dat0 <- data.table(dbGetQuery(con1, q1))
dbDisconnect(con1)
dat0 <- dat0[, c("x", "y", "z", chan), with = FALSE]
setkey(dat0, x,y,z)

xyz <- dat0[, .(x,y,z)]
dat <- dat0[, !(x:z)]


if(!file.exists("sample_locations_1e4.csv")){
  print("Selecting sample locations")
  
  set.seed(317)
  s1 <- sample(nrow(dat), 1e4)
  write.csv(xyz[s1], file = "sample_locations_1e4.csv", row.names = FALSE)

  set.seed(1030)
  s2 <- sample(nrow(dat), 1e3)
  write.csv(xyz[s2], file = "sample_locations_1e3.csv", row.names = FALSE)
} else {
  samp10k <- read.csv("sample_locations_1e4.csv", header=TRUE)
  samp1k  <- read.csv("sample_locations_1e3.csv", header=TRUE)
}


d1ko  <- dat0[samp1k, !(x:z)]
d10ko <- dat0[samp10k, !(x:z)]

d1k  <- d1ko[, lapply(.SD, scale, center=TRUE, scale=TRUE)]
d10k <- d10ko[, lapply(.SD, scale, center=TRUE, scale=TRUE)]


outdir <- "samp1k/"
out1k <- list()
i <- 1
print("Running mlocation")
out1k[[i]] <- mlocationDat <- mlocation(d1ko, ccol = ccol)
saveRDS(mlocationDat, file = paste0(outdir, "mlocation1k.rds"))

print("Running d1heat")
out1k[[i+1]] <- d1heatDat <- d1heat(d1k, ccol = ccol)
saveRDS(d1heatDat, file = paste0(outdir, "d1heat1k.rds"))

print("Running cumvar")
out1k[[i+2]] <- cumvarDat <- cumvar(d1k)
saveRDS(cumvarDat, file = paste0(outdir, "cumvar1k.rds"))

print("Running outliers")
out1k[[i+3]] <- outliersDat <- outliers(d1k)
saveRDS(outliersDat, file = paste0(outdir, "outliers1k.rds"))

print("Running pairHex")
out1k[[i+4]] <- pairHexDat <- invisible(pairhex(d1k, maxd = 6))
saveRDS(pairHexDat, file = paste0(outdir, "pairhexDat1k.rds"))

print("Running correlation")
out1k[[i+5]] <- corDat <- medacor(d1k, ccol = ccol)
saveRDS(corDat, file = paste0(outdir, "medacor1k.rds"))

print("Running hmc")
out1k[[i+6]] <- hmcDat <- hmc(scale(d1k, center = TRUE, scale = FALSE), 
                           maxDepth = 6, modelNames = "VVV", ccol = ccol)
saveRDS(hmcDat, file = paste0(outdir, "hmc1k.rds"))

system("say done")

outdir <- "samp10k/"
out10k <- list()
i <- 1
print("Running mlocation")
out10k[[i+0]] <- mlocationDat <- mlocation(d10ko, ccol = ccol)
saveRDS(mlocationDat, file = paste0(outdir, "mlocation10k.rds"))

print("Running d1heat")
out10k[[i+1]] <- d1heatDat <- d1heat(d10k, ccol = ccol)
saveRDS(d1heatDat, file = paste0(outdir, "d1heat10k.rds"))

print("Running cumvar")
out10k[[i+2]] <- cumvarDat <- cumvar(d10k)
saveRDS(cumvarDat, file = paste0(outdir, "cumvar10k.rds"))

print("Running outliers")
out10k[[i+3]] <- outliersDat <- outliers(d10k)
saveRDS(outliersDat, file = paste0(outdir, "outliers10k.rds"))

print("Running pairHex")
out10k[[i+4]] <- pairHexDat <- invisible(pairhex(d10k, maxd = 6))
saveRDS(pairHexDat, file = paste0(outdir, "pairhexDat10k.rds"))

print("Running correlation")
out10k[[i+5]] <- corDat <- medacor(d10k, ccol = ccol)
saveRDS(corDat, file = paste0(outdir, "medacor10k.rds"))

print("Running hmc")
out10k[[i+6]] <- hmcDat <- hmc(scale(d10k, center = TRUE, scale = FALSE), 
                           maxDepth = 6, modelNames = "VVV", ccol = ccol)
saveRDS(hmcDat, file = paste0(outdir, "hmc10k.rds"))

system("say done")
####

indir <- "samp1k/"
outdir <- "samp1k/"

files <- grep(".rds$", dir(indir, full.names = TRUE), value = TRUE)
#heatFile <- grep("heatmap.RData", dir(indir, full.names = TRUE), value = TRUE)

L <- list()

for(f in 1:length(files)){
  L[[f]] <- readRDS(files[f])
}

for(i in 1:length(L)){
  png(paste0(outdir, i, ".png"), height = 720, width = 720)
  print(plot(L[[i]]))
  dev.off()
}

if(any(grepl("*hmc*", files))){
  ind <- which(grepl("*hmc*", files))
  hmcDat <- L[grepl("*hmc*", files)][[1]]

  png(paste0(outdir, "dend.png"), height = 720, width = 720)
  plotDend(hmcDat)
  dev.off()

  h <-  8 +ifelse(dim(hmcDat$dat$sigma$dat)[3] %% 2==0,
                  dim(hmcDat$dat$sigma$dat)[3]- 2,dim(hmcDat$dat$sigma$dat)[3])

  png(paste0(outdir, "corr.png"), height = h, width = 8, units = "in", res = 72)
  plot(hmcDat$dat$sigma, ccol = L[[ind]]$ccol)
  dev.off()

  png(paste0(outdir, "clusterMeans.png"), height = 8, width = 18, units = "in", res = 300)
  show(clusterMeans(hmcDat, ccol = L[[ind]]$ccol))
  dev.off()

  png(paste0(outdir, "stackM.png"), height = 8, width = 8, units = "in", res = 300)
  show(stackM(hmcDat, centered = TRUE, ccol = L[[ind]]$ccol))
  dev.off()
}

###

indir <- "samp10k/"
outdir <- "samp10k/"

files <- grep(".rds$", dir(indir, full.names = TRUE), value = TRUE)
#heatFile <- grep("heatmap.RData", dir(indir, full.names = TRUE), value = TRUE)

L <- list()

for(f in 1:length(files)){
  L[[f]] <- readRDS(files[f])
}

for(i in 1:length(L)){
  png(paste0(outdir, i, ".png"), height = 720, width = 720)
  print(plot(L[[i]]))
  dev.off()
}

if(any(grepl("*hmc*", files))){
  ind <- which(grepl("*hmc*", files))
  hmcDat <- L[grepl("*hmc*", files)][[1]]

  png(paste0(outdir, "dend.png"), height = 720, width = 720)
  plotDend(hmcDat)
  dev.off()

  h <-  8 +ifelse(dim(hmcDat$dat$sigma$dat)[3] %% 2==0,
                  dim(hmcDat$dat$sigma$dat)[3]- 2,dim(hmcDat$dat$sigma$dat)[3])

  png(paste0(outdir, "corr.png"), height = h, width = 8, units = "in", res = 72)
  plot(hmcDat$dat$sigma, ccol = L[[ind]]$ccol)
  dev.off()

  png(paste0(outdir, "clusterMeans.png"), height = 8, width = 18, units = "in", res = 300)
  show(clusterMeans(hmcDat, ccol = L[[ind]]$ccol))
  dev.off()

  png(paste0(outdir, "stackM.png"), height = 8, width = 8, units = "in", res = 300)
  show(stackM(hmcDat, centered = TRUE, ccol = L[[ind]]$ccol))
  dev.off()
}


