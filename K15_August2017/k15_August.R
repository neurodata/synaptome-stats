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

ccol3 <- Syncol3[c(rep(1,11), rep(2,6), rep(3,7))]

q1 <- c("SELECT * FROM kristina15")

dat0 <- data.table(dbGetQuery(con1, q1))
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
}

samp10k <- read.csv("sample_locations_1e4.csv", header=TRUE)
samp1k  <- read.csv("sample_locations_1e3.csv", header=TRUE)


d1ko  <- dat0[samp1k, !(x:z)]
d10ko <- dat0[samp10k, !(x:z)]

d1k  <- d1ko[, lapply(.SD, scale, center=TRUE, scale=TRUE)]
d10k <- d10ko[, lapply(.SD, scale, center=TRUE, scale=TRUE)]

plot(mlocation(d1ko, ccol = ccol3))
plot(mlocation(d10ko, ccol = ccol3))

plot(d1heat(d1k, ccol = ccol3))
plot(d1heat(d10k, ccol = ccol3))

plot(outliers(d1k))
plot(outliers(d10k))

plot(medacor(d1k, ccol = ccol3))
plot(medacor(d10k, ccol = ccol3))

plot(cumvar(d1k))
plot(cumvar(d10k))

pairhex(d1k, maxd = 8)
pairhex(d10k, maxd = 8)

h1 <- hmc(d1k, ccol = ccol3)
h10 <- hmc(d10k, ccol = ccol3)

plot(h1)
plot(h10)

plotDend(h1)
plotDend(h10)

stackM(h1, centered=TRUE, ccol = ccol3)
stackM(h10, centered=TRUE, ccol = ccol3)

clusterMeans(h1, ccol = ccol3)
clusterMeans(h10, ccol = ccol3)


on.exit(dbDisconnect(con1))
