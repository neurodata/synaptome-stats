# Synapse Clustering: Y2 Progress Report
Jesse Leigh Patsolic  
`r Sys.Date()`  




[Homepage](http://docs.neurodata.io/synaptome-stats/)  
The formatted source code for this file is [here](https://github.com/neurodata/synaptome-stats/blob/gh-pages/Code/Y2progress.Rmd).  
And a [raw version here](https://raw.githubusercontent.com/neurodata/synaptome-stats/gh-pages/Code/Y2progress.Rmd).    
Previous work by Youngser Park can be found [here](http://www.cis.jhu.edu/~parky/Synapse/synapse.html).  


# Introduction

> On Fri, Dec 11, 2015 at 11:53 AM, joshua vogelstein <jovo@jhu.edu> wrote:  
> we will get n=10^6 points, each in d=25 dimensions.  
> i want to hierarchically cluster them, in a ways:  

1. recursive k-means on the data, maybe 5 levels
2. compute approximate k-neighbors, svd in d=dimensions, and then #1
3. maybe some other ways.

# Data
> This corresponds to 24 channels x 6 features per synapse, ordered like
> c0f0,c0f1,c0f2,c0f3,c0f4,c0f5,c1f0,c1f1... etc
>
><FONT COLOR=#ff0000>f0 = integrated brightness </FONT>  
>f1 = local brightness  
>f2 = distance to Center of Mass  
>f3 = moment of inertia around synapsin maxima  
>f4,f5 are features that I forget what they are.. would need to ask brad.   
>i would throw them out, I did so in my kohonen code (which you have, its in matlab).

and

> On Feb 8, 2016, at 2:00 PM, Kristina Micheva <kmicheva@stanford.edu> wrote:

* <FONT COLOR=#197300>_Excitatory presynaptic: 'Synap', 'Synap', 'VGlut1', 'VGlut1', 'VGlut2'_</FONT>,
* <FONT COLOR=#5ed155>_Excitatory postsynaptic: 'psd', 'glur2', 'nmdar1', 'nr2b', 'NOS', 'Synapo'_</FONT> (but further away than PSD, gluR2, nmdar1 and nr2b)
* <FONT COLOR=#660000>_Inhibitory presynaptic: 'gad', 'VGAT', 'PV'_</FONT>,
* <FONT COLOR=#ff3333>_Inhibitory postsynaptic: 'Gephyr', 'GABAR1', 'GABABR', 'NOS'_</FONT>,
* <FONT COLOR=#ff9933>_At a very small number of inhibitory: 'Vglut3' (presynaptic), 'CR1'(presynaptic)_</FONT>,
* <FONT COLOR="mediumblue">_Other synapses:'5HT1A', 'TH', 'VACht'_</FONT>,
* <FONT COLOR="gold">_Not at synapses: 'tubuli', 'DAPI'_</FONT>.

and 

> On March 10, 2016, 00:29:04 (UTC), Kristina Micheva <kmicheva@stanford.edu> wrote:

There are 2 different Synap channels (2 different antibodies were
used), so that part is fine.
And 2 different VGluT1 channels (same antibody but done at different
times)
The NOS channel is the same, so count it as one even though it appears
twice. It is listed two times because it can be found at both excitatory
and inhibitory synapses. This is where your count of 25 comes, even
though there are 24 channels.
I would also add the 2 Synap channels to the Inhibitory presynaptic
category - there is supposed to be synapsin there, but at lower levels
compared to excitatory presynaptic category.

- Note:  The order of the channels are given by line `227` in the `kohenen.m` file which can be found in the dropbox folder. 
- `Synap` and `Synap` have been augmented to `Synap_1` and `Synap_2` for clarity. 
- `VGlut1` and `VGlut1` have been augmented to `VGlut1_t1` and `VGlut1_t2` to distinguish between the different times of collection (which are unknown).

## Potential Filtering 

Following is a discussion on which subset of markers could be used as a
subset to explore.

> On Thu, Apr 14, 2016 at 3:05 AM, Kristina Micheva kmicheva@stanford.edu wrote:
>
>I suggest: 
>Synap, VGluT1, VGluT2, psd, gad, vgat, gephyr,
>
>Or a bit bigger:
>Synap, VGluT1, VGluT2, psd, gad, vgat, gephyr, VGlut3, CB1
>
>> On Apr 12, 2016, at 9:54 AM, Jesse L. Patsolic studiojlp@gmail.com wrote:
>> 
>> Kristina,
>> 
>> Out of the markers available, which do you think are the best to use as a subset?
>> 

This subset has not yet been explored.


```r
featFull <- fread("../Data/synapsinR_7thA.tif.Pivots.txt.2011Features.txt",showProgress=FALSE)

### Setting a seed and creating an index vector
### to select half of the data
set.seed(2^10)
half1 <- sample(dim(featFull)[1],dim(featFull)[1]/2)
feat <- featFull[half1,]
dim(feat)
```

```
# [1] 559649    144
```

```r
channel <- c('Synap_1','Synap_2','VGlut1_t1','VGlut1_t2','VGlut2','Vglut3',
              'psd','glur2','nmdar1','nr2b','gad','VGAT',
              'PV','Gephyr','GABAR1','GABABR','CR1','5HT1A',
              'NOS','TH','VACht','Synapo','tubuli','DAPI')
channel.type <- c('ex.pre','ex.pre','ex.pre','ex.pre','ex.pre','in.pre.small',
                  'ex.post','ex.post','ex.post','ex.post','in.pre','in.pre',
                  'in.pre','in.post','in.post','in.post','in.pre.small','other',
                  'ex.post','other','other','ex.post','none','none')
nchannel <- length(channel)
nfeat <- ncol(feat) / nchannel
ffchannel <- (factor(channel.type,
    levels= c("ex.pre","ex.post","in.pre","in.post","in.pre.small","other","none")
    ))
fchannel <- as.numeric(factor(channel.type,
    levels= c("ex.pre","ex.post","in.pre","in.post","in.pre.small","other","none")
    ))
ford <- order(fchannel)
Syncol <- c("#197300","#5ed155","#660000","#cc0000","#ff9933","mediumblue","gold")
ccol <- Syncol[fchannel]

exType <- factor(c(rep("ex",11),rep("in",6),rep("other",7)),ordered=TRUE)
exCol<-exType;levels(exCol) <- c("#197300","#990000","mediumblue");
exCol <- as.character(exCol)

fname <- as.vector(sapply(channel,function(x) paste0(x,paste0("F",0:5))))
names(feat) <- fname
fcol <- rep(ccol, each=6)
mycol <- colorpanel(100, "purple", "black", "green")
mycol2 <- matlab.like(nchannel)
```

## Transformations: Considering only `f0` Integrated Brightness

We will consider only the `f0` (integrated brightness) features, and will transform the feature vector data by
scaling and filtering.



```r
f0 <- seq(1,ncol(feat),by=nfeat)
featF0 <- subset(feat, select=f0)
f01e3 <- 1e3*data.table(apply(X=featF0, 2, function(x){((x-min(x))/(max(x)-min(x)))}))

fs <- f01e3

### Taking log_10 on data + 1.
log1f <- log10(featF0 + 1)
slog1f <- data.table(scale(log1f, center=TRUE,scale=TRUE))
```

We now have the following data sets:

- `featF0`: The feature vector looking only at the integrated brightness features.
- `fs`:  The feature vector scaled between $[0,1000]$.
- `logf1`: The feature vector, plus one, then $log_{10}$ is applied. 
- `slog1f`: The feature vector, plus one, $log_{10}$, then scaled by
  subtracting the mean and dividing by the sample standard deviation.

### Kernel Density Estimates of the marginals


```r
df <- melt(as.matrix(log1f))
names(df) <- c("ind","channel","value")
df$type <- factor(rep(ffchannel,each=dim(fs)[1]),levels=levels(ffchannel))

lvo <- c(1:5,7:10,19,22,11:16,6,17,18,20,21,23,24)
levels(df$channel)<-levels(df$channel)[lvo]

ts <- 22

gg1 <- ggplot(df, aes(x=value)) + 
    scale_color_manual(values=ccol[lvo]) +
    scale_fill_manual(values=ccol[lvo]) +
    geom_histogram(aes(y=..density..,group=channel,colour=channel),bins=100) +
    geom_density(aes(group=channel, color=channel),size=1.5) +
    facet_wrap( ~ channel, scale='free', ncol=6) +
    theme(plot.title=element_text(size=ts),
          axis.title.x=element_text(size=ts),
          axis.title.y=element_text(size=ts),
          legend.title=element_text(size=ts),
          legend.text=element_text(size=ts-2),
          axis.text=element_text(size=ts-2),
          strip.text=element_text(size=ts), 
          legend.position='none')+
    ggtitle("Kernel Density Estimates of `log1f` data.")

print(gg1)
```

<figure><img src="../Figures/Y2progress_figure/cc_kde1-1.png"><figcaption><b>Figure 1: Kernel density estimates for each channel, on `log` data.</b><br><br></figcaption></figure>



## Correlations


```r
tmp <- as.numeric(table(fchannel))
cmatslog1f <- cor(slog1f)
corrplot(cmatslog1f[ford,ford],method="color",tl.col=ccol[ford], tl.cex=0.8)
corrRect(tmp,col=Syncol,lwd=4)
```

<figure><img src="../Figures/Y2progress_figure/cc_corLog-1.png"><figcaption><b>Figure 2: Correlation on log_10  data, reordered by synapse type.</b><br><br></figcaption></figure>
### PCA on the Correlation Matrix


```r
pcaLog <- prcomp(cmatslog1f,scale=TRUE, center=TRUE)
elLog <- getElbows(pcaLog$sdev, plot=FALSE) 
```

We run K-means for $K=3$ on the PCA embedded in $\mathbb{R}^3$ of the correlation matrix.


```r
K1 <- c(3)  ## The set of K's.

## Run kmeans on the pca of the correlation matrix of the slog1f data
kvecCor <- foreach(i = K1) %dopar% {
    set.seed(2^4 - 1)
    kmeans(pcaLog$x[,1:3],centers=i)
}
```
### Plots of embeddings


```r
par(mfrow=c(1,2))

plot(pcaLog$x[,1:2],
     col=ccol,
     pch=as.numeric(exType)+15,
     cex=1.5,
     xlim=c(min(pcaLog$x[,1])-0.2,max(pcaLog$x[,1])+0.7),
     main="Embedding of PCA log_10 correlation data")
text(pcaLog$x[,1:2],label=abbreviate(channel),offset=1, pos=4)


plot(pcaLog$x[,1:2],
     col=as.numeric(kvecCor[[1]]$cluster)+1,
     pch=as.numeric(kvecCor[[1]]$cluster)-1,
     cex=1.5,
     xlim=c(min(pcaLog$x[,1])-0.2,max(pcaLog$x[,1])+0.7),
     main="Embedding of PCA log_10 correlation data, \ncolor based on 3-means clustering.")
text(pcaLog$x[,1:2],col='black',label=abbreviate(channel),offset=1, pos=4)
```

<figure><img src="../Figures/Y2progress_figure/cc_2dEmb-1.png"><figcaption><b>Figure 3: 2d Embeddings.</b><br><br></figcaption></figure>


```r
pca <- pcaLog$x[,1:3]
rgl::plot3d(pca[,1],pca[,2],pca[,3],type='s',col=ccol, size=1, main="Log")
rgl::rgl.texts(pca[,1],pca[,2],pca[,3],abbreviate(channel), col=ccol, adj=c(0,1.5))
subid <- currentSubscene3d()
rglwidget(elementId="plot3dLog")
```

<!--html_preserve--><div id="plot3dLog" style="width:908px;height:908px;" class="rglWebGL html-widget"></div>
<script type="application/json" data-for="plot3dLog">{"x":{"material":{"color":["#197300","#197300","#197300","#197300","#197300","#FF9933","#5ED155","#5ED155","#5ED155","#5ED155","#660000","#660000","#660000","#CC0000","#CC0000","#CC0000","#FF9933","#0000CD","#5ED155","#0000CD","#0000CD","#5ED155","#FFD700","#FFD700"],"alpha":1,"lit":true,"ambient":"#000000","specular":"#FFFFFF","emission":"#000000","shininess":50,"smooth":true,"front":"filled","back":"filled","size":3,"lwd":1,"fog":true,"point_antialias":false,"line_antialias":false,"texture":null,"textype":"rgb","texmipmap":false,"texminfilter":"linear","texmagfilter":"linear","texenvmap":false,"depth_mask":true,"depth_test":"less"},"rootSubscene":1,"objects":{"7":{"id":7,"type":"spheres","material":{"fog":false},"vertices":[[-3.38446593284607,-2.29389095306396,0.618443965911865],[-1.69793093204498,-2.60134029388428,1.39720737934113],[-4.48616218566895,-0.175301522016525,0.443882137537003],[-4.33124923706055,-1.40111243724823,-0.733580768108368],[-1.57346439361572,-0.218749970197678,-1.09146070480347],[1.93022656440735,2.07484555244446,0.889972567558289],[-3.91450786590576,-1.16529226303101,-0.820843100547791],[-2.96932363510132,-0.10941506922245,-0.754895508289337],[-1.17625892162323,1.42155623435974,0.665551960468292],[-1.74300587177277,1.85960602760315,1.69779121875763],[4.8787260055542,-1.8766348361969,0.807610869407654],[4.36625671386719,-2.14910554885864,1.11511433124542],[3.4916307926178,-1.24088978767395,0.675448060035706],[3.08042025566101,-1.27950263023376,-0.238880902528763],[1.89296567440033,-1.21416342258453,-0.310152143239975],[-2.99575209617615,0.847337067127228,1.39195263385773],[1.78843986988068,1.63774275779724,-0.901631057262421],[1.15153658390045,-0.276376575231552,1.79671812057495],[0.0864422842860222,2.58966636657715,1.82439792156219],[0.822415292263031,0.74241691827774,-1.41220259666443],[1.45771813392639,1.77055132389069,-2.44107270240784],[-0.121243894100189,2.49075078964233,-1.5885294675827],[1.30376839637756,-1.61212015151978,-3.41424798965454],[2.14281845092773,2.17942237854004,0.383405685424805]],"colors":[[0.0980392172932625,0.450980395078659,0,1],[0.0980392172932625,0.450980395078659,0,1],[0.0980392172932625,0.450980395078659,0,1],[0.0980392172932625,0.450980395078659,0,1],[0.0980392172932625,0.450980395078659,0,1],[1,0.600000023841858,0.200000002980232,1],[0.368627458810806,0.819607853889465,0.333333343267441,1],[0.368627458810806,0.819607853889465,0.333333343267441,1],[0.368627458810806,0.819607853889465,0.333333343267441,1],[0.368627458810806,0.819607853889465,0.333333343267441,1],[0.400000005960464,0,0,1],[0.400000005960464,0,0,1],[0.400000005960464,0,0,1],[0.800000011920929,0,0,1],[0.800000011920929,0,0,1],[0.800000011920929,0,0,1],[1,0.600000023841858,0.200000002980232,1],[0,0,0.803921580314636,1],[0.368627458810806,0.819607853889465,0.333333343267441,1],[0,0,0.803921580314636,1],[0,0,0.803921580314636,1],[0.368627458810806,0.819607853889465,0.333333343267441,1],[1,0.843137264251709,0,1],[1,0.843137264251709,0,1]],"radii":[[0.114702142775059]],"centers":[[-3.38446593284607,-2.29389095306396,0.618443965911865],[-1.69793093204498,-2.60134029388428,1.39720737934113],[-4.48616218566895,-0.175301522016525,0.443882137537003],[-4.33124923706055,-1.40111243724823,-0.733580768108368],[-1.57346439361572,-0.218749970197678,-1.09146070480347],[1.93022656440735,2.07484555244446,0.889972567558289],[-3.91450786590576,-1.16529226303101,-0.820843100547791],[-2.96932363510132,-0.10941506922245,-0.754895508289337],[-1.17625892162323,1.42155623435974,0.665551960468292],[-1.74300587177277,1.85960602760315,1.69779121875763],[4.8787260055542,-1.8766348361969,0.807610869407654],[4.36625671386719,-2.14910554885864,1.11511433124542],[3.4916307926178,-1.24088978767395,0.675448060035706],[3.08042025566101,-1.27950263023376,-0.238880902528763],[1.89296567440033,-1.21416342258453,-0.310152143239975],[-2.99575209617615,0.847337067127228,1.39195263385773],[1.78843986988068,1.63774275779724,-0.901631057262421],[1.15153658390045,-0.276376575231552,1.79671812057495],[0.0864422842860222,2.58966636657715,1.82439792156219],[0.822415292263031,0.74241691827774,-1.41220259666443],[1.45771813392639,1.77055132389069,-2.44107270240784],[-0.121243894100189,2.49075078964233,-1.5885294675827],[1.30376839637756,-1.61212015151978,-3.41424798965454],[2.14281845092773,2.17942237854004,0.383405685424805]],"ignoreExtent":false,"flags":3},"9":{"id":9,"type":"text","material":{"lit":false,"fog":false},"vertices":[[0.196281909942627,3.58538818359375,2.82925748825073]],"colors":[[0,0,0,1]],"texts":[["Log"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[0.196281909942627,3.58538818359375,2.82925748825073]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"flags":40},"10":{"id":10,"type":"text","material":{"lit":false,"fog":false},"vertices":[[0.196281909942627,-3.59706211090088,-4.41910743713379]],"colors":[[0,0,0,1]],"texts":[["pca[, 1]"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[0.196281909942627,-3.59706211090088,-4.41910743713379]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"flags":40},"11":{"id":11,"type":"text","material":{"lit":false,"fog":false},"vertices":[[-6.28250408172607,-0.00583696365356445,-4.41910743713379]],"colors":[[0,0,0,1]],"texts":[["pca[, 2]"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-6.28250408172607,-0.00583696365356445,-4.41910743713379]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"flags":40},"12":{"id":12,"type":"text","material":{"lit":false,"fog":false},"vertices":[[-6.28250408172607,-3.59706211090088,-0.794925034046173]],"colors":[[0,0,0,1]],"texts":[["pca[, 3]"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-6.28250408172607,-3.59706211090088,-0.794925034046173]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"flags":40},"13":{"id":13,"type":"text","material":{"lit":false},"vertices":[[-3.38446593284607,-2.29389095306396,0.618443965911865],[-1.69793093204498,-2.60134029388428,1.39720737934113],[-4.48616218566895,-0.175301522016525,0.443882137537003],[-4.33124923706055,-1.40111243724823,-0.733580768108368],[-1.57346439361572,-0.218749970197678,-1.09146070480347],[1.93022656440735,2.07484555244446,0.889972567558289],[-3.91450786590576,-1.16529226303101,-0.820843100547791],[-2.96932363510132,-0.10941506922245,-0.754895508289337],[-1.17625892162323,1.42155623435974,0.665551960468292],[-1.74300587177277,1.85960602760315,1.69779121875763],[4.8787260055542,-1.8766348361969,0.807610869407654],[4.36625671386719,-2.14910554885864,1.11511433124542],[3.4916307926178,-1.24088978767395,0.675448060035706],[3.08042025566101,-1.27950263023376,-0.238880902528763],[1.89296567440033,-1.21416342258453,-0.310152143239975],[-2.99575209617615,0.847337067127228,1.39195263385773],[1.78843986988068,1.63774275779724,-0.901631057262421],[1.15153658390045,-0.276376575231552,1.79671812057495],[0.0864422842860222,2.58966636657715,1.82439792156219],[0.822415292263031,0.74241691827774,-1.41220259666443],[1.45771813392639,1.77055132389069,-2.44107270240784],[-0.121243894100189,2.49075078964233,-1.5885294675827],[1.30376839637756,-1.61212015151978,-3.41424798965454],[2.14281845092773,2.17942237854004,0.383405685424805]],"colors":[[0.0980392172932625,0.450980395078659,0,1],[0.0980392172932625,0.450980395078659,0,1],[0.0980392172932625,0.450980395078659,0,1],[0.0980392172932625,0.450980395078659,0,1],[0.0980392172932625,0.450980395078659,0,1],[1,0.600000023841858,0.200000002980232,1],[0.368627458810806,0.819607853889465,0.333333343267441,1],[0.368627458810806,0.819607853889465,0.333333343267441,1],[0.368627458810806,0.819607853889465,0.333333343267441,1],[0.368627458810806,0.819607853889465,0.333333343267441,1],[0.400000005960464,0,0,1],[0.400000005960464,0,0,1],[0.400000005960464,0,0,1],[0.800000011920929,0,0,1],[0.800000011920929,0,0,1],[0.800000011920929,0,0,1],[1,0.600000023841858,0.200000002980232,1],[0,0,0.803921580314636,1],[0.368627458810806,0.819607853889465,0.333333343267441,1],[0,0,0.803921580314636,1],[0,0,0.803921580314636,1],[0.368627458810806,0.819607853889465,0.333333343267441,1],[1,0.843137264251709,0,1],[1,0.843137264251709,0,1]],"texts":[["Sy_1"],["Sy_2"],["VG1_1"],["VG1_2"],["VGl2"],["Vgl3"],["psd"],["glr2"],["nmd1"],["nr2b"],["gad"],["VGAT"],["PV"],["Gphy"],["GABAR"],["GABAB"],["CR1"],["5HT1"],["NOS"],["TH"],["VACh"],["Synp"],["tubl"],["DAPI"]],"cex":[[1]],"adj":[[0,1.5]],"centers":[[-3.38446593284607,-2.29389095306396,0.618443965911865],[-1.69793093204498,-2.60134029388428,1.39720737934113],[-4.48616218566895,-0.175301522016525,0.443882137537003],[-4.33124923706055,-1.40111243724823,-0.733580768108368],[-1.57346439361572,-0.218749970197678,-1.09146070480347],[1.93022656440735,2.07484555244446,0.889972567558289],[-3.91450786590576,-1.16529226303101,-0.820843100547791],[-2.96932363510132,-0.10941506922245,-0.754895508289337],[-1.17625892162323,1.42155623435974,0.665551960468292],[-1.74300587177277,1.85960602760315,1.69779121875763],[4.8787260055542,-1.8766348361969,0.807610869407654],[4.36625671386719,-2.14910554885864,1.11511433124542],[3.4916307926178,-1.24088978767395,0.675448060035706],[3.08042025566101,-1.27950263023376,-0.238880902528763],[1.89296567440033,-1.21416342258453,-0.310152143239975],[-2.99575209617615,0.847337067127228,1.39195263385773],[1.78843986988068,1.63774275779724,-0.901631057262421],[1.15153658390045,-0.276376575231552,1.79671812057495],[0.0864422842860222,2.58966636657715,1.82439792156219],[0.822415292263031,0.74241691827774,-1.41220259666443],[1.45771813392639,1.77055132389069,-2.44107270240784],[-0.121243894100189,2.49075078964233,-1.5885294675827],[1.30376839637756,-1.61212015151978,-3.41424798965454],[2.14281845092773,2.17942237854004,0.383405685424805]],"family":[["sans"]],"font":[[1]],"ignoreExtent":false,"flags":40},"5":{"id":5,"type":"light","vertices":[[0,0,1]],"colors":[[1,1,1,1],[1,1,1,1],[1,1,1,1]],"viewpoint":true,"finite":false},"4":{"id":4,"type":"background","colors":[[0.298039227724075,0.298039227724075,0.298039227724075,1]],"centers":[[0,0,0]],"sphere":false,"fogtype":"none"},"6":{"id":6,"type":"background","colors":[[1,1,1,1]],"centers":[[0,0,0]],"sphere":false,"fogtype":"none"},"8":{"id":8,"type":"bboxdeco","material":{"color":"#000000","front":"lines","back":"lines","fog":false},"vertices":[[-4,null,null],[-2,null,null],[0,null,null],[2,null,null],[4,null,null],[null,-2,null],[null,-1,null],[null,0,null],[null,1,null],[null,2,null],[null,null,-3],[null,null,-2],[null,null,-1],[null,null,0],[null,null,1]],"colors":[[0,0,0,1]],"draw_front":true,"newIds":[21,22,23,24,25,26,27]},"1":{"id":1,"type":"subscene","par3d":{"antialias":8,"FOV":30,"ignoreExtent":false,"listeners":1,"mouseMode":{"left":"trackball","right":"zoom","middle":"fov","wheel":"pull"},"observer":[0,0,30.1411647796631],"modelMatrix":[[0.734886407852173,0,0,-0.144244909286499],[0,0.453443199396133,1.23449563980103,0.983978152275085],[0,-1.24582493305206,0.449319690465927,-29.7912616729736],[0,0,0,1]],"projMatrix":[[3.73205065727234,0,0,0],[0,3.73205065727234,0,0],[0,0,-3.86370325088501,-108.655410766602],[0,0,-1,0]],"skipRedraw":false,"userMatrix":[[1,0,0,0],[0,0.342020143325668,0.939692620785909,0],[0,-0.939692620785909,0.342020143325668,0],[0,0,0,1]],"scale":[0.734886407852173,1.32577919960022,1.31372284889221],"viewport":{"x":0,"y":0,"width":1,"height":1},"zoom":1,"bbox":[-4.6422438621521,5.03480768203735,-2.68785715103149,2.67618322372437,-3.50155878067017,1.91170871257782],"windowRect":[2560,45,2816,301],"family":"sans","font":1,"cex":1,"useFreeType":true,"fontname":"/Users/JLP/R_libs/rgl/fonts/FreeSans.ttf","maxClipPlanes":6},"embeddings":{"viewport":"replace","projection":"replace","model":"replace"},"objects":[6,8,7,9,10,11,12,13,5,21,22,23,24,25,26,27],"subscenes":[],"flags":1195},"21":{"id":21,"type":"lines","material":{"lit":false,"front":"lines","back":"lines"},"vertices":[[-4,-2.76831769943237,-3.58275771141052],[4,-2.76831769943237,-3.58275771141052],[-4,-2.76831769943237,-3.58275771141052],[-4,-2.9064416885376,-3.72214937210083],[-2,-2.76831769943237,-3.58275771141052],[-2,-2.9064416885376,-3.72214937210083],[0,-2.76831769943237,-3.58275771141052],[0,-2.9064416885376,-3.72214937210083],[2,-2.76831769943237,-3.58275771141052],[2,-2.9064416885376,-3.72214937210083],[4,-2.76831769943237,-3.58275771141052],[4,-2.9064416885376,-3.72214937210083]],"colors":[[0,0,0,1]],"centers":[[0,-2.76831769943237,-3.58275771141052],[-4,-2.83737969398499,-3.65245342254639],[-2,-2.83737969398499,-3.65245342254639],[0,-2.83737969398499,-3.65245342254639],[2,-2.83737969398499,-3.65245342254639],[4,-2.83737969398499,-3.65245342254639]],"ignoreExtent":true,"origId":8,"flags":128},"22":{"id":22,"type":"text","material":{"lit":false,"front":"lines","back":"lines"},"vertices":[[-4,-3.18268990516663,-4.00093269348145],[-2,-3.18268990516663,-4.00093269348145],[0,-3.18268990516663,-4.00093269348145],[2,-3.18268990516663,-4.00093269348145],[4,-3.18268990516663,-4.00093269348145]],"colors":[[0,0,0,1]],"texts":[["-4"],["-2"],["0"],["2"],["4"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-4,-3.18268990516663,-4.00093269348145],[-2,-3.18268990516663,-4.00093269348145],[0,-3.18268990516663,-4.00093269348145],[2,-3.18268990516663,-4.00093269348145],[4,-3.18268990516663,-4.00093269348145]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"origId":8,"flags":40},"23":{"id":23,"type":"lines","material":{"lit":false,"front":"lines","back":"lines"},"vertices":[[-4.78739976882935,-2,-3.58275771141052],[-4.78739976882935,2,-3.58275771141052],[-4.78739976882935,-2,-3.58275771141052],[-5.03658390045166,-2,-3.72214937210083],[-4.78739976882935,-1,-3.58275771141052],[-5.03658390045166,-1,-3.72214937210083],[-4.78739976882935,0,-3.58275771141052],[-5.03658390045166,0,-3.72214937210083],[-4.78739976882935,1,-3.58275771141052],[-5.03658390045166,1,-3.72214937210083],[-4.78739976882935,2,-3.58275771141052],[-5.03658390045166,2,-3.72214937210083]],"colors":[[0,0,0,1]],"centers":[[-4.78739976882935,0,-3.58275771141052],[-4.91199207305908,-2,-3.65245342254639],[-4.91199207305908,-1,-3.65245342254639],[-4.91199207305908,0,-3.65245342254639],[-4.91199207305908,1,-3.65245342254639],[-4.91199207305908,2,-3.65245342254639]],"ignoreExtent":true,"origId":8,"flags":128},"24":{"id":24,"type":"text","material":{"lit":false,"front":"lines","back":"lines"},"vertices":[[-5.53495168685913,-2,-4.00093269348145],[-5.53495168685913,-1,-4.00093269348145],[-5.53495168685913,0,-4.00093269348145],[-5.53495168685913,1,-4.00093269348145],[-5.53495168685913,2,-4.00093269348145]],"colors":[[0,0,0,1]],"texts":[["-2"],["-1"],["0"],["1"],["2"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-5.53495168685913,-2,-4.00093269348145],[-5.53495168685913,-1,-4.00093269348145],[-5.53495168685913,0,-4.00093269348145],[-5.53495168685913,1,-4.00093269348145],[-5.53495168685913,2,-4.00093269348145]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"origId":8,"flags":40},"25":{"id":25,"type":"lines","material":{"lit":false,"front":"lines","back":"lines"},"vertices":[[-4.78739976882935,-2.76831769943237,-3],[-4.78739976882935,-2.76831769943237,1],[-4.78739976882935,-2.76831769943237,-3],[-5.03658390045166,-2.9064416885376,-3],[-4.78739976882935,-2.76831769943237,-2],[-5.03658390045166,-2.9064416885376,-2],[-4.78739976882935,-2.76831769943237,-1],[-5.03658390045166,-2.9064416885376,-1],[-4.78739976882935,-2.76831769943237,0],[-5.03658390045166,-2.9064416885376,0],[-4.78739976882935,-2.76831769943237,1],[-5.03658390045166,-2.9064416885376,1]],"colors":[[0,0,0,1]],"centers":[[-4.78739976882935,-2.76831769943237,-1],[-4.91199207305908,-2.83737969398499,-3],[-4.91199207305908,-2.83737969398499,-2],[-4.91199207305908,-2.83737969398499,-1],[-4.91199207305908,-2.83737969398499,0],[-4.91199207305908,-2.83737969398499,1]],"ignoreExtent":true,"origId":8,"flags":128},"26":{"id":26,"type":"text","material":{"lit":false,"front":"lines","back":"lines"},"vertices":[[-5.53495168685913,-3.18268990516663,-3],[-5.53495168685913,-3.18268990516663,-2],[-5.53495168685913,-3.18268990516663,-1],[-5.53495168685913,-3.18268990516663,0],[-5.53495168685913,-3.18268990516663,1]],"colors":[[0,0,0,1]],"texts":[["-3"],["-2"],["-1"],["0"],["1"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-5.53495168685913,-3.18268990516663,-3],[-5.53495168685913,-3.18268990516663,-2],[-5.53495168685913,-3.18268990516663,-1],[-5.53495168685913,-3.18268990516663,0],[-5.53495168685913,-3.18268990516663,1]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"origId":8,"flags":40},"27":{"id":27,"type":"lines","material":{"lit":false,"front":"lines","back":"lines"},"vertices":[[-4.78739976882935,-2.76831769943237,-3.58275771141052],[-4.78739976882935,2.75664377212524,-3.58275771141052],[-4.78739976882935,-2.76831769943237,1.99290776252747],[-4.78739976882935,2.75664377212524,1.99290776252747],[-4.78739976882935,-2.76831769943237,-3.58275771141052],[-4.78739976882935,-2.76831769943237,1.99290776252747],[-4.78739976882935,2.75664377212524,-3.58275771141052],[-4.78739976882935,2.75664377212524,1.99290776252747],[-4.78739976882935,-2.76831769943237,-3.58275771141052],[5.1799635887146,-2.76831769943237,-3.58275771141052],[-4.78739976882935,-2.76831769943237,1.99290776252747],[5.1799635887146,-2.76831769943237,1.99290776252747],[-4.78739976882935,2.75664377212524,-3.58275771141052],[5.1799635887146,2.75664377212524,-3.58275771141052],[-4.78739976882935,2.75664377212524,1.99290776252747],[5.1799635887146,2.75664377212524,1.99290776252747],[5.1799635887146,-2.76831769943237,-3.58275771141052],[5.1799635887146,2.75664377212524,-3.58275771141052],[5.1799635887146,-2.76831769943237,1.99290776252747],[5.1799635887146,2.75664377212524,1.99290776252747],[5.1799635887146,-2.76831769943237,-3.58275771141052],[5.1799635887146,-2.76831769943237,1.99290776252747],[5.1799635887146,2.75664377212524,-3.58275771141052],[5.1799635887146,2.75664377212524,1.99290776252747]],"colors":[[0,0,0,1]],"centers":[[-4.78739976882935,-0.00583696365356445,-3.58275771141052],[-4.78739976882935,-0.00583696365356445,1.99290776252747],[-4.78739976882935,-2.76831769943237,-0.794924974441528],[-4.78739976882935,2.75664377212524,-0.794924974441528],[0.196281909942627,-2.76831769943237,-3.58275771141052],[0.196281909942627,-2.76831769943237,1.99290776252747],[0.196281909942627,2.75664377212524,-3.58275771141052],[0.196281909942627,2.75664377212524,1.99290776252747],[5.1799635887146,-0.00583696365356445,-3.58275771141052],[5.1799635887146,-0.00583696365356445,1.99290776252747],[5.1799635887146,-2.76831769943237,-0.794924974441528],[5.1799635887146,2.75664377212524,-0.794924974441528]],"ignoreExtent":true,"origId":8,"flags":128}},"width":257,"height":257,"sphereVerts":{"material":[],"it":[[0,6,7,19,4,8,6,22,2,7,8,25,7,6,8,26,0,7,9,27,2,10,7,24,5,9,10,32,9,7,10,33,0,11,6,18,3,12,11,37,4,6,12,39,6,11,12,40,0,9,11,34,5,13,9,31,3,11,13,44,11,9,13,45,1,14,15,47,2,8,14,49,4,15,8,21,15,14,8,52,1,16,14,46,5,10,16,55,2,14,10,29,14,16,10,57,1,15,17,58,4,12,15,51,3,17,12,36,17,15,12,62,1,17,16,53,3,13,17,61,5,16,13,42,16,17,13,65],[18,20,19,18,21,23,22,21,24,26,25,24,20,23,26,20,19,28,27,19,29,30,24,29,31,33,32,31,28,30,33,28,34,35,18,34,36,38,37,36,22,40,39,22,35,38,40,35,27,41,34,27,42,43,31,42,37,45,44,37,41,43,45,41,46,48,47,46,25,50,49,25,51,52,21,51,48,50,52,48,53,54,46,53,32,56,55,32,49,57,29,49,54,56,57,54,47,59,58,47,39,60,51,39,61,62,36,61,59,60,62,59,58,63,53,58,44,64,61,44,55,65,42,55,63,64,65,63],[19,18,20,20,22,21,23,23,25,24,26,26,26,20,23,23,27,19,28,28,24,29,30,30,32,31,33,33,33,28,30,30,18,34,35,35,37,36,38,38,39,22,40,40,40,35,38,38,34,27,41,41,31,42,43,43,44,37,45,45,45,41,43,43,47,46,48,48,49,25,50,50,21,51,52,52,52,48,50,50,46,53,54,54,55,32,56,56,29,49,57,57,57,54,56,56,58,47,59,59,51,39,60,60,36,61,62,62,62,59,60,60,53,58,63,63,61,44,64,64,42,55,65,65,65,63,64,64]],"vb":[[-1,1,0,0,0,0,-0.707106781186548,-0.707106781186548,0,-0.707106781186548,0,-0.707106781186548,0,0,0.707106781186548,0.707106781186548,0.707106781186548,0.707106781186548,-0.934997526317783,-0.934997526317783,-0.770440042047682,0,-0.354654234120539,-0.450789386304495,-0.354654234120539,0,-0.450789386304495,-0.934997526317783,-0.770440042047682,0,-0.450789386304495,-0.354654234120539,0,-0.450789386304495,-0.934997526317783,-0.770440042047682,0,-0.354654234120539,-0.450789386304495,0,-0.450789386304495,-0.770440042047682,0,-0.450789386304495,0,-0.450789386304495,0.934997526317783,0.934997526317783,0.770440042047682,0.354654234120539,0.450789386304495,0.354654234120539,0.450789386304495,0.934997526317783,0.770440042047682,0.354654234120539,0.450789386304495,0.450789386304495,0.934997526317783,0.770440042047682,0.450789386304495,0.354654234120539,0.450789386304495,0.770440042047682,0.450789386304495,0.450789386304495],[0,0,-1,1,0,0,0,-0.707106781186548,-0.707106781186548,0,-0.707106781186548,0.707106781186548,0.707106781186548,0.707106781186548,-0.707106781186548,0,0,0.707106781186548,0,-0.354654234120539,-0.450789386304495,-0.354654234120539,0,-0.450789386304495,-0.934997526317783,-0.934997526317783,-0.770440042047682,0,-0.450789386304495,-0.934997526317783,-0.770440042047682,0,-0.354654234120539,-0.450789386304495,0.354654234120539,0.450789386304495,0.934997526317783,0.934997526317783,0.770440042047682,0.354654234120539,0.450789386304495,0.450789386304495,0.354654234120539,0.450789386304495,0.934997526317783,0.770440042047682,-0.354654234120539,0,-0.450789386304495,-0.934997526317783,-0.770440042047682,0,-0.450789386304495,0,-0.450789386304495,0,-0.450789386304495,-0.770440042047682,0.354654234120539,0.450789386304495,0.450789386304495,0.934997526317783,0.770440042047682,0.450789386304495,0.770440042047682,0.450789386304495],[0,0,0,0,-1,1,-0.707106781186548,0,-0.707106781186548,0.707106781186548,0.707106781186548,0,-0.707106781186548,0.707106781186548,0,-0.707106781186548,0.707106781186548,0,-0.354654234120539,0,-0.450789386304495,-0.934997526317783,-0.934997526317783,-0.770440042047682,0,-0.354654234120539,-0.450789386304495,0.354654234120539,0.450789386304495,0.354654234120539,0.450789386304495,0.934997526317783,0.934997526317783,0.770440042047682,0,-0.450789386304495,-0.354654234120539,0,-0.450789386304495,-0.934997526317783,-0.770440042047682,0.450789386304495,0.934997526317783,0.770440042047682,0.354654234120539,0.450789386304495,0,-0.354654234120539,-0.450789386304495,0,-0.450789386304495,-0.934997526317783,-0.770440042047682,0.354654234120539,0.450789386304495,0.934997526317783,0.770440042047682,0.450789386304495,0,-0.450789386304495,-0.770440042047682,0,-0.450789386304495,0.450789386304495,0.450789386304495,0.770440042047682]],"primitivetype":"triangle"}},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->

## K-Means Level 1

Next we run K-means with $K=3$.

** <FONT COLOR=#ff3333> Note that a seed is being set for the random initialization of K-means. </FONT> **


```r
K2 <- c(2)  ## The set of K's.

## Run kmeans on the slog1f data
kvecslog1f <- foreach(i = K2) %dopar% {
    set.seed(2^4 - 1)
    kmeans(slog1f,centers=i)
}
```


### Heat maps: scaled data.

For the following we manualy choose 2 clusters.


```r
## Formatting data for heatmap
aggslog1f <- aggregate(slog1f,by=list(lab=kvecslog1f[[1]]$cluster),FUN=mean)
aggslog1f <- as.matrix(aggslog1f[,-1])
rownames(aggslog1f) <- clusterFraction(kvecslog1f[[1]])

ford <- order(fchannel)
```




```r
heatmap.2(as.matrix(aggslog1f[,ford]),dendrogram='row',Colv=NA,trace="none", col=mycol,colCol=ccol[ford],cexRow=0.8, keysize=1.25,symkey=FALSE,symbreaks=FALSE,scale="none", srtCol=90,main="Heatmap of `slog1f` data.") 
```

```
#  [1] "#197300"    "#197300"    "#197300"    "#197300"    "#197300"   
#  [6] "#5ed155"    "#5ed155"    "#5ed155"    "#5ed155"    "#5ed155"   
# [11] "#5ed155"    "#660000"    "#660000"    "#660000"    "#cc0000"   
# [16] "#cc0000"    "#cc0000"    "#ff9933"    "#ff9933"    "mediumblue"
# [21] "mediumblue" "mediumblue" "gold"       "gold"
```

<figure><img src="../Figures/Y2progress_figure/cc_km2-heatmapSorted-1.png"><figcaption><b>Figure 4: Heatmap of the cluster means vs channels. Rows and columns are rearranged according to synapse type.</b><br><br></figcaption></figure>

Percentage of data within cluster is presented on the right side of the heatmap.

# Exploring pair-wise relationships


```r
### Sampling to reduce size
set.seed(2^13 - 2)
s1 <- sample(dim(slog1f)[1],2.5e5)
dlog1f <- data.table(log1f[s1,])
```

## `GABABR`


```r
## re-formatting data for use in lattice 
dlog1f2 <- data.table(stack(dlog1f, select=-GABABRF0))[,.(values)]
dlog1f2$GABABR <- dlog1f$GABABRF0

### Adding relationship factor variables
nd <- paste0("GABABR","~",abbreviate(channel[-which(channel=="GABABR")]))

dlog1f2$ind <- factor(rep(nd,each=dim(dlog1f)[1]), ordered=TRUE,levels=nd)
names(dlog1f2) <- c("x","y","g")

rg1 <- xyplot(y ~ x | g, data=dlog1f2,
       as.table=TRUE,
       colramp=BTC,
       pch='.',
       scales = list(y = list(relation = "free"),x = list(relation = "free")),
       panel=function(x,y,...){
           panel.hexbinplot(x,y,..., type='g')
           panel.loess(x,y,col='red', lwd=2,...)
        }
       )
```

<figure><img src="../Figures/Y2progress_figure/rg1-1.png"><figcaption><b>Figure 5: Lattice plot of pairwise regressions involving `GABABR`</b><br><br></figcaption></figure>

# ARI


Here we consider the channel types to be the "ground truth" and computer
the Adjusted Rand Index of between that and the output from k-means.


## Approximate permutation test.

```r
levels(ffchannel) <- c(rep("ex", 2), rep("in", 2), rep("other", 3))
levels(ffchannel)
```

```
# [1] "ex"    "in"    "other"
```

```r
truth <- as.numeric(ffchannel)
```





Making a data.table of the permutation data for ggplot.

```r
DT <- data.table(avd(pcaLog, elLog[2], truth),key='Embedding_Dimension') 

DT <- DT[,pval := sum(permARI>=ari)/length(permARI),by=Embedding_Dimension]

ua <- DT[,unique(ari),by=Embedding_Dimension]
arid <- data.frame(Embedding_Dimension=as.numeric(ua$Emb),
                   ARI=ua$V1,
                   pval=DT[,unique(pval),by=Embedding_Dimension]$V1)
```



```r
gg3 <- ggplot(data=DT,aes(x=permARI, y=..density..,color=Embedding_Dimension,label=pval)) + 
        #geom_histogram(binwidth=3.49*sd(DT$permARI)*length(DT$permARI)^(-1/3)) +
        geom_histogram(bins=25)+
        geom_vline(aes(xintercept=ari),colour='darkred',size=1.2)+
        #geom_text(aes(x=(ari-ari/2),y=7))+
        theme(axis.text=element_text(size=18),
                  title=element_text(size=16),  
                  strip.text.x=element_text(size=16)) + 
        facet_wrap(~Embedding_Dimension+pval,scale='free',labeller=label_both)
print(gg3)
```

<figure><img src="../Figures/Y2progress_figure/cc-gg3-1.png"><figcaption><b>Figure 6: ARI Permutation Tests</b><br><br></figcaption></figure>


```r
gg5 <- ggplot(data=arid,aes(x=Embedding_Dimension,y=ARI,label=pval,colour=pval)) + 
        geom_line(size=1.5,colour='salmon') + geom_point(size=3) +
            #geom_text(hjust='right',vjust='top',nudge_x=0.5,nudge_y=0.01,size=5)+
            theme(axis.text=element_text(size=18),
                  title=element_text(size=16)) + 
            ggtitle("ARI vs. DIM with estimated p-values")

gg6 <- 
    ggplot(data=arid, aes(x=Embedding_Dimension, y=pval, colour=pval))+
        geom_line(size=1.5, colour='salmon') + 
        geom_point(size=3) + 
       # geom_hline(yintercept=0.05, 
       #            colour='forestgreen',
       #            size=1.5) +
        scale_y_log10(breaks=c(1e-4,1.53-4,1e-3), na.value=0) +
        theme(axis.text=element_text(size=18),
                  title=element_text(size=16)) + 
        ggtitle("P-values")

grid.arrange(gg5,gg6,nrow=2)
```

<figure><img src="../Figures/Y2progress_figure/arivDim-1.png"><figcaption><b>Figure 7: ARI and P-Values</b><br><br></figcaption></figure>


<footer>
<p> [Back to top][Introduction]</p>
</footer>
