# Synapse Clustering: TRA Y3 Progress Report (JHU)
Jesse Leigh Patsolic  
`r Sys.Date()`  





# Outline

1. data  
	1. which dataset
	3. definition of markers
	1. how did puncta get to us (what computer vision)
	2. definition of features (with equations)
2. feature exploration  
	1. synapsin1 vs. synapsin2 in log domain and linear
	2. same for vglut1 & 2
	3. repeat the above 4 panel figure for each of the features
	4. kde plots of chosen transformation/feature pair
3. marker exploration  
	1. correlation matrix
	2. 2D scatter plot colored by truth
	3. 2D scatter plot colored by truth, overlay optimal voronoi diagram (fit using linear discriminant analysis, not kmeans)
	4. same as above in 3D  
	5. ARI vs. dimension using optimal voronoi diagram
4. synapse exploration
	1. kmeans (k=2) heatmap
	2. lattice plots of GABABR (fix colors)
	3. correlation matrices for each of the 2 clusters
	4. kmeans (k=2) for level 2


# Data 

## The Data Set
This report deals with the exploratory and statistical analysis of the
Kristina15 data set.  <FONT COLOR=#ff0000>Put a link to the data here and maybe a citation? </FONT>  

## Definition of Markers

From email correspondance with Kristina Micheva we have the following
definitions of the given markers.  

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
> 
> There are 2 different Synap channels (2 different antibodies were
> used), so that part is fine.
> And 2 different VGluT1 channels (same antibody but done at different
> times)
> The NOS channel is the same, so count it as one even though it appears
> twice. It is listed two times because it can be found at both excitatory
> and inhibitory synapses. This is where your count of 25 comes, even
> though there are 24 channels.
> I would also add the 2 Synap channels to the Inhibitory presynaptic
> category - there is supposed to be synapsin there, but at lower levels
> compared to excitatory presynaptic category.

- Note:  The order of the channels are given by line `227` in the `kohenen.m` file which can be found in the dropbox folder. 
- `Synap` and `Synap` have been augmented to `Synap_1` and `Synap_2` for clarity. 
- `VGlut1` and `VGlut1` have been augmented to `VGlut1_t1` and `VGlut1_t2` to distinguish between the different times of collection (which are unknown).

## Processing of puncta: computer vision

<FONT COLOR=#ff0000> How were the puncta processed?</FONT>  
Array tomography images, then what?

## Definition of features

> The [sic] corresponds to 24 channels x 6 features per synapse, ordered like
> c0f0,c0f1,c0f2,c0f3,c0f4,c0f5,c1f0,c1f1... etc
>
>f0 = integrated brightness  
>f1 = local brightness  
>f2 = distance to Center of Mass  
>f3 = moment of inertia around synapsin maxima  
>f4,f5 are features that I forget what they are.. would need to ask brad.   
>i would throw them out, I did so in my kohonen code (which you have, its in matlab).



Here we read in the data and select a random half of it for exploration. 


```r
featFull <- fread("../data/synapsinR_7thA.tif.Pivots.txt.2011Features.txt",showProgress=FALSE)

### Setting a seed and creating an index vector
### to select half of the data
set.seed(2^10)
half1 <- sample(dim(featFull)[1],dim(featFull)[1]/2)
half2 <- setdiff(1:dim(featFull)[1],half1)

feat <- featFull[half1,]
dim(feat)
```

```
# [1] 559649    144
```

```r
## Setting the channel names
channel <- c('Synap_1','Synap_2','VGlut1_t1','VGlut1_t2','VGlut2','Vglut3',
              'psd','glur2','nmdar1','nr2b','gad','VGAT',
              'PV','Gephyr','GABAR1','GABABR','CR1','5HT1A',
              'NOS','TH','VACht','Synapo','tubuli','DAPI')

## Setting the channel types
channel.type <- c('ex.pre','ex.pre','ex.pre','ex.pre','ex.pre','in.pre.small',
                  'ex.post','ex.post','ex.post','ex.post','in.pre','in.pre',
                  'in.pre','in.post','in.post','in.post','in.pre.small','other',
                  'ex.post','other','other','ex.post','none','none')

nchannel <- length(channel)
nfeat <- ncol(feat) / nchannel

## Createing factor variables for channel and channel type sorted properly
ffchannel <- (factor(channel.type,
    levels= c("ex.pre","ex.post","in.pre","in.post","in.pre.small","other","none")
    ))
fchannel <- as.numeric(factor(channel.type,
    levels= c("ex.pre","ex.post","in.pre","in.post","in.pre.small","other","none")
    ))
ford <- order(fchannel)


## Setting up colors for channel types
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

## Data transformations


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

# Feature Exploration

## Synapsin1 Vs. Synapsin2 for all features


```r
synF <- feat[, grep("Synap_", names(feat)),with=FALSE]
lsynF <- synF[, lapply(.SD, function(x){scale(log10(x+1),center=TRUE,scale=TRUE)})]

gg1 <- ggplot(data=synF,aes(x=Synap_1F0,y=Synap_2F0)) +   
        geom_point(pch='.',alpha=0.2) + 
        geom_smooth()+
        ggtitle("Untransformed Features: F0")
        
gg2 <- ggplot(data=lsynF,aes(x=Synap_1F0,y=Synap_2F0)) +
        geom_point(pch='.',alpha=0.2) +
        geom_smooth()+
        ggtitle("Scaled Logged Features: F0")

gg3 <- ggplot(data=synF,aes(x=Synap_1F1,y=Synap_2F1)) +   
        geom_point(pch='.',alpha=0.2) + 
        geom_smooth()+
        ggtitle("Untransformed Features: F1")
        
gg4 <- ggplot(data=lsynF,aes(x=Synap_1F1,y=Synap_2F1)) +   
        geom_point(pch='.',alpha=0.2) +
        geom_smooth()+
        ggtitle("Scaled Logged Features: F1")

gg5 <- ggplot(data=synF,aes(x=Synap_1F2,y=Synap_2F2)) +   
        geom_point(pch='.',alpha=0.2) + 
        geom_smooth() +
        ggtitle("Untransformed Features: F2")
        
gg6 <- ggplot(data=lsynF,aes(x=Synap_1F2,y=Synap_2F2)) +   
        geom_point(pch='.',alpha=0.2) +
        geom_smooth() +
        ggtitle("Scaled Logged Features: F2")

gg7 <- ggplot(data=synF,aes(x=Synap_1F3,y=Synap_2F3)) +   
        geom_point(pch='.',alpha=0.2) + 
        geom_smooth()+
        ggtitle("Untransformed Features: F3")
        
gg8 <- ggplot(data=lsynF,aes(x=Synap_1F3,y=Synap_2F3)) +   
        geom_point(pch='.',alpha=0.2) +
        geom_smooth()+
        ggtitle("Scaled Logged Features: F3")

gg09 <- ggplot(data=synF,aes(x=Synap_1F4,y=Synap_2F4)) +   
        geom_point(pch='.',alpha=0.2) + 
        geom_smooth()+
        ggtitle("Untransformed Features: F4")
        
gg10 <- ggplot(data=lsynF,aes(x=Synap_1F4,y=Synap_2F4)) +   
        geom_point(pch='.',alpha=0.2) +
        geom_smooth()+
        ggtitle("Scaled Logged Features: F4")

gg11 <- ggplot(data=synF,aes(x=Synap_1F5,y=Synap_2F5)) +   
        geom_point(pch='.',alpha=0.2,position='jitter') + 
        geom_smooth()+
        ggtitle("Untransformed Features: F5")
        
gg12 <- ggplot(data=lsynF,aes(x=Synap_1F5,y=Synap_2F5)) +   
        geom_point(pch='.',alpha=0.2, position='jitter') +
        geom_smooth()+
        ggtitle("Scaled Logged Features: F5")
```



```r
grid.arrange(gg1,gg2,gg3,gg4,gg5,gg6,gg7,gg8,gg09,gg10, gg11, gg12, ncol=2)
```

<figure><img src="../Figures/Y3progress_figure/cc-F1-5-1.png"><figcaption><b>Figure 2: Scatter plots of Synapsin1 and Synapsin2 on linear and log scale.</b><br><br></figcaption></figure>

## VGlut1_t1 Vs. VGlut1_t2


```r
vglutF <- feat[,grep("VGlut", names(feat)),with=FALSE]
lvglutF <- vglutF[,lapply(.SD, function(x){scale(log10(x+1),center=TRUE,scale=TRUE)})]

gg13 <- ggplot(data=vglutF,aes(x=VGlut1_t1F0,y=VGlut1_t2F0)) +   
        geom_point(pch='.',alpha=0.2) + 
        geom_smooth()+
        ggtitle("Untransformed Features: F0")
        
gg14 <- ggplot(data=lvglutF,aes(x=VGlut1_t1F0,y=VGlut1_t2F0)) +   
        geom_point(pch='.',alpha=0.2) +
        geom_smooth()+
        ggtitle("Scaled Logged Features: F0")

gg15 <- ggplot(data=vglutF,aes(x=VGlut1_t1F1,y=VGlut1_t2F1)) +   
        geom_point(pch='.',alpha=0.2) + 
        geom_smooth()+
        ggtitle("Untransformed Features: F1")
        
gg16 <- ggplot(data=lvglutF,aes(x=VGlut1_t1F1,y=VGlut1_t2F1)) +   
        geom_point(pch='.',alpha=0.2) +
        geom_smooth()+
        ggtitle("Scaled Logged Features: F1")

gg17 <- ggplot(data=vglutF,aes(x=VGlut1_t1F2,y=VGlut1_t2F2)) +   
        geom_point(pch='.',alpha=0.2) + 
        geom_smooth()+
        ggtitle("Untransformed Features: F2")
        
gg18 <- ggplot(data=lvglutF,aes(x=VGlut1_t1F2,y=VGlut1_t2F2)) +   
        geom_point(pch='.',alpha=0.2) +
        geom_smooth()+
        ggtitle("Scaled Logged Features: F2")

gg19 <- ggplot(data=vglutF,aes(x=VGlut1_t1F3,y=VGlut1_t2F3)) +   
        geom_point(pch='.',alpha=0.2) + 
        geom_smooth()+
        ggtitle("Untransformed Features: F3")
        
gg20 <- ggplot(data=lvglutF,aes(x=VGlut1_t1F3,y=VGlut1_t2F3)) +   
        geom_point(pch='.',alpha=0.2) +
        geom_smooth()+
        ggtitle("Scaled Logged Features: F3")

gg21 <- ggplot(data=vglutF,aes(x=VGlut1_t1F4,y=VGlut1_t2F4)) +   
        geom_point(pch='.',alpha=0.2) + 
        geom_smooth()+
        ggtitle("Untransformed Features: F4")
        
gg22 <- ggplot(data=lvglutF,aes(x=VGlut1_t1F4,y=VGlut1_t2F4)) +   
        geom_point(pch='.',alpha=0.2) +
        geom_smooth()+
        ggtitle("Scaled Logged Features: F4")

gg23 <- ggplot(data=vglutF,aes(x=VGlut1_t1F5,y=VGlut1_t2F5)) +   
        geom_point(pch='.',alpha=0.2, position='jitter') + 
        geom_smooth()+
        ggtitle("Untransformed Features: F5")
        
gg24 <- ggplot(data=lvglutF,aes(x=VGlut1_t1F5,y=VGlut1_t2F5)) +   
        geom_point(pch='.',alpha=0.2, position='jitter') +
        geom_smooth()+
        ggtitle("Scaled Logged Features: F5")
```


```r
grid.arrange(gg13, gg14, gg15, gg16, gg17, gg18, gg19, gg20,gg21,gg22,gg23,gg24, ncol=2)
```

<figure><img src="../Figures/Y3progress_figure/cc-vglutF1-5-1.png"><figcaption><b>Figure 4: Scatter plots of VGlut1_t1 and VGlut1_t2 on linear and log scale for all features.</b><br><br></figcaption></figure>

## KDE plots of chosen transformation/feature pair


















## Distance Covariance Test





# Marker Exploration 

## Correlation Matrix

## 2D scatter plot colored by truth

# Synapse Exploration 



<footer>
<p> [Back to Top][Outline]</p>
</footer>

