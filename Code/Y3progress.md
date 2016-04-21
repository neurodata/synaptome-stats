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
Kristina15 data set.  <FONT COLOR=#ff0000>Put a link to the data here. </FONT>  


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

## How did the puncta get to us (what computer vision)

The data were gathered as array tomogray images and processed.  
<FONT COLOR=#ff0000> How were the puncta processed.</FONT>  

> The [sic] corresponds to 24 channels x 6 features per synapse, ordered like
> c0f0,c0f1,c0f2,c0f3,c0f4,c0f5,c1f0,c1f1... etc
>
>f0 = integrated brightness  
>f1 = local brightness  
>f2 = distance to Center of Mass  
>f3 = moment of inertia around synapsin maxima  
>f4,f5 are features that I forget what they are.. would need to ask brad.   
>i would throw them out, I did so in my kohonen code (which you have, its in matlab).


<footer>
<p> [Back to Top][Outline]</p>
</footer>

