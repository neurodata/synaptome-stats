# Background on Kristina15 data.

[**Forrest Collman et al.**](http://www.jneurosci.org/content/35/14/5792.full)


##### Synapse detection
> To begin identifying synapses from GABA-negative axons, 2D
> segmentation of PSD-95 IF was performed by applying a threshold, set
> to include all IF puncta above the background autofluorescence (see
> Fig. 9a).  For quantitative consistency, this threshold was defined by
> lowering the threshold until the resulting segmentation produced a
> median PSD-95 punctum size of > 0.09 µm². We merged 2D PSD-95 puncta
> from adjacent sections into a single 3D PSD-95 punctum when their IF
> weighted centroids were < 400 nm apart (see Fig. 9c).  The threshold
> was set at 400 nm, because analysis of the conjugate AT data show that
> ~90% of all such merges resulted in linking PSD-95 puncta that overlap
> with the same synapse, without merging distinct synapses (see Fig.
> 10a).



#### Feature meanings



### Correspondence with Forrest

On Tue, May 17, 2016 at 11:48 AM, Forrest Collman <forrestc@alleninstitute.org> wrote:
> Hi Jesse, 
> 
> 
> so a clarification...
> 
> 
> you are speaking about the Kristina15 dataset, and the puncta that I sent you.
> 
> The method of detecting those puncta and quantifying them is fundamentally different
> 
> than the method that was used in my paper.
> 
> 
> This is because I wanted to get you guys something fast, that included millions of synapses of a variety of types, for you to begin to exercise your taxonomy chops. 
> 
> 
> The data from my paper is much smaller volumes, with far fewer synapses and my detection scheme was aimed at excitatory synapses only.  
> 
> 
> The puncta from Kristina15 I sent you were 3d local maxima of synapsin 1, with the features i spelled out.  Integrated brightness is just the sum of the IF intensity within a small cube. Localized Brightness is the same thing with a 1/x weight attached to it. The distance to the center of mass, is exactly what its described as.  Take the local maxima, take the COM within the cube, calculate the distance. Similar with the moment of inertia. 
> 
> 
> These are spelled out in the following paper
> 
> http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002976
> 
> Automated Analysis of a Diverse Synapse Population - PLOS
> journals.plos.org
> Abstract. Synapses of the mammalian central nervous system are highly diverse in function and molecular composition. Synapse diversity per se may be critical to brain ...
> 
> 
> With respect to the synapse detection in my paper, it was as described in that paragraph... i'm a bit confused as to what is not described there?  We are thresholding, and using connected components above that threshold to detect puncta.   I don't understand how you could use the mean intensity within puncta in order to detect puncta.. that seems a bit circular to me.  
> 
> 
> Once puncta are identified, there is a process of filtering them out.  I describe two different approaches in my paper, a simple version based upon integrated intensity within the puncta, and a more sophisticated one using machine learning.
> 
> 
> Here is my matlab code i used to detect puncta
> 
> 
> function object_file=segment_vstack_2d_v2(targetpath,GPfile,signalthresh)
> 
> currdir=pwd;
> cd(targetpath);
> files=dir('*.tif');
> 
> % Validate input 
> files = validateFolderFiles(files); 
> 
> numfiles=length(files);
> 
> info=imfinfo(files(1).name); % #FIXME_read
> N=info.Height;
> M=info.Width;
> 
> data=zeros(N,M,length(files),'uint8');
> threshdata=logical(data);
> %clear data;
> %%BB={};
> for slice=1:numfiles 
>     frame=uint8(imread(files(slice).name,'tif'));
>     gpFrame=uint8(imread(GPfile,'tif',slice));
>     badpixels=gpFrame==0;
>     frame(badpixels)=0;
>     data(:,:,slice)=frame;
>     %threshdata(:,:,slice)=frame>signalthresh;
>  %%   thisbound=bwboundaries(threshdata(:,:,slice));   
>  %%   BB=[BB;thisbound];
>     disp([slice  length(files)]);
> end
> 
> % 
> %size(threshdata)
> CC=bwconncomp(data>signalthresh,4);
> %clear threshdata
> % %BB=bwboundaries(data>signalthresh,4);
> %stats={};
> stats=regionprops(CC,data,'BoundingBox','MeanIntensity','Area','MaxIntensity','WeightedCentroid');
> % 
>  cd(currdir);
> % 
>  slashes=strfind(targetpath,filesep);
>  channelname=targetpath(slashes(end)+1:end);
> object_file=[targetpath(slashes(end)+1:end) '_objfilev2_thresh' num2str(signalthresh) '.mat'];
>  save(object_file,'CC','stats','targetpath','signalthresh','GPfile','numfiles','channelname');
> 
> hope this helps..
> 
> Forrest
>
>> From: Jesse L. Patsolic <studiojlp@gmail.com>
>> Sent: Monday, May 16, 2016 8:10:25 AM
>> To: Forrest Collman
>> Cc: joshua vogelstein
>> Subject: Questions about puncta detection and features in the Kristina15 data.
>>  
>> Hello Forrest,
>> 
>> Regarding the Kristina15 data set that Joshua and I have been working with for synaptome taxonomy:  
>> 
>> We found the following paragraph in your Journal of Neuroscience paper,
>> 
>> To begin identifying synapses from GABA-negative axons, 2D segmentation of PSD-95 IF was performed by applying a threshold, set to include all IF puncta above the background autofluorescence (see Fig. 9a). For quantitative consistency, this threshold was defined by lowering the threshold until the resulting segmentation produced a median PSD-95 punctum size of > 0.09 µm². We merged 2D PSD-95 puncta from adjacent sections into a single 3D PSD-95 punctum when their IF weighted centroids were < 400 nm apart (see Fig. 9c). The threshold was set at 400 nm, because analysis of the conjugate AT data show that ~90% of all such merges resulted in linking PSD-95 puncta that overlap with the same synapse, without merging distinct synapses (see Fig. 10a).
>> 
>> We are wondering how the puncta were detected; a link to code with a brief explanation would most likely suffice for this.  
>> 
>> What do you use as your threshold: for example, is it the mean intensity within the puncta?
>> 
>> In the paper, I have not found definitions for computing the features
>> 
>> f0 = integrated brightness  
>> f1 = local brightness  
>> f2 = distance to Center of Mass  
>> f3 = moment of inertia around synapsin maxima  
>> f4 = ?  
>> f5 = ?,  
>> 
>> would you be able to point us to a place where those are defined? 
>> 
>> 
>> Thank you.
>> 
>> Best,
>> 
>> Jesse
