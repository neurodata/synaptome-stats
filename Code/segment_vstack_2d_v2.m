%%%
%%% Matlab function for puncta detection from Forrest Collman
%%% in email correspondence with J. Patsolic and J. T. Vogelstein 20160517
%%%

function object_file=segment_vstack_2d_v2(targetpath,GPfile,signalthresh)

currdir=pwd;
cd(targetpath);
files=dir('*.tif');

% Validate input 
files = validateFolderFiles(files); 

numfiles=length(files);

info=imfinfo(files(1).name); % #FIXME_read
N=info.Height;
M=info.Width;

data=zeros(N,M,length(files),'uint8');
threshdata=logical(data);
%clear data;
%%BB={};
for slice=1:numfiles 
    frame=uint8(imread(files(slice).name,'tif'));
    gpFrame=uint8(imread(GPfile,'tif',slice));
    badpixels=gpFrame==0;
    frame(badpixels)=0;
    data(:,:,slice)=frame;
    %threshdata(:,:,slice)=frame>signalthresh;
 %%   thisbound=bwboundaries(threshdata(:,:,slice));   
 %%   BB=[BB;thisbound];
    disp([slice  length(files)]);
end

% 
%size(threshdata)
CC=bwconncomp(data>signalthresh,4);
%clear threshdata
% %BB=bwboundaries(data>signalthresh,4);
%stats={};
stats=regionprops(CC,data,'BoundingBox','MeanIntensity','Area','MaxIntensity','WeightedCentroid');
% 
 cd(currdir);
% 
 slashes=strfind(targetpath,filesep);
 channelname=targetpath(slashes(end)+1:end);
object_file=[targetpath(slashes(end)+1:end) '_objfilev2_thresh' num2str(signalthresh) '.mat'];
 save(object_file,'CC','stats','targetpath','signalthresh','GPfile','numfiles','channelname');

