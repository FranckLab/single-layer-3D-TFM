function [beadParam,findParam] = getBeadParams(deconvName,maxhist,beadParam)
% Preprocess images with the particle locator and localization to to test and
% optimize bear parameters for TPT. WARNING: TPT can be sensitive to variations
% in bead images - this step is particularly important for lower-quality images
%
%--- INPUTS ---
%  deconvName : the list of files to the deconvolved images (store on disk rather
%               than RAM, since the deconv process requires a lot of RAM)
%  maxhist    : Max size (# of voxels) of beads shown in histogram
%  beamParam  : test set of parameters to use
%
%--- OUTPUTS ---
% beamParam : final set of parameters to used
% findParam : flag indicating if found parameters are okay (0 if final, 1 to rerun)
%
% NOTES
% ----------------------------------------------------------------------
% Nov, 2019; Alex Landauer, Hadley Wiit, Lauren Hazlett, Mohak Patel
% Franck Lab, Brown Univerisity and University of Wisc - Madison
%
% If used please cite:
%

%% ~~~~~~~~~~~ Compute bead stats and prompt for new params ~~~~~~~~~~~~


%load the example (first) deconvolved image
load(deconvName{1});

% Parameters
minPixels = beadParameter.minSize;  %Minimum pixel count in blob for bead
maxPixels = beadParameter.maxSize;  %Maximum pixel count in blob for bead

%normalize and binarize (input I in locateParticles is similarly normalized)
vol = vol/max(vol(:));
BW = vol>beadParam{1}.thres;

%use same process as the locateParticles script
CC = bwconncomp(BW);
numPixels = cellfun(@numel,CC.PixelIdxList);
beadBlob = numPixels>minPixels & numPixels<maxPixels;

numPixels = regionprops(beadBlob, 'Area');
numPixels = double(struct2dataset(numPixels));

%show the resultant bead stats
figure;
hist(numPixels(numPixels<maxhist),30)
title('Histogram of detected particle size')
xlabel('Detected particle size (in voxels)')
ylabel('Number of occurances')

%prompt for updated bead parameters: threshold, min radius, and max radius
disp('Enter updated bead parameters; leave blank to keep same value')
thresh = input('Enter new threshold: ');
min_r = input('Enter new min bead size: ');
max_r = input('Enter new max bead size: ');
disp('Recomputing bead identification... ')
close;

if ~isempty(thresh)
    beadParam{1}.thres = thresh;
else
    thresh = beadParam{1}.thres;
end
if ~isempty(min_r)
    beadParam{1}.minSize = min_r;
end
if ~isempty(max_r)
    beadParam{1}.maxSize = max_r;
end


%rerun to check thresh param
BW = vol>beadParam{1}.thres;

%use same process as the locateParticles script
CC = bwconncomp(BW);
numPixels = cellfun(@numel,CC.PixelIdxList);
beadBlob = numPixels>minPixels & numPixels<maxPixels;

numPixels = regionprops(beadBlob, 'Area');
numPixels = double(struct2dataset(numPixels));

%show the resultant bead stats
figure;
hist(numPixels(numPixels<maxhist),30)
title('Histogram of detected particle size')
xlabel('Detected particle size (in voxels)')
ylabel('Number of occurances')
drawnow

YN = [];
while isempty(YN)
    YN = input('\nFinalize parameters? (Y/N): ','s');
end

if YN(1) == 'N' ||  YN(1) == 'n'
    findParam = 1;
else
    findParam = 0;
end

end
