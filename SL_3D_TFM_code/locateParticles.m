function [x] = locateParticles(I, beadParameter)
% [x] = locateParticles(I, beadParameter) locates particles in the image using
% a global fixed thresholding-based blob detection
%
% INPUTS
% -------------------------------------------------------------------------
%   I:              Input volumetric image
%   beadParameter:  Parameters to detect particle position in images
%
% OUTPUTS
% -------------------------------------------------------------------------
%   x:              Voxel-level estimate of particle center in MxNxO format
%

% Parameters
thres = beadParameter.thres;    %Threshold value
minPixels = beadParameter.minSize;  %Minimum pixel count in blob for bead
maxPixels = beadParameter.maxSize;  %Maximum pixel count in blob for bead

% Image thresholding
BW = I>thres;

% Find bead blobs
CC = bwconncomp(BW);
numPixels = cellfun(@numel,CC.PixelIdxList);
beadBlob = numPixels>minPixels & numPixels<maxPixels;

% Find centroid of all the eligible blobs
S = regionprops(CC,'Centroid');
BW2dProj = sum(BW,3);

% ====== JY !!! Detect bead eccentricity ======
try
    S2 = regionprops3(CC,'Centroid','Area','Eccentricity');
    S2EquatorialEccentricity = [];
    S2MeridionalEccentricity = [];
    for tempi = 1:length(S2)
        S2EquatorialEccentricity(tempi) = S2(tempi).EquatorialEccentricity;
        S2MeridionalEccentricity(tempi) = S2(tempi).MeridionalEccentricity;
    end

    thrs1 = 0; thrs2 = 0; % 0.4, 0.7
    beadBlob2 = S2EquatorialEccentricity > thrs1 & S2MeridionalEccentricity > thrs2;
    blobPts = round(double(struct2dataset(S)));
    blobPts = blobPts(and(beadBlob,1*beadBlob2) ,:);
catch
    blobPts = round(double(struct2dataset(S)));
    blobPts = blobPts(beadBlob ,:);
end

temp = blobPts;

% Convert to m,n,o coordinates
blobPts(:,1) = temp(:,2);
blobPts(:,2) = temp(:,1);
x = blobPts;

end
