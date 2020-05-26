function [beadParam,findParam] = getBeadParams(deconvName,maxhist,beadParam)

%load the example (first) deconv'd image
load(deconvName{1});

%normalize and binarize
vol = vol/max(vol(:));
BW = vol>beadParam{1}.thres;

CC = bwconncomp(BW);
numPixels = cellfun(@numel,CC.PixelIdxList);
idx = numPixels<3;
idx = find(idx);
for i = 1:length(idx)
    BW(CC.PixelIdxList{idx(i)}) = 0;
end
I = imgaussfilt3(vol,0.75);
Im = -I;
Im(~BW) = Inf;
L = watershed(Im);
L(~BW) = 0;

numPixels = regionprops(L, 'Area');
numPixels = double(struct2dataset(numPixels));

figure;
hist(numPixels(numPixels<maxhist),30)
title('Histogram of detected particle size')
xlabel('Detected particle size (in voxels)')
ylabel('Number of occurances')

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

YN = [];
while isempty(YN)
    YN = input('\nFinalize parameters? (Y/N): ','s');
end

if YN(1) == 'N' ||  YN(1) == 'n'
    findParam = 1;
else
    findParam = 0;
end


%rerun to check thresh param
BW = vol>thresh;
CC = bwconncomp(BW);
numPixels = cellfun(@numel,CC.PixelIdxList);
idx = numPixels<3;
idx = find(idx);
for i = 1:length(idx)
    BW(CC.PixelIdxList{idx(i)}) = 0;
end
I = imgaussfilt3(vol,0.75);
Im = -I;
Im(~BW) = Inf;
L = watershed(Im);
L(~BW) = 0;

numPixels = regionprops(L, 'Area');
numPixels = double(struct2dataset(numPixels));

% figure;
% hist(numPixels(numPixels<maxhist),30)

end