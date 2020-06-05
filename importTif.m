function [Vol] = importTif(filename,nSlice,alternateImg,saveFlag)
% Import multi-channel OR multi-time point images. This function does not
% handle multi-channel AND multi-time point tiff files. Seperate them in
% FIJI into separate files.

%--- INPUTS ---
% filename      : Tif filename
% nSlice        : Number of z slices for the 3D image for a single channel
% alternateImg  : Defines how images are stored in the tif files
%               If alternateImg == 0
%               111111111222222222333333333
%               If alternateImg == 1
%               123123123123123123123123123
% saveFlag      : (1) Save file or (0) do not save file
%
%--- OUTPUTS ---
% vol           : stores the file in a cell array
%                 vol{fileNumber}{Channel or time point}

%% ~~~~~~~~~~~ Read in stacked tif images ~~~~~~~~~~~~

imagefiles = dir(filename);
nfiles = length(imagefiles);    % Number of files found
Vol = cell(nfiles,1);



for j=1:nfiles

    %Import the whole image
    filename = fullfile([imagefiles(j).folder,filesep,imagefiles(j).name]);

    nImages = numel(imfinfo(filename))/nSlice;
    sizeI = [size(imread(filename,1)),nSlice];

    Vol{j} = cell(nImages,1);

    idx = 1:nSlice;
    for i = 1:nImages
        Vol{j}{i} = zeros(sizeI);
    end

    for i = 1:nImages

        % Read through images in sequence
        if alternateImg == 0

            for k =1:length(idx)
                Vol{j}{i}(:,:,k) = imread(filename,idx(k));
            end
            idx = idx+nSlice;
        end

        % Read through alternate image
        if alternateImg == 1
            for k = 1:nSlice
                Vol{j}{i}(:,:,k) = imread(filename,2*k-2+i);
            end
        end

    end


    %Save file
    if saveFlag == 1
        [~,name,~] = fileparts(filename);
        name = [name '.mat'];

        vol = Vol{j};
        save(name,'vol','-v7.3');
    end
end

end
