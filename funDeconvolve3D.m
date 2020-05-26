function [savename,sizeI] = funDeconvolve3D(folder,savefolder,imgName, imgNum,PSF)
%~~~~~~~~~~~~~ Single-layer TPT-based Traction Force Microscopy ~~~~~~~~~~~~~~~~
%
%Function to deconvolve a volume image using the experimentally gathered
%PSF from the previous step (funCropPSF).  More sophisticated deconvolution
%technique are available (e.g. see the "deconvlucy" documentation in Matlab)
%but we have found this to be sufficent
%
%--- INPUTS ---
%  folder: the path to the data storage folder (default: "./data/")
%  savefolder: the path to the main folder where deconvolved data will be
%              stored
%   imgName:   path to correct images, including image names
%   imgNum:    number of image of interest within data structure 
%   PSF   : (optional) the PSF to use for deconvolution, if not given the
%          script will attempt to load PSF.mat from the cd
%
%--- OUTPUTS ---
%   savename:   name of deconvolved images 
%   sizeI:      final image size
%
% June, 2019; Mohak Patel, Alex Landauer, Lauren Hazlett
% Franck Lab, Brown Univerisity and University of Wisc - Madison

%% Set default params
if nargin < 1
    folder = ['.',filesep,'data',filesep];
    filename = [folder,'*.mat'];
    file = dir(filename);
    if isempty(file) %print warning if no files are found
        fprintf('No .mat files found under: %s \n',folder)
    end
end
if nargin < 2
    imgNum = 1;
end
if nargin < 3
    load(['.',filesep,'PSF.mat'],'PSF');
end

deconv_iter = 10; %number of deconv iterations to use
prefilter = false; %true/false gaussian prefilter option

%% Set up for loading and saving images
[~,name,~] = fileparts(folder{imgNum});
savedir = [savefolder,'deconv',filesep, imgName, filesep]; %save in subdir "deconv"
if exist(savedir,'dir') ~= 7 %make a new output folder if none exists
    mkdir(savedir);
end
savename = [savedir,name,'_deconv.mat']; %append "deconv" to the output image


%% Load image
I = load(folder{imgNum});
fieldName = fieldnames(I);
I = getfield(I,fieldName{1});
if iscell(I)
    if numel(I)==1, I = I{1};
    else
        I = I{fileInfo.datachannel};
    end
end


%Do the deconvolution
vol = deconvlucy(I, PSF, deconv_iter);
% [vol,psf] = deconvblind(I, PSF, deconv_iter); %another possible deconv technique
% vol = deconvreg(I, PSF); %a third possible deconv

sizeI = size(vol); %record the final image size

figure; imagesc3D(vol);

%% Save out the result
% filename = file(img_num).name;
save(savename, 'name', 'vol', '-v7.3');

% disp(['Deconv: ',num2str(img_num)])

end

