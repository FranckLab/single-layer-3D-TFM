function [savename,PSF] = funCropPSF(folder,imgName, savefolder,center_loc,sizePSF)
% Function to crop out a single bead volume image from a complete volume, to
% be used as an experimental point spread function (PSF) for deconvolution.
% This assumes that a bead is idealized a point source of light in the
% image. PSF size is defined as 1/2 the edge length of a cube centered at
% 'center_loc'
%
%--- INPUTS ---
%  folder    : cell array list of file names for the time series at the
%              current multipoint
%  center_loc: (optional) the predefined center of bead in xyz format
%  sizePSF   : size of the PSF to grab, should be around [10,10,10] for TPT
%              beads (a 21x21x21 bounding box)
%
%--- OUTPUTS ---
% PSF: the cropped-out single bead image to be used for deconvolution
%
% NOTES
% ----------------------------------------------------------------------
% Aug, 2019; Alex Landauer, Mohak Patel, Lauren Hazlett
% Franck Lab, Brown Univerisity and University of Wisc - Madison
%
% If used please cite:
%


i = 1; %use the first data file found, this should have the variable "I",
%the full experimental volume image containing many beads

if nargin < 1
    folder = ['.',filesep,'data',filesep];
    filename = [folder,'*.mat'];
    file = dir(filename);
    if isempty(file) %print warning if no files are found
        fprintf('No .mat files found under: %s \n',folder)
    end

end

if nargin < 4
    sizePSF = 10;
end

%% Set up for loading and saving images
[~,name,~] = fileparts(folder{1});
savedir = [savefolder,'PSF',filesep, imgName, filesep]; %save in subdir "PSF"
if exist(savedir,'dir') ~= 7 %make a new output folder if none exists
    mkdir(savedir);
end
savename = [savedir,name,'_PSF.mat']; %append "PSF" to the output image

%% Load image
I = load(folder{1});
fieldName = fieldnames(I);
I = getfield(I,fieldName{1});
if iscell(I)
    if numel(I)==1, I = I{1};
    else
        I = I{fileInfo.datachannel};
    end
end

%region around the centerpoint where the bead shows up
cropReg = [sizePSF,sizePSF,sizePSF]; %adjust as needed based on relative size of bead

%convert and normalize input image
I = double(I);
I = I/max(I(:));

if ~exist('center_loc','var') || prod(center_loc) == 0
    %note the position of an isolated bead
    figure
    imagesc3D(I)
    drawnow

    %find a single, isolated bead bead and record the center-point
    x_pos = input('enter the x-coord center point: ');
    y_pos = input('enter the y-coord center point: ');
    z_pos = input('enter the z-coord center point: ');

    % The point is in yxz
    pt = [y_pos,x_pos,z_pos];
else
    pt = [center_loc(2),center_loc(1),center_loc(3)];
end
%% Crop PSF

PSF = I(pt(1)-cropReg(1):pt(1)+cropReg(1),pt(2)-cropReg(2):pt(2)+cropReg(2),...
    pt(3)-cropReg(3):pt(3)+cropReg(3));

%normalize PSF
PSF = PSF/max(PSF(:));
close;
%% Save to file
save(savename,'PSF')
