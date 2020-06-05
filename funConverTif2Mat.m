function [savename] = funConverTif2Mat(folder,imgDir,imgName,imgNum,zHeight,isCrop,Icrop)
% Function to read in a image stack in .tif format and save to a .mat
% datafile for further processing (see importTiff for more details)
%
%--- INPUTS ---
%  folder : the path to the data storage folder (default: "./data/tifs")
%  imgDir : subdirectory for raw image data, (default: "tifs/")
%  imgName: name of files to be used for creating image names (default:
%  blank)
%  imgNum : image number in the data directory
%  zHeight: the number of z-slice for stack
%  isCrop: binary flag to crop (1) or not (0)
%  Icrop: corner locations for cropping
%
%--- OUTPUTS ---
% savename: the path and name of the .mat data file
%
% NOTES
% ----------------------------------------------------------------------
% Aug, 2019; Alex Landauer, Mohak Patel, Lauren Hazlett
% Franck Lab, Brown Univerisity and University of Wisc - Madison
%
% If used please cite:
%

%% Set default params
alternateImg = 0; %How images are stored in the tif, 1 => "123123"

if nargin < 1
    folder = ['.',filesep,'data',filesep,'tifs',filesep];
end
if nargin < 2
    imgDir = 'tifs';
end
if nargin <3
    imgName = '';
end
if nargin < 4
    imgNum = 1;
end
if nargin < 5
   zHeight = 101;
end
if nargin < 6
   isCrop = 0;
end
if nargin < 7
   Icrop = [768,568,10; 1791,1591,40];
end

%% Set up for loading and saving images
filename = [folder,imgDir,'*tif'];
file = dir(filename);
if isempty(file) %print warning if no files are found
    fprintf('No .tif files found under: %s \n',[folder,imgDir])
end
[~,name,~] = fileparts(file(imgNum).name);
savedir = [folder,'mat files',filesep,imgName,filesep]; %save in subdir "mat files"-->multipoint
if exist(savedir,'dir') ~= 7 %make a new output folder if none exists
    mkdir(savedir);
end
savename = [savedir,name,'.mat'];


%% Read in the tif stack
[vol] = importTif(fullfile([file(imgNum).folder,filesep,file(imgNum).name]), zHeight, alternateImg, 0);

%% Save the volume to a .mat
% Crop file
if isCrop
    vol = vol{1}{1}(Icrop(1,2):Icrop(2,2), Icrop(1,1):Icrop(2,1), ...
        Icrop(1,3):Icrop(2,3));
else
    vol = vol{1}{1}; %only save image 1, channel 1
end

save(savename,'vol','-v7.3');

% disp(['Img num: ',num2str(imgNum)])
% disp('Tif convert complete')

end
