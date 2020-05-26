%~~~~~~~~~~~~~ Single-layer TPT-based Traction Force Microscopy ~~~~~~~~~~~~~~~~
%
%Caller function to run each (Matlab) step in the pre-FEniCS part of the
%SL-TPT-TFM process.  In sequence this code: converts a tif stack to .mat,
%crops a single bead image from the data to get a PSF (with user input),
%deconvolves the dataset with this experimental PSF, feeds the deconvolved
%volumes to TPT to compute 3D displacements, regularizes and remeshes TPT
%displacements onto a regular grid using regularizeNd and a smoothing step,
%and output the displacements for the FEniCS finite element computation.  Most
%subfunctions also have tunable parameters that may be changed as necessary.

% ====== Outputs for FEA ======
% x0{multipoint}{timepoint}
% x1{mulitpoint}{timepoint}
% ====== Interface fitting ======
% zPlaneCoeffRef0{m}{t}
% zPlaneCoeffRef1{m}{t}
% ====== TPT tracked beads ======
% x0tracked{m}{t}
% x1tracked{m}{t}
% dispx     \
% dispy      > renamed u_scplane{m}{t}{1-3} (for scattered, planarized)
% dispz     /
% ====== Regular mesh with regularized results ======
% xGrid         \
% yGrid         / renamed gridPts{m}{t}{1-2}
% dispxGrid     \
% dispyGrid      > renamed u_plane{m}{t}{1-3}
% dispzGrid     /
%
% Lauren Hazlett, Alex Landauer, Jin Yang, Mohak Patel
% March, 2020
% Franck Lab, Brown Univerisity and University of Wisc - Madison

clear all

%% Set up input data

%folder to contain data structure (start with a 'tifs' subfolder with images)
data_name = 'traction_error_analysis'; %folder that contains subfolder 'tifs'
save_file_descriptor = 'err_ana_imgs';  % change for EACH day/data set to prevent files from saving over each other


% data_name = 'data'; %default
cur_dir = dir();
data_dir = [cur_dir(1).folder,filesep,data_name,filesep];

multipoint_names = {'force_scaling','width_scaling'}; %if only one multipoint, leave '' in the brackets
cellBWfilename = 'neutrophil_celldata_BW.mat';

%image sequence set up
start_img = {1,1};    %start image in the image squence for each multipoint
total_images = {12,12}; %total number of images to operate over for each multipoint 
z_height = {42,42};    %number of z-slices acquired for each multipoint - used for parsing the .tif

runMode = 'c';      %'i' = incremental, 'c' = cumulative
isCrop = {0,0}; % 1 or 0; if you want to crop tif images, set isCrop to 1 
Icrop = {[1 1 59; 1024 1024 102],[1 1 58; 1024 1024 101]}; % crop format is {[y1 x1 z1; y2 x2 z1], ...} per multipoint. 
isfigCrop = {1}; % 1 or 0; if you want to crop tif images, set isCrop to 1
figCrop = {[350 380; 630 650],[65 50; 980 970]}; % fig crop format is [ymin xmin, ymax xmax] per multipoint

%micrometer to pixel ratios in each dimension, [x,y,z]
um2px = [0.16,0.16,0.5]; %n um per 1 vxl

%gel material properties (linear elastic, NH free energy function)
E = 1500.0;
nu = 0.45;
thickness = 70; %approximate thickness of the gel in um,
%enter 0 for automatic 'thick gel' assumption/meshing
%assume gel thickness is approx constant across multipoints

save([data_dir,data_name,'_',save_file_descriptor,'_','all','_settings.mat'], '-v7.3')

%% Convert tif series/stack to .mat
%     Convert .tif images from the microscope software into variables in a .mat
%     binary container. Comment out or skip if you already have .mat files to use
%     in the data directory.
disp('Running .tif conversion')
matName = cell(length(multipoint_names),1);
for multipoint = 1:length(multipoint_names)
    fprintf('\nWorking on multipoint %i of %i\n',multipoint,length(multipoint_names))

    %image sub-folder name
    image_dir = ['tifs',filesep,multipoint_names{multipoint},filesep];

    matName{multipoint} = cell(total_images{multipoint},1);
    for img_num = start_img{multipoint}:start_img{multipoint}+total_images{multipoint}-1
        matName{multipoint}{img_num} = ...
            funConverTif2Mat(data_dir,image_dir,multipoint_names{multipoint},img_num,z_height{multipoint},isCrop{multipoint},Icrop{multipoint});
        fprintf('Image number: %i of %i completed\n',img_num,total_images{multipoint})
    end

end
disp('Tif convert complete')

save([data_dir,data_name,'_',save_file_descriptor,'_','all','_tif'])

%% Create PSF - use only one per complete experiment (time pts + multi pts)

disp('Running PSF generation')

for multipoint = 1:length(multipoint_names)

    % Grab a window around a bead to define the PSF of the 'scope

    if ~exist('matName', 'var') %Load matNames if using converted matFiles
        for i = 1:length(multipoint_names)
            mat_dir{i} = [data_dir, 'mat files', filesep, multipoint_names{i}, filesep];
            matFiles{i} = dir(fullfile(mat_dir{i}, '*.mat'));
            for j = start_img{multipoint}:start_img{multipoint}+total_images{multipoint}-1
                matName{i}{j} = [mat_dir{i}, matFiles{i}(j).name];
            end
        end
    end

    %predefined center of a bead in xyz format only use if known, (e.g. rerunning an image)
    center_loc  = [0 0 0]; %Enter zero(s) or remove from function call to ignore
    sizePSF = 10;
    [PSFname,PSF] = funCropPSF(matName{multipoint},multipoint_names{multipoint},data_dir,center_loc,sizePSF);

    figure; imagesc3D(PSF) %displays PSF, re-run section is PSF is not acceptable
    fprintf('\nPSF %i of %i complete\n',multipoint,length(multipoint_names))

end

save([data_dir,data_name,'_',save_file_descriptor,'_','all','_tif_psf'])


%% Deconvolve images with the PSF
% Using the PSF and .mat images from above, deconvolve the 3D images
% (using Lucy-Richardson) to reconstruct the point-source locations of the
% fluorecent beads
disp('Running deconvolution')
% load([data_name,'_',save_file_descriptor,'_tif_psf'])
deconvName = cell(length(multipoint_names),1);
img_size = cell(length(multipoint_names),1);
for multipoint = 1:length(multipoint_names)
    fprintf('\nWorking on multipoint %i of %i\n',multipoint,length(multipoint_names))

    if ~exist('matName', 'var') || ~exist('PSF', 'var') %Load matNames and PSF if using converted matFiles with existing PSF
        for i = 1:length(multipoint_names)
            mat_dir{i} = [data_dir, 'mat files', filesep, multipoint_names{i}, filesep];
            matFiles{i} = dir(fullfile(mat_dir{i}, '*.mat'));
            for j = start_img{multipoint}:start_img{multipoint}+total_images{multipoint}-1
                matName{i}{j} = [mat_dir{i}, matFiles{i}(j).name];
            end
            PSF_dir{i} = [data_dir, 'PSF', filesep, multipoint_names{i}, filesep];
            PSFfiles{i} = dir(fullfile(PSF_dir{i},'*_PSF.mat'));
            PSF_{i} = [PSF_dir{i}, PSFfiles{i}.name];
            load(PSF_{i});
        end
    end

    deconvName{multipoint} = cell(total_images{multipoint},1);
    img_size{multipoint} = cell(total_images{multipoint},1);
    %     load(PSF_{multipoint})  %uncomment line if skipping 'Create PSF' section above
    for img_num = start_img{multipoint}:start_img{multipoint}+total_images{multipoint}-1
        [deconvName{multipoint}{img_num},img_size{multipoint}{img_num}] = ...
            funDeconvolve3D(matName{multipoint},data_dir,multipoint_names{multipoint},img_num,PSF);
        fprintf('Deconv image: %i of %i completed\n',img_num,total_images{multipoint})
    end
end
disp('Deconvolution complete')
save([data_dir,data_name,'_',save_file_descriptor,'_','all','_tif_psf_deconv'])

%% Particle localization and tracking
disp('Running localization and particle tracking')
for multipoint = 1:length(multipoint_names)
    

    fprintf('\nWorking on multipoint %i of %i\n',multipoint,length(multipoint_names))

    if ~exist('deconvName', 'var') %Load names of deconvolved files
        for i = 1:length(multipoint_names)
            %input folder within data_dir that contains mat files
            deconv_dir{i} = [data_dir, 'deconv', filesep, multipoint_names{i}, filesep];
            deconvFiles{i} = dir(fullfile(deconv_dir{i},'*.mat'));
            for j = start_img{i}:start_img{i}+total_images{i}-1
                deconvName{i}{j} = [deconv_dir{i}, deconvFiles{i}(j).name];
            end
        end
    end

    % Bead Parameter
    beadParam{1}.thres = 0.125;
    beadParam{1}.minSize = 2;
    beadParam{1}.maxSize = 200;
    beadParam{1}.forloop = 1;

    % TPT Parameters
%     tptParam{1}.knnFM = 10;
%     tptParam{1}.knnFD = 16;
%     tptParam{1}.fmThres = 2;
%     tptParam{1}.outlrThres = 10;
%     tptParam{1}.nSpheres = 2;

    tptParam{1}.knnFM = 10;
    tptParam{1}.knnFD = 16;
    tptParam{1}.fmThres = 2;
    tptParam{1}.outlrThres = 7;
    tptParam{1}.nSpheres = 2;    

    % Check bead parameters
    %     if multipoint == 1
    %         maxhist = 100;
    %         findParams = 1;
    %         while findParams == 1
    %             [beadParam,findParams] = getBeadParams(...
    %                 deconvName{multipoint},maxhist,beadParam);
    %         end
    %     end

    % Track Particles
    [x0{multipoint}, x1{multipoint}, x{multipoint}, track{multipoint}, u{multipoint}] = ...
        funRunTPT(deconvName{multipoint}, beadParam, tptParam, runMode, um2px, multipoint_names{multipoint});

end

disp('Localization and particle tracking complete')
save([data_dir,data_name,'_',save_file_descriptor,'_','localizationtrackingdispresults'], 'x', 'track', 'x0', 'x1', 'u', 'beadParam', 'tptParam')



%% Regularize scattered data and save for FEniCS
disp('Planarizing displacement data')

for multipoint = 1:length(multipoint_names)

    fprintf('\nWorking on multipoint %i of %i\n',multipoint,length(multipoint_names))

    if ~exist('x', 'var') %Load names of deconvolved files
        load([data_dir,data_name,'_',save_file_descriptor,'_','localizationtrackingdispresults'])
    end

    [zPlaneCoeff0{multipoint}, zPlaneCoeff1{multipoint}, x0tracked{multipoint}, ...
        x1tracked{multipoint}, u_scplane{multipoint}, u_plane{multipoint}, gridPts{multipoint}] ...
        = funPlanarizeDispData(x0{multipoint}, x1{multipoint}, x{multipoint},...
        track{multipoint}, u{multipoint}, multipoint_names{multipoint}, runMode);

end

disp('Planarization complete')
save([data_dir,data_name,'_',save_file_descriptor,'_planarizedDataForFEniCS'], ...
    'x0', 'x1', 'u', 'zPlaneCoeff0', 'zPlaneCoeff1', 'x0tracked', 'x1tracked',...
    'u_scplane', 'gridPts', 'u_plane', 'multipoint_names')

%% Plot displacement 

disp('Plotting displacement')

for multipoint = 1:length(multipoint_names)
    
    fprintf('\nWorking on multipoint %i of %i\n',multipoint,length(multipoint_names))
    
        if ~exist('u', 'var') %Load displacement information
            load([data_dir,data_name,'_',save_file_descriptor,'_planarizationdataforFEniCS']);
        end
        
    % For visualization: density and size of cones for displacement coneplot
    density = 6;            % smaller number = higher density
    coneSize = 0.02;       % smaller number = smaller cone size
                 
    load([data_dir,cellBWfilename])
        
    [cellcentroid{multipoint}] = funPlotDisplacement(u_plane{multipoint}, ...
        gridPts{multipoint}, BW{multipoint}, um2px, density, coneSize, ...
        multipoint_names{multipoint}, isfigCrop{multipoint}, figCrop{multipoint}, data_dir);
    
end

disp('Displacement plotting complete')
save([data_dir,data_name,'_',save_file_descriptor,'_displacementplottingdata'], '-v7.3')

%% Output for FEA step
%For each displacement field shift onto the range 0 to max(x,y,z) and
%define scaling factors for FEA mesh manipulation in FEniCS
disp('Saving for FEniCS')
for multipoint = 1:length(multipoint_names)
    fprintf('\nWorking on multipoint %i of %i\n',multipoint,length(multipoint_names))

    %set up Python runscript for the FEniCS call

    load_steps = 10;

    %input and output in the CD for FEniCS I/O
    saveNameIn = [cur_dir(1).folder,filesep,data_name,'_',...
        multipoint_names{multipoint},'_fenicsIn','.mat']; %Input for FEniCS
    saveNameOut = [cur_dir(1).folder,filesep,data_name,'_',...
        multipoint_names{multipoint},'_fenicsOut']; %Output from FEniCS

    for img_num = 1:length(cellcentroid{multipoint})
        x_centroid(img_num) = cellcentroid{multipoint}{img_num}(1);
        y_centroid(img_num) = cellcentroid{multipoint}{img_num}(2);
    end
    x_centroid_norm = mean(x_centroid);
    y_centroid_norm = mean(y_centroid);

%     x_centroid_norm = (x_centroid + min(gridPts{multipoint}{dataset_number}{1}))/...
%         (max(gridPts{multipoint}{dataset_number}{1}) - min(gridPts{multipoint}{dataset_number}{1}));
%     y_centroid_norm = (y_centroid + min(gridPts{multipoint}{dataset_number}{1}))/...
%         (max(gridPts{multipoint}{dataset_number}{1}) - min(gridPts{multipoint}{dataset_number}{1}));
    
    %update the runscript
    updatePyRun(saveNameIn,saveNameOut,multipoint_names{multipoint},...
        E,nu,load_steps,length(u_plane{multipoint}),thickness,x_centroid_norm,y_centroid_norm);

    %save out the data for FEniCS
    for dataset_number = 1:length(u_plane{multipoint})

        %filter edges to blend 'softly' to zeros in a boarder of 1/10th the
        %total with of the image with a blend width (sigma) of 1/50th
        clear dispdata coor
        
        filt_wd = round(size(u_plane{multipoint}{dataset_number}{1},1)/10);
        filt_str = round(filt_wd/5);
        u_x = spatial_rolloff_filt(u_plane{multipoint}{dataset_number}{1},filt_str,filt_wd);
        u_y = spatial_rolloff_filt(u_plane{multipoint}{dataset_number}{2},filt_str,filt_wd);
        u_z = spatial_rolloff_filt(u_plane{multipoint}{dataset_number}{3},filt_str,filt_wd);

        dispdata(:,1) = u_x(:);
        dispdata(:,2) = u_y(:);
        dispdata(:,3) = u_z(:);
        coor(:,1) = um2px(2)*gridPts{multipoint}{dataset_number}{1}(:);
        coor(:,2) = um2px(2)*gridPts{multipoint}{dataset_number}{2}(:);
%       coor(:,3) = gridPts{multipoint}{dataset_number}{3}(:);
        coor(:,3) = ones(size(coor(:,1)))*0;%*thickness;
        save(sprintf('%s_%03d.mat',saveNameIn(1:end-4),dataset_number-1),'coor','dispdata','um2px','-v7')
    end

    figure
    quiver3(coor(:,1),coor(:,2),coor(:,3),dispdata(:,1),dispdata(:,2),dispdata(:,3),1)

    save([data_dir,data_name,'_',multipoint_names{multipoint},'_tif_psf_deconv_localization_tpt_fenics'])

end

disp('Save for FEniCS complete')

