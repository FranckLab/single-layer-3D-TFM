%~~~~~~~~~~~~~ Single-layer TPT-based Traction Force Microscopy ~~~~~~~~~~~~~~~~
%
% Caller function to run each (Matlab) step in the pre-FEniCS part of the
% SL-TPT-TFM process.
% In sequence this script:
%   - converts a tif stack to .mat
%   - crops a single bead image from the data to get a PSF (with user input)
%   - deconvolves the dataset with this experimental PSF
%   - feeds the deconvolved volumes to TPT to compute 3D displacements
%   - regularizes and remeshes TPT displacements using regularizeNd
%   - outputs the displacements for the FEniCS finite element computation
%
% Most subfunctions also have tunable parameters that may be changed as
% necessary. the most important of these are prompted for during the run,
% but other adjustable parameters may be tunable if result are not saticfactory
%
%Set up the INPUT SETUP block with the information for your dataset
%
% Lauren Hazlett, Alex Landauer, Mohak Patel, Hadley Witt, Jin Yang
% June, 2020
% Franck Lab and Reichner Lab; Brown Univerisity, RI Hospital, and UW - Madison

clear all

%% --------------------------- INPUT SET UP ------------------------------------
data_name = 'example_data'; %data folder that contains subfolder 'tifs', with tif stacks stored in multipoint subfolders
save_file_descriptor = 'example_saved_data';  % change for EACH day/data set to prevent files from saving over each other

cur_dir = dir();
data_dir = [cur_dir(1).folder,filesep,data_name,filesep];

multipoint_names = {'XYpoint_name_1','XYpoint_name_2'}; %if only one multipoint, leave '' in the brackets
cellBWfilename = 'exampleBWfile.mat';


%% ------------------------ IMAGE SEQUENCE SET UP ------------------------------
start_img = {1,1};    %start image in the image squence for each multipoint
total_images = {3,3}; %total number of images to operate over for each multipoint
z_height = {42,42};    %number of z-slices acquired for each multipoint - used for parsing the .tif

runMode = 'c';      %'i' = incremental, 'c' = cumulative
isCrop = {0,0}; % 1 or 0; if you want to crop tif images, set isCrop to 1
Icrop = {[256, 256 1; 767 767 40],[256, 256 1; 767 767 40]}; % crop format is {[y1 x1 z1; y2 x2 z2], ...} per multipoint.
isfigCrop = {0,0}; % 1 or 0; if you want to crop tif images, set isCrop to 1
figCrop = {[256 256; 767 767],[256 256; 767 767]}; % fig crop format is [ymin xmin, ymax xmax] per multipoint

%micrometer to pixel ratios in each dimension, [x,y,z]
um2px = [0.16,0.16,0.5]; %n um per 1 vxl

% gel properties set up
%gel material properties (linear elastic, NH free energy function)
E = 1500.0;
nu = 0.45;
thickness = 70; %approximate thickness of the gel in um,
%enter 0 for automatic 'thick gel' assumption/meshing
%assume gel thickness is approx constant across multipoints

save([data_dir,data_name,'_',save_file_descriptor,'_','all','_settings.mat'], '-v7.3')


%% -------------------- CONVERT TIFs STACKS TO .MAT ----------------------------
%    Convert .tif images from the microscope software into variables in a .mat
%    binary container. Comment out or skip if you already have .mat files to use
%    in the data directory.

disp('Running .tif conversion')
matName = cell(length(multipoint_names),1);
for multipoint = 1:length(multipoint_names)
    fprintf('\nWorking on multipoint %i of %i\n',multipoint,length(multipoint_names))

    %image sub-folder name
    image_dir = ['tifs',filesep,multipoint_names{multipoint},filesep];

    matName{multipoint} = cell(total_images{multipoint},1);
    for img_num = start_img{multipoint}:start_img{multipoint}+total_images{multipoint}-1
        matName{multipoint}{img_num} = ...
            funConverTif2Mat(data_dir,image_dir,multipoint_names{multipoint},img_num,...
                z_height{multipoint},isCrop{multipoint},Icrop{multipoint});
        fprintf('Image number: %i of %i completed\n',img_num,total_images{multipoint})
    end

end
disp('Tif conversion complete')

save([data_dir,data_name,'_',save_file_descriptor,'_','all','_tif'])

%% ----------------------------- GENERATE PSF ----------------------------------
%   Crop a single bead image out of the reference (first) volume to use as a
%   point-source estimator of the PSF. Use only one per multipoint, if the PSF
%   changes between timepoints this likely indicates an problem with data collection

disp('Running PSF generation')
for multipoint = 1:length(multipoint_names)

    % Grab a window around a bead to define the PSF of the microscope
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

    figure
    imagesc3D(PSF) %displays PSF, re-run section is PSF is not acceptable
    fprintf('\nPSF %i of %i complete\n',multipoint,length(multipoint_names))

end
save([data_dir,data_name,'_',save_file_descriptor,'_','all','_tif_psf'])
disp('PSF generation complete')

%% ------------------ DECONVOLVE IMAGES WITH THE PSF ---------------------------
%   Using the PSF and .mat images from above, deconvolve the 3D images
%   (using Lucy-Richardson) to reconstruct the point-source locations of the
%   fluorecent beads

disp('Running deconvolution')
% load([data_name,'_',save_file_descriptor,'_tif_psf']) %uncomment to re-run the block
deconvName = cell(length(multipoint_names),1);
img_size = cell(length(multipoint_names),1);
for multipoint = 1:length(multipoint_names)
    fprintf('\nWorking on multipoint %i of %i\n',multipoint,length(multipoint_names))

    if ~exist('matName', 'var') || ~exist('PSF', 'var') %Load if using converted matFiles with existing PSF
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
save([data_dir,data_name,'_',save_file_descriptor,'_','all','_tif_psf_deconv'])
disp('Deconvolution complete')

%% ------------------ PARTICLE LOCALIZATION AND TRACKING -----------------------
%   Use the tools from the Topology-based Particle Tracking (TPT) package to
%   locate, localize, and track particles (a.k.a beads) in the deconvolved volumes.

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
    beadParam{multipoint}.thres = 0.125;
    beadParam{multipoint}.minSize = 2;
    beadParam{multipoint}.maxSize = 200;
    beadParam{multipoint}.forloop = 1;

    % TPT Parameters
    tptParam{1}.knnFM = 10;
    tptParam{1}.knnFD = 16;
    tptParam{1}.fmThres = 2;
    tptParam{1}.outlrThres = 7;
    tptParam{1}.nSpheres = 2;

    % Check (and refine as needed) bead parameters
    if multipoint == 1
       maxhist = 100;
       findParams = 1;
       while findParams == 1
         [beadParam{multipoint},findParams] = getBeadParams(deconvName{multipoint},maxhist,beadParam{multipoint});
       end
    end
    
    cur_beadParam{1} = beadParam{multipoint};
    % Track Particles with TPT
    [x0{multipoint}, x1{multipoint}, x{multipoint}, track{multipoint}, u{multipoint}] = ...
        funRunTPT(deconvName{multipoint}, cur_beadParam, tptParam, runMode, um2px, multipoint_names{multipoint});

end

save([data_dir,data_name,'_',save_file_descriptor,'_','localizationtrackingdispresults'], 'x',...
    'track', 'x0', 'x1', 'u', 'beadParam', 'tptParam')
disp('Localization and particle tracking complete')


%% --------------------- REGULARIZE SCATTERED DATA -----------------------------
%   Use an nD regularizer to smooth tracked displacement data and construct
%   a gridded data representation from the scattered point given by TPT.
%   Resample on a single best-fit plane to the single-layer of beads in the
%   deconvolved volume

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
save([data_dir,data_name,'_',save_file_descriptor,'_planarizedDataForFEniCS'], ...
    'x0', 'x1', 'u', 'zPlaneCoeff0', 'zPlaneCoeff1', 'x0tracked', 'x1tracked',...
    'u_scplane', 'gridPts', 'u_plane', 'multipoint_names')

disp('Planarization complete')


%% ---------------------------- PLOT DISPLACEMENTS -----------------------------
%   Plot the displacement, cell outline and images on a single set of axes to
%   visualize the tracking result.  Also, use the BW image to compute the cell
%   centroid to estimate the best location for mesh refinement (roughly)

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
save([data_dir,data_name,'_',save_file_descriptor,'_displacementplottingdata'], '-v7.3')
disp('Displacement plotting complete')


%% --------------------------- OUTPUT FOR FEA ----------------------------------
%   Save out displacement data for FEniCS and update .py runscripts accordingly
%   The runscripts and data are save in the current working directory - make
%   sure only the latest version of each that you wish to run is included

%for reference:
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
        coor(:,1) = um2px(1)*gridPts{multipoint}{dataset_number}{1}(:);
        coor(:,2) = um2px(2)*gridPts{multipoint}{dataset_number}{2}(:);
%       coor(:,3) = gridPts{multipoint}{dataset_number}{3}(:);
        coor(:,3) = zeros(size(coor(:,1)));
        save(sprintf('%s_%03d.mat',saveNameIn(1:end-4),dataset_number-1),'coor','dispdata','um2px','-v7')
    end

    % uncomment to show a diagnostic plot fo the displacments as saved for FEniCS
    % figure
    % quiver3(coor(:,1),coor(:,2),coor(:,3),dispdata(:,1),dispdata(:,2),dispdata(:,3),1)

    save([data_dir,data_name,'_',multipoint_names{multipoint},'_tif_psf_deconv_localization_tpt_fenics'])

end

disp('Save for FEniCS complete')

% ------------------------------ END OF FILE -----------------------------------
