function [x0, x1, x, track, u] = funRunTPT(varargin)
%[x0, x1, x, track, u] = funRunTPT(deconvName, beadParam, tptParam, runMode, um2px)
% Find and localize particle positions in an image and input these to TPT for tracking
%
%--- INPUTS ---
%  deconvName : the list of files to the deconvolved images (store on disk rather
%               than RAM, since the deconv process requires a lot of RAM)
%  beadParam: initial bead parameters to use when localizing particles
%  tptParam: parameters to use when tracking particles
%  runMode: incremental or cumulative displacement (either 'inc' or 'cum')
%  um2px: micrometer to pixel ratio for bead images (format: [x y z])
%
%--- OUTPUTS ---
% x0      : particle positions in reference image
% x1      : particle positions in displaced image
% x       : particle locations before and after displacement
% track   : the output tracking results from TPT
% u       : particle displacements
%
% NOTES
% ----------------------------------------------------------------------
% March, 2020; Lauren Hazlett, Jin Yang, Alex Landauer
% Franck Lab, Brown Univerisity and University of Wisc - Madison
%
% If used please cite: 
%

%% ~~~~~~~~~~~ Locate, localize and track particles in 3D ~~~~~~~~~~~~

%% Parse inputs
[fileInfo, beadParameter, tptParameter, runMode, um2px, mpname] = parseInputs(varargin{:});

for t = 2:length(fileInfo)
    tStart = tic;
    disp(['Current time point: t = ' num2str(t-1)])
    
    % Detect and Localize Particles ---------------------------------
    tPP = tic;
    disp(['  Current Filename: ' fileInfo{t}.name])
    
    if t == 2 % timePoint = 1
        I = loadFile(fileInfo{1});
        x{1}{1} = locateParticles(I,beadParameter{1});
        x{1}{1} = radialcenter3dvec(I,x{1}{1},beadParameter{1});
    end
    
    I = loadFile(fileInfo{t}); %Load image
    x{t}{1} = locateParticles(I,beadParameter{1}); % Detect particles
    x{t}{1} = radialcenter3dvec(I,x{t}{1},beadParameter{1}); % Localize particles
    J{t} = I;
    
    disp(['    Time to localize particles = ', num2str(toc(tPP)),' seconds']);
    disp(['    Number of  particles = ', num2str(size(x{t}{1},1))]);
    
    % Remove outliers -----------------------------------------------
    % -- BLANK --
    
    % Particle tracking ---------------------------------------------
    
    predictor.flag = false;
    predictor.x0 = [];
    predictor.u = [];
    
    tptParameter{1}.sizeI = size(I);
    
    
    %Run TPT
    if runMode(1) == 'i' || runMode(1) == 'I'
        track{t-1}{1} = TPT(x{t-1}{1},x{t}{1},tptParameter{1},predictor);
    elseif runMode(1) == 'c' || runMode(1) == 'C'
        track{t-1}{1} = TPT(x{1}{1},x{t}{1},tptParameter{1},predictor);
    end
    
    disp(['Total Elaspsed Time = ', num2str(toc(tStart)),' seconds']);fprintf('\n');
    
end
%% Compute displacements
if runMode(1) == 'i' || runMode(1) == 'I'
    for i = 1:length(track)
        tempT = track{i}{1};
        x0{i} = [x{i}{1}(tempT>0,:)];
        x1{i} = [x{i+1}{1}(tempT(tempT>0),:)];
        % Compute displacement
        u{i} = x1{i}-x0{i};
        u{i} = u{i}-median(u{i}); %remove rigid drift
        % Convert voxel measurement to um
        u{i}(:,1) = u{i}(:,1) * um2px(1);
        u{i}(:,2) = u{i}(:,2) * um2px(2);
        u{i}(:,3) = u{i}(:,3) * um2px(3);
        % Compute displacement magnitude
        u{i}(:,4) = sqrt(u{i}(:,1).^2 + u{i}(:,2).^2 + u{i}(:,3).^2);
    end
elseif runMode(1) == 'c' || runMode(1) == 'C'
    for i = 1:length(track)
        tempT = track{i}{1};
        x0{i} = [x{1}{1}(tempT>0,:)];
        x1{i} = [x{i+1}{1}(tempT(tempT>0),:)];
        % Compute displacement
        u{i} = x1{i}-x0{i};
        u{i} = u{i}-median(u{i}); %remove rigid drift
        % Convert voxel measurement to um
        u{i}(:,1) = u{i}(:,1) * um2px(1);
        u{i}(:,2) = u{i}(:,2) * um2px(2);
        u{i}(:,3) = u{i}(:,3) * um2px(3);
        % Compute displacement magnitude
        u{i}(:,4) = sqrt(u{i}(:,1).^2 + u{i}(:,2).^2 + u{i}(:,3).^2);
    end
end


% Visualize detected beads overlaid with original image
figure
for i = 1:length(track)
    if runMode(1) == 'i' || runMode(1) == 'I'
        J = loadFile(fileInfo{i}); %Load image
    elseif runMode(1) == 'c' || runMode(1) == 'C'
        if i == 1
            J = loadFile(fileInfo{1}); %Load image
        end
    end
    J_mip = max(J, [], 3);
    subplot(1,length(track),i)
    imshow(J_mip, [])
    hold on
    plot(x0{i}(:,2), x0{i}(:,1), 'go');
    title(['Detected beads overlaid with original image, multipoint: ' mpname, ', timepoint: ',num2str(i)]);
end

figure
for i = 1:length(track)
    if runMode(1) == 'i' || runMode(1) == 'I'
        J = loadFile(fileInfo{i}); %Load image
    elseif runMode(1) == 'c' || runMode(1) == 'C'
        if i == 1
            J = loadFile(fileInfo{1}); %Load image
        end
    end
    J_mip = max(J, [], 3);
    
    % Visualize measured displacement fields
    subplot(1,length(track),i)
    imshow(J_mip, []); axis on; axis image
    hold on
    quiver(x0{i}(:,2), x0{i}(:,1), u{i}(:,2), u{i}(:,1));
    title(['Displacement overlaid on original image, multipoint: ' mpname, ', timepoint: ',num2str(i)]);
end

end


function varargout = parseInputs(varargin)
% varargout = parseInputs(filename, beadParameter, TPTParameter, runMode, um2px)


%%% Parse filenames
filename = varargin{1};

for i = 1:length(filename)
    filedata(i) = dir(fullfile(filename{i}));
    fileInfo{i}.folder = filedata(i).folder;
    fileInfo{i}.name = filedata(i).name;
    loadfile = [filedata(i).folder, filesep, filedata(i).name];
    fileInfo{i}.loadfile = loadfile;
end

%%% Detection and Localization Parameters
beadParameter = varargin{2};

% Define default values
thres = 0.5;
minSize = 5;
maxSize = 100;
winSize = [7,7,7];
dccd = [1,1,1];
abc = [1,1,1];
forloop = 1;
randNoise = 1/10^7; % Something small
xy = 1; % Assume symmetrical if not given
z = 1;  % Assume symmetrical if not givenf
diameter = 5;   % Assume if not given

for i = 1:length(beadParameter)
    
    p = inputParser;
    addParameter(p,'thres',thres);
    addParameter(p,'minSize',minSize);
    addParameter(p,'maxSize',maxSize);
    addParameter(p,'winSize',winSize);
    addParameter(p,'dccd',dccd);
    addParameter(p,'abc',abc);
    addParameter(p,'forloop',forloop);
    addParameter(p,'randNoise',randNoise);
    addParameter(p,'xy',xy);
    addParameter(p,'z',z);
    addParameter(p,'diameter',diameter);
    
    parse(p,beadParameter{i})
    
    beadParameter{i} = p.Results;
    
end

%%% TPT Parameters
tptParameter = varargin{3};

% Define default values
knnFD =16;
knnFM = 10;
fmThres = 2;
maxIter = 14;
nSpheres = 2;
outlrThres = 10;

for i = 1:length(beadParameter)
    
    p = inputParser;
    addParameter(p,'knnFD',knnFD);
    addParameter(p,'knnFM',knnFM);
    addParameter(p,'fmThres',fmThres);
    addParameter(p,'maxIter',maxIter);
    addParameter(p,'nSpheres',nSpheres);
    addParameter(p,'outlrThres',outlrThres);
    
    parse(p,tptParameter{i})
    
    tptParameter{i} = p.Results;
    
end

runMode = varargin{4};
um2px = varargin{5};
mpname = varargin{6};


%%% Outputs
varargout{1} = fileInfo;
varargout{2} = beadParameter;
varargout{3} = tptParameter;
varargout{4} = runMode;
varargout{5} = um2px;
varargout{6} = mpname;

end

function I = loadFile(fileInfo)
I = load(fileInfo.loadfile);
I = I.vol;
I = double(I);
I = I/max(I(:));
end
