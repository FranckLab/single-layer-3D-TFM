function [x0, x1, x, track, u] = funRunTPT(varargin)

%[x0, x1, x, track, u] = funRunTPT(deconvName, beadParam, tptParam, runMode, um2px)
%~~~~~~~~~~~~~~~~~ TPT-based Traction Force Microscopy ~~~~~~~~~~~~~~~~
%
%%Find particle position in the image based on a watershed method. 
%
%--- INPUTS ---
%  filename : the file path to the data (image) storage folder (default:
%               "./data/"), inlcuding file name
%  beadParam: the bead parameters to use when localizing particles
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
% March, 2020; Lauren Hazlett, Jin Yang
% Franck Lab, Brown Univerisity and University of Wisc - Madison

%% Parse inputs
[fileInfo, beadParameter, tptParameter, runMode, um2px, mpname] = parseInputs(varargin{:});

for t = 2: length(fileInfo)
    tStart = tic;
    disp(['Current time point: t = ' num2str(t-1)])
    
% Detect and Localize Particles ---------------------------------
    tPP = tic;
    disp(['  Current Filename: ' fileInfo{t}.name]) 
    
    if t == 2 % timePoint = 1
        I{t} = loadFile(fileInfo{1},beadParameter{1}.randNoise);
        x{1}{1} = locateParticles(I{t},beadParameter{1});
        x{1}{1} = radialcenter3dvec(I{t},x{1}{1},beadParameter{1});
    end
    
    I = loadFile(fileInfo{t},beadParameter{1}.randNoise); %Load image 
    x{t}{1} = locateParticles(I,beadParameter{1}); % Detect particles
    x{t}{1} = radialcenter3dvec(I,x{t}{1},beadParameter{1}); % Localize particles
    J{t} = I;
    
    disp(['    Time to localize particles = ', num2str(toc(tPP)),' seconds']);
    disp(['    Number of  particles = ', num2str(size(x{t}{1},1))]);
   
% Remove outliers -----------------------------------------------


% Particle tracking ---------------------------------------------

    predictor.flag = false;
    predictor.x0 = [];
    predictor.u = [];
    
    tptParameter{1}.sizeI = size(I);


%Run TPT
    if runMode == 'i' || runMode == 'I'
        track{t-1}{1} = TPT(x{t-1}{1},x{t}{1},tptParameter{1},predictor);
    elseif runMode == 'c' || runMode == 'C'
        track{t-1}{1} = TPT(x{1}{1},x{t}{1},tptParameter{1},predictor);
    end
       
    disp(['Total Elaspsed Time = ', num2str(toc(tStart)),' seconds']);fprintf('\n');

end
%% Compute displacements
if runMode == 'i' || runMode == 'I'
    for i = 1:length(track)
        tempT = track{i}{1};
        x0{i} = [x{i}{1}(tempT>0,:)];
        x1{i} = [x{i+1}{1}(tempT(tempT>0),:)];
        % Compute displacement
        u{i} = x1{i}-x0{i};
        u{i} = u{i}-median(u{i});
        % Convert voxel measurement to um
        u{i}(:,1) = u{i}(:,1) * um2px(1);
        u{i}(:,2) = u{i}(:,2) * um2px(2);
        u{i}(:,3) = u{i}(:,3) * um2px(3);
        % Compute displacement magnitude
        u{i}(:,4) = sqrt(u{i}(:,1).^2 + u{i}(:,2).^2 + u{i}(:,3).^2);
    end
elseif runMode == 'c' || runMode == 'C'
    for i = 1:length(track)
        tempT = track{i}{1};
        x0{i} = [x{1}{1}(tempT>0,:)];
        x1{i} = [x{i+1}{1}(tempT(tempT>0),:)];
        % Compute displacement
        u{i} = x1{i}-x0{i};
        u{i} = u{i}-median(u{i});
        % Convert voxel measurement to um
        u{i}(:,1) = u{i}(:,1) * um2px(1);
        u{i}(:,2) = u{i}(:,2) * um2px(2);
        u{i}(:,3) = u{i}(:,3) * um2px(3);
        % Compute displacement magnitude
        u{i}(:,4) = sqrt(u{i}(:,1).^2 + u{i}(:,2).^2 + u{i}(:,3).^2);
    end
end

    
    % Visualize detected beads overlaid with original image
for i = 1:length(track)
    if runMode == 'i' || runMode == 'I'
        J = loadFile(fileInfo{i},beadParameter{1}.randNoise); %Load image 
    elseif runMode == 'c' || runMode == 'C'
        if i == 1
        J = loadFile(fileInfo{1},beadParameter{1}.randNoise); %Load image
        end
    end
    J_mip = max(J, [], 3);
        figure; imshow(J_mip, [])
        hold on
        plot(x0{i}(:,2), x0{i}(:,1), 'go');
        title(['Detected beads overlaid with original image, multipoint: ' mpname, ', timepoint: ',num2str(i)]);

        % Visualize measured displacement fields
        figure; imshow(J_mip, []); axis on
        hold on
        quiver(x0{i}(:,2), x0{i}(:,1), u{i}(:,2), u{i}(:,1));
        title(['Measured displacement field overlaid with original image, multipoint: ' mpname, ', timepoint: ',num2str(i)]);

%         figure; scatter3(x0{i}(:,2), x0{i}(:,1), x0{i}(:,3));  
%         hold on; 
%         quiver3(x0{i}(:,2), x0{i}(:,1), x0{i}(:,3), u{i}(:,2), u{i}(:,1), u{i}(:,3)); 
%         title(['Measured 3D displacement field overlaid with 3D particle positions, multipoint: ' mpname, ', timepoint: ',num2str(i)]);

end

end


function varargout = parseInputs(varargin)
% varargout = parseInputs(filename,beadParameter, TPTParameter, runMode, um2px)


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

function I = loadFile(fileInfo,randNoise)
I = load(fileInfo.loadfile);
I = I.vol;
I = double(I);
I = I/max(I(:));
end
