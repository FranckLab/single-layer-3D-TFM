function [cellcentroid] = funPlotTraction(trac, points, BW, density, coneSize, mpname, isfigCrop, figCrop, savefolder)
% Plots and saves figures of a cone plot of FE-computed traction data using
% colormap Turbo (Copyright 2020, Daniel Fortunato. See:
%  https://www.mathworks.com/matlabcentral/fileexchange/74662-turbo) along
%  with the cell.
%
%
%--- INPUTS ---
% trac      : surface traction 3D vector components per each time stack.
%             Output format - trac{timePoint}{component}. trac{t}{4} - magnitude
% points    : locations at which traction in sampled
% BW        : 3D matrix of black and white (double format) cell image
% density   : density of cones for coneplot
% coneSize 	: size of cones for coneplot
% mpname    : name of multipoint data
% isfigCrop : 0 or 1 to crop or not crop the data
% figCrop   : [xmin xmax; ymin ymax] format borders for figure cropping
% saveFolder: output folder to which to save figures
%
%
%--- OUTPUTS ---
% cellcentroid  : location of cell centroid from BW image
%
%
% NOTES
% ----------------------------------------------------------------------
% May, 2020; Lauren Hazlett, Alex Landauer, Mohak Patel
% Franck Lab, Brown Univerisity and University of Wisc - Madison
%
% If used please cite:
%

%% ~~~~~~~~~~~ Plot Traction data ~~~~~~~~~~~~

%% Format Displacements and produce cone plot
% folder = ['.',filesep,'data',filesep];
% [~,name,~] = fileparts(folder{1});
savedir = [savefolder,'Figures',filesep, mpname, filesep]; %save in subdir "PSF"
if exist(savedir,'dir') ~= 7 %make a new output folder if none exists
    mkdir(savedir);
end

for timePt = 1 : length(trac)
    fprintf('\nPlotting timepoint %i of %i\n',timePt,length(trac))

    % Initialize cell border information
    im{timePt} = squeeze(sum(BW{timePt},3));
    im{timePt} = im{timePt}/max(im{timePt}(:));
    sizeBW = size(BW{timePt});

    if isfigCrop == 1
        xmin = figCrop(3);
        ymin = figCrop(1);
        xmax = figCrop(4);
        ymax = figCrop(2);
    elseif isfigCrop == 0
        xmin = 1;
        ymin = 1;
        xmax = sizeBW(1);
        ymax = sizeBW(2);
    end

    BW_mip{timePt} = max(BW{timePt},[],3);
    cc = bwconncomp(BW_mip{timePt});
    BW_centroid = regionprops(cc, {'centroid'});
    BW_centroid = BW_centroid.Centroid;
    cellcentroid{timePt} = [BW_centroid(1)/sizeBW(1) BW_centroid(2)/sizeBW(2)];

    % Initializie traction information
    x(:,1) = points{timePt}{1};
    x(:,2) = points{timePt}{2};
    x(:,3) = ones(size(x(:,2)));
    ti(:,1) = trac{timePt}{1};
    ti(:,2) = trac{timePt}{2};
    ti(:,3) = trac{timePt}{3};
%     ti(:,4) = sqrt(ti(:,1).^2 + ti(:,2).^2 + ti(:,3).^2);

    n{1} = 1:sizeBW(1);
    n{2} = 1:sizeBW(1);
    [n{:}] = ndgrid(n{:});

    % Interpolate traction on the deformed grid
    for i = 1:3
        F = scatteredInterpolant(x(:,1), x(:,2), ti(:,i), 'natural', 'none');
       Ti{i} = F(n{1},n{2});
    end


    % Density of cone plot
    plotdensity = 1/density^2;

    % Density_mask;
    cone_mask = rand(size(BW_mip{1}));
    cone_mask(cone_mask>plotdensity) = 0;
    cone_mask(cone_mask>0) = 1;

    % Finding index of cone position;
    idx = find(cone_mask);

    [vert{2},vert{1},vert{3}] = ind2sub(size(cone_mask),idx);

    % Calibration for x,y,z directions
%     v{1}{1} = ti(idx,2);
%     v{1}{2} = ti(idx,1);
%     v{1}{3} = ti(idx,3);
    v{1}{1} = Ti{2}(idx);
    v{1}{2} = Ti{1}(idx);
    v{1}{3} = Ti{3}(idx);
    mag = sqrt(v{1}{1}.^2 + v{1}{2}.^2 + v{1}{3}.^2);
    idx = vert{1} > xmin & vert{1} < xmax & vert{1} > ymin & vert{1} < ymax ;
    timag_max = max(mag(idx));
%     timag_max = max(mag);


    % Plot
    figure;
    tt = zeros(128,3);
    ttt = [179,136,255];
    for i = 1:128
        tt(i,:) = [0,i/128,0];
        tt(i,:) = ttt/255*i/128/3+ttt/255/3*2;
    end
    tt(1,:) = [0,0,0];
    for i = 2:40
        tt(i,:) = tt(41,:)*i/41;
    end

    colormap([tt;turbo(128)]);
    IM = timag_max.*double(im{timePt})-timag_max;

    [~,hC]= contourf((IM),128);
    set(hC,'LineStyle','none');
    axis image
    hold on

    hc = coneplot(vert{1},vert{2},vert{3},v{1}{1},v{1}{2},v{1}{3},coneSize,'nointerp');
%     hc = coneplot(x(:,1),x(:,2),x(:,3),ti(:,1),ti(:,2),ti(:,3),0.03,'nointerp');
    caxis([-timag_max,timag_max])
    fvc = repmat(mag(:)',[42 1]);
%     fvc = repmat(ti(:,4)',[42 1]);
    set(hc, 'facecolor', 'flat', 'facevertexcdata', fvc(:));
    hc.EdgeColor = 'none';
    hc.AmbientStrength = 0.6;
    hc.DiffuseStrength = 0.75;
    hc.SpecularStrength = 0.4;
    hc.SpecularExponent = 3;
    h_color = colorbar;
    set(h_color, 'ylim', [0, timag_max]);
    axis off;
    lighting phong;

    xlim([xmin, xmax])
    ylim([ymin, ymax])
    ylabel(h_color, 'Traction Magnitude (pN)');

    title({'Cell Traction'; ['Multipoint: ' mpname, ', Timepoint: ' num2str(timePt)]})
    savefigname = [savedir, mpname, '_trac_tp_' num2str(timePt), '.png'];
    % Optional scale bar
%     scalelength = 25; % desired scale bar length in um
%      quiver(xmin + 15,ymin + 15,scalelength/um2px(1),0,'ShowArrowHead','off', 'LineWidth', 2, 'Color', [1 1 1])

    % Save figure
    saveas(hc, savefigname);

end


end
