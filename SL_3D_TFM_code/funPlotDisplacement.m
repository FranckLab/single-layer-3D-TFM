function [cellcentroid] = funPlotDisplacement(u_plane, gridPts, BW, um2px, density, coneSize, mpname, isfigCrop, figCrop, savefolder)
% Plots and saves figures of a cone plot of the planarized displacement
% data using colormap Turbo (Copyright 2020, Daniel Fortunato.
% See: https://www.mathworks.com/matlabcentral/fileexchange/74662-turbo) along
% with the cell outline. Uses a previously user-generated and saved binarized
% cell outline in 'BW'.
%
%
%--- INPUTS ---
% u_plane   : planarized displacement data
% gridPts   : gridded plane on which displacement data is computed
% BW        : 3D matrix of black and white (double format) cell image
% um2px     : um to pixel ratio of images
% density   : density of cones for coneplot
% coneSize 	: size of cones for coneplot
% mpname    : name of multipoint data
% isfigCrop : 0 or 1 to crop or not crop the data
% figCrop   : [xmin xmax; ymin ymax] format borders for figure cropping
% saveFolder: output folder to which to save figures
%
%
%--- OUTPUTS ---
% cellcentroid  : location of cell centroid (from BW) for FEA
%
% NOTES
% ----------------------------------------------------------------------
% May, 2020; Lauren Hazlett, Alex Landauer
% Franck Lab, Brown Univerisity and University of Wisc - Madison
%
% If used please cite:
%


%% ~~~~~~~~~~~ Plot displacement data ~~~~~~~~~~~~

%% Format Displacements and produce cone plot
% folder = ['.',filesep,'data',filesep];
% [~,name,~] = fileparts(folder{1});
savedir = [savefolder,'Figures',filesep, mpname, filesep]; %save in subdir "PSF"
if exist(savedir,'dir') ~= 7 %make a new output folder if none exists
    mkdir(savedir);
end

for timePt = 1 : length(u_plane)
    fprintf('\nPlotting timepoint %i of %i\n',timePt,length(u_plane))
    
    % Initialize cell border information
    im{timePt} = squeeze(sum(BW{timePt},3));
    im{timePt} = im{timePt}/max(im{timePt}(:));
    sizeBW = size(BW{timePt});
    
    if isfigCrop == 1
        xmin = figCrop(2);
        ymin = figCrop(1);
        xmax = figCrop(4);
        ymax = figCrop(3);
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
    
    % Initializie displacement information
    xygrid = gridPts{timePt};
    xygrid{3} = ones(size(xygrid{1}));
    U = u_plane{timePt};
    U{4} = sqrt(U{1}.^2 + U{2}.^2 + U{3}.^2);
    
    % Density of cone plot
    plotdensity = 1/density^2;
    
    % Density_mask;
    cone_mask = rand(size(U{1}));
    cone_mask(cone_mask>plotdensity) = 0;
    cone_mask(cone_mask>0) = 1;
    
    % Finding index of cone position;
    idx = find(cone_mask);
    
    [vert{2},vert{1},vert{3}] = ind2sub(size(cone_mask),idx);
    
    % Calibration for x,y,z directions
    v{1}{1} = U{2}(idx);
    v{1}{2} = U{1}(idx);
    v{1}{3} = U{3}(idx);
    mag = sqrt(v{1}{1}.^2 + v{1}{2}.^2 + v{1}{3}.^2);
    idx = vert{1} > xmin & vert{1} < xmax & vert{1} > ymin & vert{1} < ymax ;
    umag_max = max(mag(idx));
    
    % Plot
    figure
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
    IM = umag_max*double(im{timePt})-umag_max;
    
    [~,hC]= contourf((IM),128);
    set(hC,'LineStyle','none');
    axis image
    hold on
    
    hc = coneplot(vert{1},vert{2},vert{3},v{1}{1},v{1}{2},v{1}{3},coneSize,'nointerp');
    caxis([-umag_max,umag_max])
    fvc = repmat(mag(:)',[42 1]);
    set(hc, 'facecolor', 'flat', 'facevertexcdata', fvc(:));
    hc.EdgeColor = 'none';
    hc.AmbientStrength = 0.6;
    hc.DiffuseStrength = 0.75;
    hc.SpecularStrength = 0.4;
    hc.SpecularExponent = 3;
    h_color = colorbar;
    set(h_color, 'ylim', [0, umag_max]);
    axis off;
    lighting phong;
    
    xlim([xmin, xmax])
    ylim([ymin, ymax])
    ylabel(h_color, 'Displacement (um)');
    
    title({'Cell Displacement'; ['Multipoint: ' mpname, ', Timepoint: ' num2str(timePt)]},'interpreter','none')
    savefigname = [savedir, mpname, '_disp_tp_' num2str(timePt), '.png'];
    % Optional scale bar
    %     scalelength = 25; % desired scale bar length in um
    %     quiver(xmin + 15,ymin + 15,scalelength/um2px(1),0,'ShowArrowHead','off', 'LineWidth', 2, 'Color', [1 1 1])
    
    % Save figure
    saveas(hc, savefigname);
    drawnow
    
    repick = input('Pick a location of max displacement (Y/N)? [Default: cell centroid] ','s');
    
    if strcmp(repick,'Y') || strcmp(repick,'y')
        figure
        
        colormap([tt;turbo(128)]);
        [~,hC]= contourf((IM),128);
        set(hC,'LineStyle','none');
        axis image
        hold on
        
        hc = coneplot(vert{1},vert{2},vert{3},v{1}{1},v{1}{2},v{1}{3},coneSize,'nointerp');
        caxis([-umag_max,umag_max])
        fvc = repmat(mag(:)',[42 1]);
        set(hc, 'facecolor', 'flat', 'facevertexcdata', fvc(:));
        hc.EdgeColor = 'none';
        hc.AmbientStrength = 0.6;
        hc.DiffuseStrength = 0.75;
        hc.SpecularStrength = 0.4;
        hc.SpecularExponent = 3;
        h_color = colorbar;
        set(h_color, 'ylim', [0, umag_max]);
        axis off;
        lighting phong;
        
        xlim([xmin, xmax])
        ylim([ymin, ymax])
        ylabel(h_color, 'Displacement (um)');
        
        %now add centroid marker
        Sc = scatter3(sizeBW(1)*cellcentroid{timePt}(1),sizeBW(2)*cellcentroid{timePt}(2),15);
        Sc.SizeData = 200;
        Sc.Marker = 'x';
        Sc.MarkerEdgeColor = '#00FF00';
        Sc.LineWidth = 3;
        text(0,-25,' Green ''x'' marks computed cell centroid','FontSize',16)
        drawnow
        
        title('Pick point of max displacement for mesh refinement: enter x- and y-coords','interpreter','none')
        xpoint = input('Enter x-coordinate: ');
        ypoint = input('Enter y-coordinate: ');
        
        cellcentroid{timePt} = [ypoint/sizeBW(1) xpoint/sizeBW(2)];
        
    end
    
end
end
