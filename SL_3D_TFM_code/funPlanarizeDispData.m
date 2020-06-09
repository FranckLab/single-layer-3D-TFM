function [zPlaneCoeff0, zPlaneCoeff1, x0tracked, x1tracked, u_scplane, u_plane, gridPts]...
    = funPlanarizeDispData(x0, x1, x, track, u, mpname, runMode)
% Finds best fit plane for single layer bead images by removing out-of-plane
% particles from the image and fitting a bicubic function to the planarized
% particle data. Regularizes displacement data onto the grid using
% regularizeNd code from Matlab file exchange (Jason Nicholson (2020). regularizeNd
% (https://www.mathworks.com/matlabcentral/fileexchange/61436-regularizend),
% MATLAB Central File Exchange.)
%
%
%--- INPUTS ---
% x0      : particle positions in reference image
% x1      : particle positions in displaced image
% x       : particle locations before and after displacement
% track   : the output tracking results from TPT
% u       : particle displacements
% mpname  : name of multipoint data
% runMode : incremental or cumulative time steps in TPT
%
%--- OUTPUTS ---
% zPlaneCoeff0  : bicubic fitting function coefficients for reference z-plane
% zPlaneCoeff1  : bicubic fitting function coefficients for deformed z-plane
% x0tracked     : tracked particles in plane in reference image
% x1tracked     : tracked particles in plane in deformed image
% u_scplane     : scattered displacement data with outliers from plane removed
% u_plane       : planarized displacement data
% gridPts       : gridded plane on which displacement data is computed
%
% NOTES
% ----------------------------------------------------------------------
% March, 2020; Lauren Hazlett, Jin Yang, Alex Landauer
% Franck Lab, Brown Univerisity and University of Wisc - Madison
%
% If used please cite:
%


% ~~~~~~~~~~~ Planarization of scattered TPT displacment data ~~~~~~~~~~~~
%% Initialize for plots

disp('Removing plane outliers')

% ====== Fit gel interface shape z=z(x,y) for FEA boundary conditions ======
% FUTURE: zero strain warping of initial FEA mesh to match top surface of gel
for t = 1:length(track)
    xList{t} = linspace(floor(min(x0{t}(:,1))),ceil(max(x0{t}(:,1))),1e2);
    yList{t} = linspace(floor(min(x0{t}(:,2))),ceil(max(x0{t}(:,2))),1e3);
    [gridPts{t}{2},gridPts{t}{1}] = meshgrid(yList{t},xList{t});

end

% Select beads near interface: median +/- std
% Select beads near interface for ref image
zPlane_std_alpha = 3; %starting 'default' value of stds to use
change_param = 0;
while change_param == 0
    for t = 1:length(track)
        %get the top and bottom extents
        zPlaneFitMin{t} = min(x0{t}(:,3));
        zPlaneFitMax{t} = max(x0{t}(:,3));
        for tempi = 1:3
            [x0tracked{t}] = find(x0{t}(:,3)>zPlaneFitMin{t} & x0{t}(:,3)<zPlaneFitMax{t} );
            indx0_3_ref_median = median(x0{t}(x0tracked{t},3));
            indx0_3_ref_std = std(x0{t}(x0tracked{t},3));
            zPlaneFitMin{t} = indx0_3_ref_median - zPlane_std_alpha*(4-tempi)*indx0_3_ref_std;
            zPlaneFitMax{t} = indx0_3_ref_median + zPlane_std_alpha*(4-tempi)*indx0_3_ref_std;
        end
        figure; plot3(x0{t}(:,1),x0{t}(:,2),x0{t}(:,3),'k.')
        title(['Tracked beads in reference image, multipoint: ',mpname,', timepoint: ', num2str(t)],'interpreter','none')
        hold on
        plot3(x0{t}(x0tracked{t},1),x0{t}(x0tracked{t},2),x0{t}(x0tracked{t},3),'ro')
    end
    prompt = '  Do you want to adjust outlier removal from the plane? Y/N?: ';
    YN = input(prompt, 's');
    if YN == 'Y' || YN == 'y'
        prompt = '    How many standard deviations away from the mean to include points in the plane (default: 3)?: ';
        zPlane_std_alpha = input(prompt);
    elseif YN == 'N' || YN == 'n'
        change_param = 1;
    end
    close all;
end

% Select beads near interface for deformed image using the same cutoff
% parameters, but the deformed positions
for t = 1:length(track)
    zPlaneFitMin1{t} = min(x1{t}(:,3));
    zPlaneFitMax1{t} = max(x1{t}(:,3));
    for tempi = 1:3
        [x1tracked{t}] = find(x1{t}(:,3)>zPlaneFitMin1{t} & x1{t}(:,3)<zPlaneFitMax1{t} );
        indx1_3_def_median = median(x1{t}(x1tracked{t},3)); indx1_3_def_std = std(x1{t}(x1tracked{t},3));
        zPlaneFitMin1{t} = indx1_3_def_median - zPlane_std_alpha*(4-tempi)*indx1_3_def_std;
        zPlaneFitMax1{t} = indx1_3_def_median + zPlane_std_alpha*(4-tempi)*indx1_3_def_std;
    end
end
disp('Plane outliers removed')

% ====== Fit interface shape using bicubic function ======
% Interface height Bicubic function has the formula
% z(x,y) = a0 + a1*x + a2*y + a3*x^2 + a4*x*y + a5*y^2 + a6*x^3 + a7*x^2*y + a8*x*y^2 + a9*y^3;

cMap = hsv(255);
disp('Generating plane fits')
for t = 1:length(x0)
    figure; subplot(1,2,1); plot3(x0{t}(:,1),x0{t}(:,2),x0{t}(:,3),'b.');
    hold on
    plot3(x0{t}(x0tracked{t},1),x0{t}(x0tracked{t},2),x0{t}(x0tracked{t},3),'ko');
    title(['Plane fit in reference image, multipoint: ',mpname,...
        ', timepoint: ', num2str(t)],'fontweight','normal','interpreter','none');

    % Cubic interface shape model
    Atemp0 = [ones(length(x0tracked{t}),1), x0{t}(x0tracked{t},1), x0{t}(x0tracked{t},2), ...
        x0{t}(x0tracked{t},1).^2, x0{t}(x0tracked{t},1).*x0{t}(x0tracked{t},2), x0{t}(x0tracked{t},2).^2 , ...
        x0{t}(x0tracked{t},1).^3, x0{t}(x0tracked{t},1).^2.*x0{t}(x0tracked{t},2), ...
        x0{t}(x0tracked{t},1).*x0{t}(x0tracked{t},2).^2, x0{t}(x0tracked{t},2).^3];

    % solve for the cubic surface coeffs for the reference image
    zPlaneCoeff0{t} = Atemp0\x0{t}(x0tracked{t},3);
    zGrid0{t} = zPlaneCoeff0{t}(1) + zPlaneCoeff0{t}(2)*gridPts{t}{1} + zPlaneCoeff0{t}(3)*gridPts{t}{2} + ...
        zPlaneCoeff0{t}(4)*gridPts{t}{1}.^2 + zPlaneCoeff0{t}(5)...
        *gridPts{t}{1}.*gridPts{t}{2} + zPlaneCoeff0{t}(6)*gridPts{t}{2}.^2 + ...
        zPlaneCoeff0{t}(7)*gridPts{t}{1}.^3 + zPlaneCoeff0{t}(8)*gridPts{t}{1}.^2.*gridPts{t}{2} + ...
        zPlaneCoeff0{t}(9)*gridPts{t}{1}.*gridPts{t}{2}.^2 + zPlaneCoeff0{t}(10)*gridPts{t}{2}.^3;

    % plot for user evaluation
    surf(gridPts{t}{1},gridPts{t}{2},zGrid0{t},'edgecolor','none' );
    cb=colorbar; colormap(cMap);
    cb.Label.String = 'Height of plane (pixels)';

    subplot(1,2,2); plot3(x1{t}(:,1),x1{t}(:,2),x1{t}(:,3),'r.');
    hold on
    plot3(x1{t}(x1tracked{t},1),x1{t}(x1tracked{t},2),x1{t}(x1tracked{t},3),'ko');
    title(['Plane fit in deformed image, multipoint: ',mpname,', timepoint: ',...
        num2str(t)],'fontweight','normal','interpreter','none');

    % Cubic interface shape model
    Atemp1 = [ones(length(x1tracked{t}),1), x1{t}(x1tracked{t},1), x1{t}(x1tracked{t},2), ...
        x1{t}(x1tracked{t},1).^2, x1{t}(x1tracked{t},1).*x1{t}(x1tracked{t},2), x1{t}(x1tracked{t},2).^2 , ...
        x1{t}(x1tracked{t},1).^3, x1{t}(x1tracked{t},1).^2.*x1{t}(x1tracked{t},2), ...
        x1{t}(x1tracked{t},1).*x1{t}(x1tracked{t},2).^2, x1{t}(x1tracked{t},2).^3];

    % solve for deformed plane fit coeffs
    zPlaneCoeff1{t} = Atemp1\x1{t}(x1tracked{t},3);
    zGrid1{t} = zPlaneCoeff1{t}(1) + zPlaneCoeff1{t}(2)*gridPts{t}{1} + zPlaneCoeff1{t}(3)*gridPts{t}{2} + ...
        zPlaneCoeff1{t}(4)*gridPts{t}{1}.^2 + zPlaneCoeff1{t}(5)...
        *gridPts{t}{1}.*gridPts{t}{2} + zPlaneCoeff1{t}(6)*gridPts{t}{2}.^2 + ...
        zPlaneCoeff1{t}(7)*gridPts{t}{1}.^3 + zPlaneCoeff1{t}(8)*gridPts{t}{1}.^2.*gridPts{t}{2} + ...
        zPlaneCoeff1{t}(9)*gridPts{t}{1}.*gridPts{t}{2}.^2 + zPlaneCoeff1{t}(10)*gridPts{t}{2}.^3;

    % plot deformed surface for evaluation
    surf(gridPts{t}{1},gridPts{t}{2},zGrid1{t},'edgecolor','none');
    cb=colorbar; colormap(cMap);
    cb.Label.String = 'Height of plane (pixels)';
end
disp('Plane fitting complete')

%% Interpolate disp to regular mesh and apply regularization
% figure, plot3(x0(x0tracked,1),x0(x0tracked,2),x0(x0tracked,3),'r.');
disp('Regularizing Displacement Data')

for t = 1:length(track)
    u_scplane{t} = u{t}(x0tracked{t},:);

    winstepsize = [1,1,1];
    xList{t} = ceil(min(x0{t}(x0tracked{t},1)))-winstepsize(1):winstepsize(1):...
        floor(max(x0{t}(x0tracked{t},1)))+winstepsize(1);
    yList{t} = ceil(min(x0{t}(x0tracked{t},2)))-winstepsize(2):winstepsize(2)...
        :floor(max(x0{t}(x0tracked{t},2)))+winstepsize(2);
    [gridPts{t}{2},gridPts{t}{1}] = meshgrid(yList{t},xList{t});
end
% ====== Compute disp x and y with regularization ======
smoothness_xy = 1e-5;
smoothnessChangeOrNot = 0;
while smoothnessChangeOrNot == 0
    for t = 1:length(track)
        % ------------- Regularization ----------------
        u_plane{t}{1} = regularizeNd([x0{t}(x0tracked{t},1),x0{t}(x0tracked{t},2)],...
            u_scplane{t}(:,1), {xList{t},yList{t}}, smoothness_xy);
        u_plane{t}{2} = regularizeNd([x0{t}(x0tracked{t},1),x0{t}(x0tracked{t},2)],...
            u_scplane{t}(:,2), {xList{t},yList{t}}, smoothness_xy);

        % ------------- Visualization -----------------
        figure; subplot(1,2,1); surf(gridPts{t}{1},gridPts{t}{2},u_plane{t}{1},'edgecolor','none');
        hold on
        plot3(x0{t}(x0tracked{t},1),x0{t}(x0tracked{t},2),u_scplane{t}(:,1),'ro');
        set(gcf,'color','w');
        set(gca,'fontsize',18);
        colormap(cMap);
        box on; grid on; grid minor;
        title({'x-Displacement (um) Smoothness Check'; ['Multipoint: ',mpname,...
            ', Timepoint: ', num2str(t)]},'fontweight','normal','interpreter','none');
        ax = gca;
        ax.Position = ax.Position - [0 0 0.2 0.2];
        cb = colorbar('Location','eastoutside');
        set(cb,'fontsize',16);
        cb.Position = cb.Position + [0.24 0 0 0];

        % ------------- Visualization -----------------
        subplot(1,2,2); surf(gridPts{t}{1},gridPts{t}{2},u_plane{t}{2},'edgecolor','none');
        hold on
        plot3(x0{t}(x0tracked{t},1),x0{t}(x0tracked{t},2),u_scplane{t}(:,2),'ro')
        set(gca,'fontsize',18);  colormap(cMap); box on; grid on; grid minor;
        title({'y-Displacement (um) Smoothness Check'; ['Multipoint: ',mpname,...
            ', Timepoint: ', num2str(t)]},'fontweight','normal','interpreter','none');
        ax = gca; ax.Position = ax.Position - [0 0 .2 .2];
        cb = colorbar('Location','eastoutside');  set(cb,'fontsize',16);
        cb.Position = cb.Position + [0.24 0 0 0]; %caxis([-1.2,0.6]);

    end

    prompt = '  Do you want to change xy regularization smoothness? Y/N?: ';
    YN_xysmooth = input(prompt, 's');
    if YN_xysmooth == 'Y' || YN_xysmooth == 'y'
        prompt = '    Input new xy smoothness parameter (Default: 1e-5): ';
        smoothness_xy = input(prompt);
    elseif YN_xysmooth == 'N' || YN_xysmooth == 'n'
        smoothnessChangeOrNot = 1;
        disp(' xy Regularization complete')
    end
    close all;
end

% ====== Compute disp z with regularization ======
smoothness_z = 1e-5;
smoothnessChangeOrNot = 0;
while smoothnessChangeOrNot == 0
    for t = 1:length(track)
        % ------------- Regularization ----------------
        u_plane{t}{3} = regularizeNd([x0{t}(x0tracked{t},1),x0{t}(x0tracked{t},2)],...
            u_scplane{t}(:,3), {xList{t},yList{t}}, smoothness_z);

        % ------------- Visualization -----------------
        figure; surf(gridPts{t}{1},gridPts{t}{2},u_plane{t}{3},'edgecolor','none');
        hold on
        plot3(x0{t}(x0tracked{t},1),x0{t}(x0tracked{t},2),u_scplane{t}(:,3),'ro')
        set(gca,'fontsize',18);  colormap(cMap); box on; grid on; grid minor;
        title({'z-Displacement (um) Smoothness Check'; ['Multipoint: ',mpname,...
            ', Timepoint: ', num2str(t)]},'fontweight','normal','interpreter','none');
        ax = gca; ax.Position = ax.Position - [0 0 .2 .2];
        cb = colorbar('Location','eastoutside');  set(cb,'fontsize',16);
        cb.Position = cb.Position + [0.24 0 0 0]; % caxis([-1.2,0.6]);
    end

    prompt = '  Do you want to change z regularization smoothness? Y/N?: ';
    YN_zsmooth = input(prompt, 's');
    if YN_zsmooth == 'Y' || YN_zsmooth == 'y'
        prompt = '    Input new z smoothness parameter (Default: 1e-5): ';
        smoothness_x = input(prompt);
    elseif YN_zsmooth == 'N' || YN_zsmooth == 'n'
        smoothnessChangeOrNot = 1;
        disp(' z Regularization complete')
    end
    close;
end

% % ====== Finally output all the figures ======
close all;
disp('Regularization Complete')
for t = 1 : length(x0)

    figure, surf(gridPts{t}{1},gridPts{t}{2},u_plane{t}{1},'edgecolor','none');
    %hold on, plot3(x0(x0tracked,1),x0(x0tracked,2),dispx,'ro'); set(gcf,'color','w');
    set(gca,'fontsize',18);  colormap(cMap); box on; grid on; grid minor;
    title({'x-Displacement (um)'; ['Multipoint: ',mpname,', Timepoint: ',...
        num2str(t)]},'fontweight','normal','interpreter','none');
    ax = gca; ax.Position = ax.Position - [0 0 .2 .2];
    cb = colorbar('Location','eastoutside');  set(cb,'fontsize',16);
    cb.Position = cb.Position + [0.24 0 0 0]; %caxis([-1.2,0.6]);

    figure, surf(gridPts{t}{1},gridPts{t}{2},u_plane{t}{2},'edgecolor','none');
    %hold on, plot3(x0(x0tracked,1),x0(x0tracked,2),dispy,'ro'); set(gcf,'color','w');
    set(gca,'fontsize',18);  colormap(cMap); box on; grid on; grid minor;
    title({'y-Displacement (um)'; ['Multipoint: ',mpname,', Timepoint: ',...
        num2str(t)]},'fontweight','normal','interpreter','none');
    ax = gca; ax.Position = ax.Position - [0 0 .2 .2];
    cb = colorbar('Location','eastoutside');  set(cb,'fontsize',16);
    cb.Position = cb.Position + [0.24 0 0 0]; %caxis([-1.2,0.6]);

    figure, surf(gridPts{t}{1},gridPts{t}{2},u_plane{t}{3},'edgecolor','none');
    %hold on, plot3(x0(x0tracked,1),x0(x0tracked,2),dispz,'ro'); set(gcf,'color','w');
    set(gca,'fontsize',18);  colormap(cMap); box on; grid on; grid minor;
    title({'z-Displacement (um)'; ['Multipoint: ',mpname,', Timepoint: ',...
        num2str(t)]},'fontweight','normal','interpreter','none');
    ax = gca; ax.Position = ax.Position - [0 0 .2 .2];
    cb = colorbar('Location','eastoutside');  set(cb,'fontsize',16);
    cb.Position = cb.Position + [0.24 0 0 0]; % caxis([-1.2,0.6]);
end


end
