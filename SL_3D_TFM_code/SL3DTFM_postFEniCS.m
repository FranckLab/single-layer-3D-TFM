%~~~~~~~~~~~~~ Single-layer TPT-based Traction Force Microscopy ~~~~~~~~~~~~~~~~
%
% Post-processing function to read in FEniCS results, compute surface normals and
% tractions, and do basic plotting
%
% June, 2019; contributors:  A Landauer,  L Hazlett, M Patel, H Witt
% Franck Lab, Brown Univerisity and University of Wisc - Madison
%
%

clear all

% Set up input data-- copy/paste this block from SL3DTFM_runfile

%% --------------------------- INPUT SET UP ------------------------------------
data_name = 'example_data'; %data folder that contains subfolder 'tifs', with tif stacks stored in multipoint subfolders
save_file_descriptor = 'example_saved_data';  % change for EACH day/data set to prevent files from saving over each other

cur_dir = dir();
data_dir = [cur_dir(1).folder,filesep,data_name,filesep];

multipoint_names = {'XYpoint_name_1','XYpoint_name_2'}; %if only one multipoint, leave '' in the brackets
cellBWfilename = 'exampleBWfile.mat';



%% ---------------------------- LOAD DISP DATA -----------------------------------
%   Load the prior setup and displacement data from the SL3DTFM_runfile.m script

setup_data = load([data_dir,data_name,'_',save_file_descriptor,'_','all','_settings.mat']);
disp_data = load([data_dir,data_name,'_',save_file_descriptor,'_displacementplottingdata.mat']);


%% --------------------------- LOAD FENICS DATA ---------------------------------
% Load the FEniCS output data files for each multipoint and time point (saved in
% the base working directory).

disp('Loading FEniCS Data')
for multipoint = 1:length(setup_data.multipoint_names)
    for timepoint = 1:(setup_data.total_images{multipoint}-1)
        tp = sprintf('%03d',timepoint-1);
        pts_fenics{multipoint}{timepoint} = load([data_name,'_',...
            setup_data.multipoint_names{multipoint},'_fenicsOut_',tp,'_pts_fenics.mat']);
        u_fenics{multipoint}{timepoint} = load([data_name,'_',...
            setup_data.multipoint_names{multipoint},'_fenicsOut_',tp,'_u_fenics.mat']);
        stress_fenics{multipoint}{timepoint} = load([data_name,'_',...
            setup_data.multipoint_names{multipoint},'_fenicsOut_',tp,'_stress.mat']);
    end
end
disp('Done Loading Data')

%% ------------------ COMPUTE SURFACE NORMALS AND STRESSES ---------------------
%   Use a Delauney triangularization of the deformed top surface of the mesh
%   to compute surface normals (for traction computation) at each vertex point.
%   Then use the Cauchy Stress Theorem to compute the surface tractions from the
%   3D volumetric stresses at the top surface, saved from FEniCS.

for multipoint = 1:length(setup_data.multipoint_names)

    fprintf('\n Computing surface normals and tractions for multipoint %i of %i\n',...
        multipoint,length(setup_data.multipoint_names))

    for timepoint = 1:(setup_data.total_images{multipoint}-1)
        vertices = [pts_fenics{multipoint}{timepoint}.x;...
            pts_fenics{multipoint}{timepoint}.y;pts_fenics{multipoint}{timepoint}.z];
        displacements = [u_fenics{multipoint}{timepoint}.u1;...
            u_fenics{multipoint}{timepoint}.u2;u_fenics{multipoint}{timepoint}.u3];

        % Find the top surface (undeformed config)
        tol = 1e-9;
        max_z = max(pts_fenics{multipoint}{timepoint}.z) - tol;

        top_vert = vertices(3,:)>max_z;

        for ii = 1:3
            top_surface{multipoint}{timepoint}{ii} = vertices(ii,top_vert);
            top_surface_disp{multipoint}{timepoint}{ii} = displacements(ii,top_vert);
            top_surface_deformed{multipoint}{timepoint}{ii} = ...
                top_surface{multipoint}{timepoint}{ii} + top_surface_disp{multipoint}{timepoint}{ii};
        end
    end
    % Find surface normals on the deformed configuration
    [n_i{multipoint},tri{multipoint}] = calculateNormals(top_surface_deformed{multipoint});

    %              |{A11}, {A12}, {A13}|
    %  Aij{time} = |{A21}, {A22}, {A23}|
    %              |{A11}, {A12}, {A13}|
    % For example, to access the A12 for the 1st time increment you would
    % call Aij{1}{1,2}.
    for timepoint = 1:(setup_data.total_images{multipoint}-1)
        for ii = 1:3
            pts_normals{multipoint}{timepoint}{ii} = tri{multipoint}{timepoint}.Points(:,ii);
        end

        % extract components of r the stress tensor from FEniCS outputs
        Sij{multipoint}{timepoint}{1,1} =  stress_fenics{multipoint}{timepoint}.s11(top_vert);
        Sij{multipoint}{timepoint}{2,2} =  stress_fenics{multipoint}{timepoint}.s22(top_vert);
        Sij{multipoint}{timepoint}{3,3} =  stress_fenics{multipoint}{timepoint}.s33(top_vert);

        Sij{multipoint}{timepoint}{1,2} =  stress_fenics{multipoint}{timepoint}.s12(top_vert);
        Sij{multipoint}{timepoint}{2,1} =  stress_fenics{multipoint}{timepoint}.s21(top_vert);
        Sij{multipoint}{timepoint}{1,3} =  stress_fenics{multipoint}{timepoint}.s13(top_vert);
        Sij{multipoint}{timepoint}{3,1} =  stress_fenics{multipoint}{timepoint}.s31(top_vert);
        Sij{multipoint}{timepoint}{2,3} =  stress_fenics{multipoint}{timepoint}.s23(top_vert);
        Sij{multipoint}{timepoint}{3,2} =  stress_fenics{multipoint}{timepoint}.s32(top_vert);
    end

    % Calcuate surface tractions at tri-points (Cauchy Relation)
    [ti{multipoint}, tiPN{multipoint}] = funCalculateTractions(Sij{multipoint},...
        pts_normals{multipoint},top_surface_deformed{multipoint},n_i{multipoint});
end

save([data_name,'_',save_file_descriptor,'_fe_tractions.mat'],'ti','top_surface','top_surface_disp','pts_normals')

disp('Surface normal and traction computations complete')


%% ----------------------------- PLOTTING TRACTION -----------------------------
%   Plot the tractions, displacements, and cell outline to visualize the results.
%   Several type of plots are useful, this is an example.  In addition to those
%   found here, we also like simple filled contour plots of traction (magnitude
%   or components) with a overlaid cell outline, using the Matlab contourf function

%% Traction and displacement visualization (simple quiver plot)
for multipoint = 1 : length(setup_data.multipoint_names)
    fprintf('\n Showing displacement and traction for multipoint %i of %i\n',...
        multipoint,length(setup_data.multipoint_names))
    for timepoint = 1:(setup_data.total_images{multipoint}-1)
    figure;
    quiver3(pts_normals{multipoint}{timepoint}{1}',pts_normals{multipoint}{timepoint}{2}',...
        pts_normals{multipoint}{timepoint}{3}',...
        top_surface_disp{multipoint}{timepoint}{ii},top_surface_disp{multipoint}{timepoint}{ii},...
        top_surface_disp{multipoint}{timepoint}{ii},1)
    title(['Surface Displacement Vectors, Multipoint: ', setup_data.multipoint_names{multipoint},...
        ', Time point: ' num2str(timepoint)],'interpreter','none')
    xlabel('x, \mum')
    ylabel('y, \mum')
    zlabel('z, \mum')
    
    figure
    trisurf(tri{multipoint}{timepoint},'FaceColor',[0.8 0.8 1.0]);
    axis equal

    hold on
    quiver3(pts_normals{multipoint}{timepoint}{1}',pts_normals{multipoint}{timepoint}{2}',...
        pts_normals{multipoint}{timepoint}{3}', ...
        ti{multipoint}{timepoint}{1},ti{multipoint}{timepoint}{2},ti{multipoint}{timepoint}{3},1,'Color','b');
    title(['Traction Vectors on FE Mesh, Multipoint: ', setup_data.multipoint_names{multipoint}, ...
        ', Time point: ' num2str(timepoint)],'interpreter','none')
    xlabel('x, \mum')
    ylabel('y, \mum')
    zlabel('z, \mum')
    
%     plot_trac_simple(tri{multipoint}{timepoint},ti{multipoint}{timepoint},top_surface_deformed{multipoint}{timepoint})
%     title({'Tractions shown as a cone plot, displacement as a surface contour plot';...
%     ['Multipoint: ', setup_data.multipoint_names{multipoint}, ', Time point: ' num2str(timepoint)]})
    end
end

%% Final traction visualization (cone plot + cell outline)
for multipoint = 1 : length(setup_data.multipoint_names)
    fprintf('\n Plotting traction for multipoint %i of %i\n',multipoint,length(setup_data.multipoint_names))

    % For visualization: density and size of cones for displacement coneplot
    density = 2;            % smaller number = higher density
    coneSize = 0.008;       % smaller number = smaller cone size

    funPlotTraction(ti{multipoint},top_surface_deformed{multipoint}, setup_data.um2px,...
        disp_data.BW{multipoint}, density, coneSize, setup_data.multipoint_names{multipoint},...
        setup_data.isfigCrop{multipoint}, setup_data.figCrop{multipoint}, data_dir)

    fprintf('\n Plotting traction complete \n')

end

save([data_dir,data_name,'_',save_file_descriptor,'_postFEniCssettings.mat'], '-v7.3');

% ------------------------------ END OF FILE -----------------------------------