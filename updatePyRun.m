function [] = updatePyRun(data_name_in,data_name_out,multipoint_name,E,nu,...
                           load_steps,numSteps,thickness,x_centroid_norm,y_centroid_norm)
%Function to update the Python runscript for FEniCS with new info for file
%names and material properties
%
%--- INPUTS ---
%  data_name_in   : input data file for FEniCS (coords, disps, params)
%  data_name_out  : the *prefix* for the output .mat data from FEniCS
%  multipoint_name: name of the current multipoint image series
%  E              : Young's modulus of the gel substrate (default: 1 kPa)
%  nu             : Poisson's ratio of the gel (usually ~0.49, fully
%                   incompressible may cause issues) (default: 0.49)
%  load_steps     : number of load steps used (default: 5)
%                  (fewer => faster; greater => more robust convergence)
%  numSteps       : number of time steps at the current multipoint
%                  (fewer => faster; greater => more robust convergence)
%  thickness      : thickness (estimate) of the gel, 0 => "thick gel (20x of
%                   max z-disp) assumption (default: 0)
%  x_center_norm  : center location of the cell in normalized coords (default: 0.5)
%  y_center_norm  : center location of the cell in normalized coords (default: 0.5)
%
%--- OUTPUTS ---
%  - none: write the updated .py file each run needed
%
% June, 2019; Alex Landauer, updated May 2020
% Franck Lab, Brown Univerisity and University of Wisc - Madison

if nargin < 4
    E = 1000.0; %1 kPa
end
if nargin < 5
    nu = 0.49; %nearly incomp
end
if nargin < 6
    load_steps = 5;
end
if nargin < 7
    numSteps = 2;
end
if nargin < 8
    thickness = 0;
end
if nargin < 9
    x_centroid_norm = 0.5;
end
if nargin < 10
    y_centroid_norm = 0.5;
end

%get names only
[~,nameIn,~] = fileparts(data_name_in);
[~,nameOut,~] = fileparts(data_name_out);

%% Read the old .py text file
fid = fopen('sl_tfm_call.py','r');
ii = 1;
tline = fgetl(fid);
PYfile{ii} = tline;
while ischar(tline)
    ii = ii+1;
    tline = fgetl(fid);
    PYfile{ii} = tline;
end
fclose(fid);

%% Change cells
PYfile{23} = sprintf('for stepNum in range(%i):',numSteps);
PYfile{25} = sprintf('\tdata_name_in = "%s_%%03d" %% (stepNum)',nameIn);
PYfile{26} = sprintf('\tdata_name_out = "%s_%%03d" %% (stepNum)',nameOut);
PYfile{27} = sprintf('\tload_steps = %i',load_steps);
PYfile{28} = sprintf('\tE = %0.5f',E);
PYfile{29} = sprintf('\tnu = %0.5f',nu);
PYfile{30} = sprintf('\tthickness = %0.5f',thickness);
PYfile{31} = sprintf('\tx_center_norm = %0.5f',x_centroid_norm);
PYfile{32} = sprintf('\ty_center_norm = %0.5f',y_centroid_norm);

%% Write the new .py text file
fid = fopen(['sl_tfm_call_',multipoint_name,'.py'], 'w');
for ii = 1:numel(PYfile)
    if PYfile{ii+1} == -1
        fprintf(fid,'%s', PYfile{ii});
        break
    else
        fprintf(fid,'%s\n', PYfile{ii});
    end
end
