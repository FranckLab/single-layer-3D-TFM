function [n, tri] = calculateNormals(s)
% [n, tri] = calculateNormals(s) is a function that calculates

% surface normals from the surface grid,s
% 
% INPUTS
% ------------------------------------------------------------------------
%   s: cell array containing the surface coordinate grid defined by a
%      rectangular meshgrid (can be non-uniform)
%      s{time}{1} = x-coordinates of nodes
%      s{time}{2} = y-coordinates of nodes
%      s{time}{3} = z-coordinates of nodes
%             
% OUTPUTS
% -------------------------------------------------------------------------
%   tri:  delaunay tiranguluzation of the surface.
%   n: cell array containing the normal vectors for each point on the 
%      surface
%       n{time}{1} = x-component of normal vector
%       n{time}{2} = y-component of normal vector
%       n{time}{3} = z-component of normal vector

% NOTES
% -------------------------------------------------------------------------
% none
% 
% If used please cite:
% Toyjanova J., Bar-Kochba E., López-Fagundo C., Reichner J., 
% Hoffman-Kim D, Franck, C. (2014) High Resolution, Large Deformation 3D
% Traction Force Microscopy. PLoS ONE 9(4): e90976. 
% doi:10.1371/journal.pone.0090976



if ~iscell(s)
    error('surface must have at least one time point within a cell array'); 
end

maxTime = length(s);
n = cell(maxTime,1);
tri = cell(maxTime,1);
for i = 1:maxTime, [n{i}, tri{i}] = funCalculateNormals(s{i}); end

end

function [n, tri] = funCalculateNormals(s)

% tri:  delaunay tiranguluzation of the surface. see TriRep() or
% triangulation in newer releases of MATLAB
% s{1}, s{2}: meshgrid that the tractions are defined on

sSize = size(s{1});
dt = delaunay(s{1}(:),s{2}(:));
tri = triangulation(dt,[s{1}(:) s{2}(:) s{3}(:)]);
n_ = vertexNormal(tri);

n = cell(1,3);
for i = 1:3, n{i} = reshape(n_(:,i),sSize); end


end