function [ti, tiPN] = funCalculateTractions(varargin)
% [ti, tiPN] = calculateTractions(Sij,Sij_pts,s,n); is used to calculate
% surface traction 3D vectors for a given stress state and deformed surface
% configuration
%
% INPUT: Sij,Sij_pts, s, n
% ------------------------------------------------------------------------
% Sij: Cauchy Stress tensor
% Sij_pts: vertices where Caushy stress is computed
% s: surface coordinates where normals are defined
% n: surface normals
%
% OUTPUT: ti, tiPN
% -----------------------------------------------------------------------
% ti: surface traction 3D vector components per each time stack.
%     Output format - ti{timePoint}{component}. ti{t}{4} - magnitude
% tiPN: in-plane and out-of-plane components of surface traction vector per
%       each time stack.
%       Output format: tiPN{t}{1}: In-plane or tiPN{t}{2}: Out-of-plane}
%
% NOTES
% ----------------------------------------------------------------------
% Adapted from Toyjanova et al. 2014 (DOI:10.1371/journal.pone.0090976),
% see: https://github.com/FranckLab/LD-3D-TFM
%
% For more information please see Section "Estimating 3D LD Cellular
% Tractions." in the aforementioned paper
%
% If used please cite:
%

[Sij,Sij_pts,s,n] = parseInputsTi(varargin{:});

interpMethod = 'natural';
maxTime = length(Sij);

ti = cell(maxTime,1);
for i = 1:maxTime
    ti{i} = funCalculateTractionsSub(Sij{i},Sij_pts{i},s{i},n{i},interpMethod);
end


tiPN = cell(maxTime,1);
for k = 1:maxTime
    tiNormal{k} = abs(ti{k}{1}.*n{k,1}{1}+ti{k}{2}.*n{k,1}{2}+ti{k}{3}.*n{k,1}{3});
    for m = 1:3
        tiParallel_{k}{m} = ti{k}{m}-tiNormal{k}.*n{k,1}{m};
    end
    tiParallel{k,1} = sqrt(tiParallel_{k}{1}.^2+tiParallel_{k}{2}.^2+tiParallel_{k}{3}.^2);
    tiPN{k}{1} = tiParallel{k};
    tiPN{k}{2} = tiNormal{k};
end


end

%% ========================================================================
function ti = funCalculateTractionsSub(Sij,Spts,s,n,interpMethod)

s33 = sampleData3(Sij{3,3},Spts,s,interpMethod);
s11 = sampleData3(Sij{1,1},Spts,s,interpMethod);
s22 = sampleData3(Sij{2,2},Spts,s,interpMethod);

s13 = sampleData3(Sij{1,3},Spts,s,interpMethod);
s23 = sampleData3(Sij{2,3},Spts,s,interpMethod);
s12 = sampleData3(Sij{1,2},Spts,s,interpMethod);

% Time to calculate those tractions
ti{1} = s11.*n{1} + s12.*n{2} + s13.*n{3};
ti{2} = s12.*n{1} + s22.*n{2} + s23.*n{3};
ti{3} = s13.*n{1} + s23.*n{2} + s33.*n{3};

for i = 1:3
    ti{i} = inpaint_nans(ti{i});
    ti{i} = real(ti{i});
end

ti{4} = sqrt(ti{1}.^2 + ti{2}.^2 + ti{3}.^2);


end

%% ========================================================================
function f1 = sampleData3(f0,pts,s,method)
%scattered vertex locations of stress and scatter tri-points of normals
sc_interp = scatteredInterpolant(pts{1}(:),pts{2}(:),pts{3}(:),f0',method);
%return the stress sampled at the tri-points
f1 = sc_interp(s{1},s{2},s{3});
end

%% ========================================================================
function varargout = parseInputsTi(varargin)
%  [Sij,m,s,n,OS] = parseInputs(Sij,pts,s,n,OS)
Sij = varargin{1};
pts = varargin{2};
s = varargin{3};
n = varargin{4};

varargout{      1} = Sij;
varargout{end + 1} = pts;
varargout{end + 1} = s;
varargout{end + 1} = n;

end
