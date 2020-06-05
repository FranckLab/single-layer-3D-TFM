function [output_field] = spatial_rolloff_filt(scalar_field,filt_str,edge_wd)
%Function to modify scalar field data to rolloff at the edges of a domain,
%according to a gaussian decay of user-defined strenght and starting
%a user-defined distance from the edge. This function takes in
% - vect_field: a scalar field (e.g. image, vector components)
% - filt_str: the strength of the Gaussian to apply, the roll-off width
% - edge_wd: the width of the zero-domain (pre-filt) at the edge
%
% Alex Landauer, March 2020

field_sz = size(scalar_field);

mask = ones(field_sz - edge_wd);
mask = padarray(mask,[floor(edge_wd/2),floor(edge_wd/2)],0,'pre');
mask = padarray(mask,[ceil(edge_wd/2),ceil(edge_wd/2)],0,'post');

mask = imgaussfilt(mask,filt_str);

output_field = scalar_field.*mask;
