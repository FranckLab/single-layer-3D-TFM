function [match, iA, uMatch] = match121(ssmatch,idx0,idx1)
% Find neighbors of unique match in original image idx

if length(ssmatch)>0

    [u,ia,~] = unique(ssmatch);
    n = histc(ssmatch,u);
    n = n==1;
    ia = nonzeros(ia.*n);
    iA = idx0(ia);

    match = ssmatch(ia);
else
    match = [];
    iA = [];
end

uMatch = 0;

end
