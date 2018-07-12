%% This code clusters patches using uniform region sampling
% Input : 
%       ls - labelled patches in image
%       sp - patch color and location
%       z - number of row and columns of the regions
% Output :
%       lc - labelled clusters in image
%       csp - clustered patches
%       sc - labels of patches
%%
function [lc, csp, sc] = regionUniformSampling(ls, sp, z)
% read color values of superpixels
no_of_sp = length(sp);    % Number of superpixels
no_of_reg = z^2;    % Number of uniform regions
super_pixel_position = zeros([no_of_sp,2]); % Creating a vector to represent position of each superpixel (x,y)
for i = 1 : no_of_sp
    super_pixel_position(i,1) = nearest(sp(1,i).r);    super_pixel_position(i,2) = nearest(sp(1,i).c);    % Assigning (x,y) position of each superpixel
end;
% Assigning boundaries of unform regions
region_position = zeros([no_of_reg,4]);
rows = size(ls,1)/z;
cols = size(ls,2)/z;
k = 1;
for i = 1 : z
    for j = 1 : z
        region_position(k,1) = nearest((i*rows)-(rows-1));
        region_position(k,2) = nearest(i*rows);
        region_position(k,3) = nearest((j*cols)-(cols-1));
        region_position(k,4) = nearest(j*cols);
        k = k+1;
    end
end
% assigning each patch into a region
sc = zeros([no_of_sp,1]);   % Clustered superpixels
for i = 1 : no_of_sp
    for j = 1 : no_of_reg
        if(super_pixel_position(i,1)>=region_position(j,1) && super_pixel_position(i,1)<=region_position(j,2) && ...
                super_pixel_position(i,2)>=region_position(j,3) && super_pixel_position(i,2)<=region_position(j,4))
            sc(i,1) = j ;
        end
    end
end

% assigning labels to the superpixel clusters
lc=ls;
for i = 1 : no_of_sp
lc(lc == i) = sc(i);
end
% save the superpixels in each cluster using cell array
csp = {};
for i = 1 : no_of_reg
    [temp] = find(sc==i);
    temp = transpose(temp); % inversion to make it row matrix
    csp(i) = {temp};
end
end
