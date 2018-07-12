%% This code segments an image into n*n non-overlapping patches
%
% Input:    im - image
%           n - height and width of the patch
% Output:   lp - labelled image of patches
%           am - Adjeacency matrix of patches
%           pt - Patch attribute structure array with fields:
%                   L  - Mean L* value
%                   a  - Mean a* value
%                   b  - Mean b* value
%                   r  - Mean row value
%                   c  - Mean column value
%%
function [lp, am, pt] = patchUniformSampling(im, n)
if nargin<2
    n = 10; % Default height and width of the patch
end
% Constructing labelled non-overlapping patches
rows = round(size(im,1)/n);
cols = round(size(im,2)/n);
no_of_pt = rows*cols;
im = imresize(im,[rows*n,cols*n]);   % Resizing an image into dimension that is divisible by n

lp = zeros([size(im,1), size(im,2)]);   % Initializing labelled patches
pt = struct('L', cell(1,no_of_pt), 'a', cell(1,no_of_pt), 'b',  cell(1,no_of_pt), 'r', cell(1,no_of_pt), 'c', cell(1,no_of_pt));

k = 1;  % Label of the patches
for i = 1 : rows
    for j = 1 : cols
        lp(((i*n)-(n-1)):(i*n),((j*n)-(n-1)):(j*n)) = k;
        pt(k).r = (((i*n)-(n-1))+(i*n))/2;
        pt(k).c = (((j*n)-(n-1))+(j*n))/2;
        k = k+1;
    end
end

% Constructing adjacency matrix using labelled patches
[am] = regionadjacency(lp, 8);  

% Abstracting each patch to calculate the average Lab and position xy
lab_image = rgb2lab(im);    % Corresponding Lab image of the input image
L_channel = lab_image(:,:,1);
a_channel = lab_image(:,:,2);
b_channel = lab_image(:,:,3);
for i = 1 : no_of_pt
    pt(i).L = mean(mean(L_channel(lp==i)));
    pt(i).a = mean(mean(a_channel(lp==i)));
    pt(i).b = mean(mean(b_channel(lp==i)));
end
end
