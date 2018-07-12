%% This code segments an image into n superpixels using SLIC superpixel segmentation
%
% Input:    im - image
%           n - number of desired superpixels
% Output:   lp - labelled image of superpixels
%           am - Adjeacency matrix of superpixels
%           sp - superpixels attribute structure array with fields:
%                   L  - Mean L* value
%                   a  - Mean a* value
%                   b  - Mean b* value
%                   r  - Mean row value
%                   c  - Mean column value
%%
function [ls, am, sp] = patchSuperpixel(im, n)
if nargin<2
    n = 500; % Default height and width of the patch
end
addpath( './Slic/' );
m = 20; % Regularity parameter
Img = double(im);  % Converting an image into double
height = size(Img,1);
width = size(Img,2);
no_of_px = height*width;
R_channel = reshape( Img(:,:,1)', no_of_px, 1);
G_channel = reshape( Img(:,:,2)', no_of_px, 1);
B_channel = reshape( Img(:,:,3)', no_of_px, 1);
Img_attr=[ height ,width, n, m, no_of_px ];    % Image attribute

[ Label_line, Sup1, Sup2, Sup3, n ] = SLIC( R_channel, G_channel, B_channel, Img_attr );    % Calling SLIC superpixels algorithm
label=reshape(Label_line,width,height);
label = label';     % the superpixels' labels ranging from 0 to n-1
ls = label+1;    % Changing the label range to 1 to n and storing it to 'ls'

sp = struct('L', cell(1,n), 'a', cell(1,n), 'b',  cell(1,n), 'r', cell(1,n), 'c', cell(1,n));

% Constructing adjacency matrix using labelled superpixels
[am] = regionadjacency(ls, 8);

% Abstracting each superpixel to calculate the average Lab and position xy
lab_image = rgb2lab(im);    % Corresponding Lab image of the input image
[X,Y] = meshgrid(1:width, 1:height);

L_channel = lab_image(:,:,1);
a_channel = lab_image(:,:,2);
b_channel = lab_image(:,:,3);
for i = 1:n
    mask = ls==i;
    nm = sum(mask(:));  
    sp(i).L = sum(L_channel(mask))/nm;
    sp(i).a = sum(a_channel(mask))/nm;
    sp(i).b = sum(b_channel(mask))/nm;
    
    sp(i).r = sum(Y(mask))/nm;
    sp(i).c = sum(X(mask))/nm;
end
end