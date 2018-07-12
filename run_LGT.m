%% This code demonstrates the proposed temporal region detection algorithm LGT
% This code reads images from 'Input' folder, computes the temporal saliency and
% saves the saliency maps to 'Output' folder
clear all;clc;
fprintf('Input images are in the directory "Input". \n');
input_dir = dir('Input\');
tic;
im1 = imread('Input\frame_51.png'); % Video frame 1
im2 = imread('Input\frame_52.png'); % Video frame 2
height = size(im1,1); width = size(im1,2);
% Resize the image into maximum dimension of 400
if(max(size(im1))~=400)
   im1 = imresize(im1,(400/max(size(im1))));    im2 = imresize(im2,(400/max(size(im2))));
end
[ls, am, sp] = patchSuperpixel(im1, 500);   % Superpixel segmentation using SLIC superpixels
addpath('Optical Flow');
[u,v] = Flow(im1,im2); % Optical flow estimation
[salMap] = temporalSaliency(im1,im2,ls,am,u,v);   % Temporal Saliency estimation
salMap = imresize(salMap,[height,width]);
imwrite(salMap, 'Output\frame_51_LGT.png');  % Writig the output saliency maps
toc;
fprintf('Saliency Map can be found in the directory "Output". \n');