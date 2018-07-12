function [vx,vy] = Flow(im1,im2)

addpath('Optical Flow\mex');

% load the two frames
im1 = im2double(im1);
im2 = im2double(im2);

% set optical flow parameters (see Coarse2FineTwoFrames.m for the definition of the parameters)
alpha = 0.012;
ratio = 0.75;
minWidth = 20;
nOuterFPIterations = 7;
nInnerFPIterations = 1;
nSORIterations = 30;

para = [alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations];

% this is the core part of calling the mexed dll file for computing optical flow
% it also returns the time that is needed for two-frame estimation
[vx,vy] = Coarse2FineTwoFrames(im1,im2,para);
end