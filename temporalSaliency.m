%% This code implements temporal saliency detection using patch level optical flow abstraction
% Input:
%       im1 - Video frame 1
%       im2 - Video frame 2
%       ls - location of superpixels in the image
%       am - adjacency matrix of superpixels
%       u - Horizondal motion
%       v - Vertical motion
% Output:
%       temporal_saliency_map - temporal saliency map
%%
function [temporal_saliency_map] = temporalSaliency(im1,im2,ls,am,u,v)
%% Magnitude and flow Calculation of optical flow vectors
px_flow_magnitude   = (u.^2 + v.^2 ).^0.5;  % Magnitude of optical flow vector of each pixel
px_flow_orientation =  atan2d(v,u);  % Orientation of optical flow vector of each pixel (in degrees)

% Orientation quantization
orientation_quantized = px_flow_orientation;    % Matrix containing the quantized orientation
% Assigning each pixel an orientation label
orientation_quantized(orientation_quantized<=22.5 & orientation_quantized>0) = 1; % orientation 1
orientation_quantized(orientation_quantized<=0 & orientation_quantized>-22.5) = 1; % orientation 1
orientation_quantized(orientation_quantized<=67.5 & orientation_quantized>22.5) = 2; % orientation 2
orientation_quantized(orientation_quantized<=112.5 & orientation_quantized>67.5) = 3; % orientation 3
orientation_quantized(orientation_quantized<=157.5 & orientation_quantized>112.5) = 4; % orientation 4
orientation_quantized(orientation_quantized<=180 & orientation_quantized>157.5) = 5; % orientation 5
orientation_quantized(orientation_quantized>=-180 & orientation_quantized<-157.5) = 5; % orientation 5
orientation_quantized(orientation_quantized>=-157.5 & orientation_quantized<-112.5) = 6; % orientation 6
orientation_quantized(orientation_quantized>=-112.5 & orientation_quantized<-67.5) = 7; % orientation 7
orientation_quantized(orientation_quantized>=-67.5 & orientation_quantized<=-22.5) = 8; % orientation 8
orient_lebel_degree = [0, 45, 90, 135, 180, 225, 270, 315]; % Corresponding degree values of each orientation
%% Dominant Flow Estimation at Superpixel Level

no_of_sp =  length(am);    % Number of superpixels
sp_flow_magnitude = zeros([no_of_sp,1]);    % Magnitude of dominant flow of superpixel
sp_flow_orientation = zeros([no_of_sp,1]);  % Orientation of dominant flow of superpixel

sp_overall_magnitude = zeros([8,1]);    % Temporary variable stores magnitude of each orientation of superpixel
for i = 1 : no_of_sp
   for j = 1 : 8
       sp_overall_magnitude(j) = sum(px_flow_magnitude(ls==i & orientation_quantized == j));  % Single & is used here caz ls and orientation_quantized are tow matrices
   end
   sp_flow_orientation(i) = find(sp_overall_magnitude==max(sp_overall_magnitude));  % Orientation of the dominant flow
   sp_flow_magnitude(i) = sum(px_flow_magnitude(ls==i & orientation_quantized == sp_flow_orientation(i)));    % Magnitude of dominant flow vectors
   sp_flow_orientation(i) = orient_lebel_degree(sp_flow_orientation(i));    % Assigning degree for each orientation
end
% Normalizing magnitude of each superpixel
max_mag = max(sp_flow_magnitude); min_mag = min(sp_flow_magnitude);
sp_flow_magnitude(:) = (sp_flow_magnitude(:)-min_mag)/(max_mag-min_mag);

%% Local Temporal Saliency Estimation using multi-level center-surround motion contrast
surround_level = 3; % Number of surrounding levels for local temporal saliency estimation
local_temporal_saliency = zeros([no_of_sp,1]);
% Finding neighbours of each superpixel
[nei sp_id] = find(am);   % 'nei' stores neighbours and 'sp_id' stores superpixel ids
for i = 1 : no_of_sp
    sp_nei = {};    % Cell containing multi-level neighbours
    arr = nei(sp_id==i);    % First level neighbours
    sp_nei{1,1} = arr;    % Storing first level neighbours in cell
    for j = 1 : (surround_level-1)
        arr = sp_nei{1,j};  % Superpixels whose neighbours to be found in j+1th level
        arr_nei = [];
        for k = 1 : length(arr)
            surr_sps = nei(sp_id==arr(k));  % Neighbour of each superpixel
            arr_nei = [arr_nei;surr_sps];   % Concatenate all the neighbours of each level of superpixel
        end
        arr_nei = unique(arr_nei);  % Unique neighbours
        for l = 1 : j % loop to remove superpixels from previous levels
            arr = sp_nei{1,l};
            arr_nei(ismember(arr_nei,arr))=[];
        end
        sp_nei{1,j+1} = arr_nei;    % Neighbours of ith superpixel at all the levels
    end
    % Calculation of local rarity of temporal saliency
    tmp_tempsal = 0;
    for j = 1 : surround_level
        surr_sps = sp_nei{1,j};
        surr_length = length(surr_sps);
        for k = 1 : surr_length
            norm_orient_contrast  = abs(sp_flow_orientation(i)-sp_flow_orientation(surr_sps(k)));   % Absolute angular difference
            if(norm_orient_contrast>180)
                norm_orient_contrast = 360-norm_orient_contrast;
            end
            norm_orient_contrast = norm_orient_contrast/180;    % Normalized angular difference
            tmp_tempsal = tmp_tempsal+((1/surround_level)*(1/surr_length)*(norm_orient_contrast));
        end
    end
    local_temporal_saliency(i) = sp_flow_magnitude(i)*tmp_tempsal;
end
% Normalizing local temporal saliency
max_lctm = max(local_temporal_saliency); min_lctm = min(local_temporal_saliency);
local_temporal_saliency(:) = (local_temporal_saliency(:)-min_lctm)/(max_lctm-min_lctm);
%% Global Temporal Saliency Estimation using self-information of flow orientation

global_temporal_saliency = zeros([no_of_sp,1]);
for i = 1 : no_of_sp
    global_temporal_saliency(i) = sp_flow_magnitude(i)*(-log(length(find(sp_flow_orientation==sp_flow_orientation(i)))/no_of_sp)); % Global temporal saliency computation
end
% Normalizing global temporal saliency
max_gltm = max(global_temporal_saliency); min_gltm = min(global_temporal_saliency);
global_temporal_saliency(:) = (global_temporal_saliency(:)-min_gltm)/(max_gltm-min_gltm);
%% Temporal Saliency Aggregation
temporal_saliency  = zeros([no_of_sp,1]);
temporal_saliency(:) = local_temporal_saliency.*global_temporal_saliency;
temporal_saliency(:) = (temporal_saliency(:)-min(temporal_saliency))/(max(temporal_saliency)-min(temporal_saliency));
%% Patch level Temporal Saliency Assignment
local_temporal_map = ls;  global_temporal_map = ls;  % Initializing contrast and distribution map
for i = 1 : no_of_sp
    local_temporal_map(local_temporal_map == i) = local_temporal_saliency(i); % Saliency assignment to superpixels
    global_temporal_map(global_temporal_map == i) = global_temporal_saliency(i); % Saliency assignment to superpixels    
end
temporal_saliency_map = ls; % Initializing temporal saliency map
for i = 1 : no_of_sp
    temporal_saliency_map(temporal_saliency_map == i) = temporal_saliency(i); % Saliency assignment to superpixels  
end
end