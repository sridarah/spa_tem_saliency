%% This code implements saliency detection using patch and region image abstractions
%
% Input:    sp - Patch attribute structure array with fields:
%                   L  - Mean L* value
%                   a  - Mean a* value
%                   b  - Mean b* value
%                   r  - Mean row value
%                   c  - Mean column value
%           csp - clustered superpixels
%           ls - labelled superpixels in image
%           sc - labels of patches
%           am - adjacency matrix of neighbouring superpixels
%
% Output:   saliency_map - Saliency map
%           
function [saliency_map] = spatialSaliency(sp, csp, ls, sc, am)
no_of_sp = length(sp);    % Number of superpixels
no_of_reg = length(csp); % Number of regions
superpixel_color = zeros([no_of_sp,3]); % Creating a vector to represent  color each superpixel (l,a,b)
superpixel_location = zeros([no_of_sp,2]); % Creating a vector to represent  color each superpixel (x,y)

% Assigning color values and location to the superpixel vector
for i = 1 : no_of_sp
    superpixel_color(i,1) = sp(1,i).L;    superpixel_color(i,2) = sp(1,i).a;     superpixel_color(i,3) = sp(1,i).b;
    superpixel_location(i,1) = sp(1,i).r;     superpixel_location(i,2) = sp(1,i).c;
end;

% Region prototyping using mediation
% Color and location of region prototypes
region_proto_color = zeros([no_of_reg,3]);
region_proto_location = zeros([no_of_reg,2]);

region_color_sp = int32(no_of_reg);   % indicates the color of prototype superpixel of a region
region_location_sp = int32(no_of_reg);   % indicates the center of prototype superpixel of a region
for i = 1 : no_of_reg
    col_distance = zeros([length(csp{1,i}),length(csp{1,i})]);
    loc_distance = zeros([length(csp{1,i}),length(csp{1,i})]);
    col_dist_selection = zeros([length(csp{1,i}),2]);
    loc_distance_selection = zeros([length(csp{1,i}),2]);
    % calculate distance
    for j = 1 : length(csp{1,i})
        for k = 1 : length(csp{1,i})
        col_distance(j,k) = norm(superpixel_color(csp{1,i}(j),:) - superpixel_color(csp{1,i}(k),:)); % color distance
        loc_distance(j,k) = norm(superpixel_location(csp{1,i}(j),:) - superpixel_location(csp{1,i}(k),:)); % location distance
        end
    end
    temp1=sum(col_distance,2); % Summing up the rows to calculate the overall color distance of each superpixel
    cum_contrast=sum(loc_distance,2); % Summing up the rows to calculate the overall location distance of each superpixel
    for j = 1 : length(csp{1,i})
        col_dist_selection(j,1) = csp{1,i}(j);
        col_dist_selection(j,2) = temp1(j);
        loc_distance_selection(j,1) = csp{1,i}(j);
        loc_distance_selection(j,2) = cum_contrast(j);
    end
    col_dist_selection = sortrows(col_dist_selection, 2); % Sorting in acending order to find the superpixel with minimum color distance
    region_color_sp(i)=int32(col_dist_selection(1,1));   % Select color prototype superpixel that is the first element
    loc_distance_selection = sortrows(loc_distance_selection, 2); % Sorting in acending order to find the superpixel with minimum location distance
    region_location_sp(i)=int32(loc_distance_selection(1,1));   % Select location prototype superpixel that is the first element
end
for i = 1 : no_of_reg
    region_proto_color(i,1)=superpixel_color(region_color_sp(i),1); % Color L
    region_proto_color(i,2)=superpixel_color(region_color_sp(i),2); % Color a
    region_proto_color(i,3)=superpixel_color(region_color_sp(i),3); % Color b
    region_proto_location(i,:)=superpixel_location(region_location_sp(i),:); % Location Coordinates [x,y]
end

% Color Contrast Cue Estimation
sp_rg_global_contrast = zeros([no_of_sp, no_of_reg-1]);
sp_rg_global_distance = zeros([no_of_sp, no_of_reg-1]);
global_contrast = zeros([no_of_sp,1]);

for i = 1 : no_of_sp
    k = 1;
    for j = 1 : no_of_reg
        if(sc(i)~=j)
            sp_rg_global_contrast(i,k) = norm(superpixel_color(i,:) - region_proto_color(j,:)); % color contrast to other region centers
            sp_rg_global_distance(i,k) = norm(superpixel_location(i,:) - region_proto_location(j,:)); % Distance to other region centers
            k = k+1;
        end
    end
end
% Normalizing the global contrast and distance values
    cum_contrast = zeros([no_of_sp,1]);
    minval=min(sp_rg_global_contrast);   maxval=max(sp_rg_global_contrast);
    maxdim=max(size(ls));   % Maximum dimension of an image
beta1 = 2;
for i = 1 : no_of_sp
    for j = 1 : no_of_reg-1
        sp_rg_global_contrast(i,j) = (sp_rg_global_contrast(i,j)-minval)/(maxval-minval);
        sp_rg_global_distance(i,j) = sp_rg_global_distance(i,j)/maxdim;
        temp1(i,j) = (length(csp{1,j})/no_of_sp)*(exp(-sp_rg_global_distance(i,j)*beta1))*sp_rg_global_contrast(i,j);
    end
end
    cum_contrast=sum(temp1,2); % Cumulative contrast of each superpixel
    global_contrast= cum_contrast;

global_contrast(:) = (global_contrast(:)-min(global_contrast))/(max(global_contrast)-min(global_contrast));
% Color Distribution Cue Estimation
sp_rg_mean_position =  zeros([no_of_sp, 2]); % Store mean position of color of a superpixel
global_distribution = zeros([no_of_sp,1]);

% Compute mean position of each superpixel wrt region prototypes
beta2=8; tempx=[]; tempy=[];
for i = 1 : no_of_sp
    k = 1;
    for j = 1 : no_of_reg
        if(sc(i)~=j)
        tempx(k) = (exp(-sp_rg_global_contrast(i,k)*beta2))*region_proto_location(j,1);
        tempy(k) = (exp(-sp_rg_global_contrast(i,k)*beta2))*region_proto_location(j,2);
        k = k+1;
        end
    end
    sp_rg_mean_position(i,1) = sum(tempx(:))/(no_of_reg-1); % Mean x position
    sp_rg_mean_position(i,2) = sum(tempy(:))/(no_of_reg-1); % Mean y position
end

% Compute color distribution of each superpixel wrt region prototypes
temp = zeros(no_of_sp,no_of_reg-1);
for i = 1 : no_of_sp
    k = 1;
    for j = 1 : no_of_reg
        if(sc(i)~=j)
        temp(i,k) = (length(csp{1,j})/no_of_sp)*(norm(region_proto_location(j,:)-sp_rg_mean_position(i,:))/maxdim)*(exp(-sp_rg_global_contrast(i,k)*beta2)); % distribution
        k=k+1;
        end
    end
end
global_distribution = sum(temp,2); % Cumulative color distribution of each superpixel

global_distribution(:) = (global_distribution(:)-min(global_distribution))/(max(global_distribution)-min(global_distribution));
global_distribution(:) = (1-global_distribution(:));
% Saliency Refinement
[nei sp_id] = find(am);   % 'nei' stores neighbours and 'sp_id' stores superpixel ids
    mu = 0.2;
    max_thr = 1-mu; min_thr = mu;

% Superpixel Saliency Assignment
saliency_map = ls;  % Initializing saliency map
global_saliency = zeros([no_of_sp,1]);
global_saliency = global_contrast.*global_distribution; % Fusion
% Normalization
global_saliency(:) = (global_saliency(:)-min(global_saliency))/(max(global_saliency)-min(global_saliency));

% Saliency refinement
refined_global_saliency = zeros([no_of_sp,1]);
for i = 1 : no_of_sp
    avg_nei_sal = mean(global_saliency(nei(sp_id==i)));
    
    if(avg_nei_sal>=max_thr && global_saliency(i)<=max(global_saliency(nei(sp_id==i))))
        refined_global_saliency(i) = max(global_saliency(nei(sp_id==i)));
    elseif(avg_nei_sal<=min_thr && global_saliency(i)>=min(global_saliency(nei(sp_id==i))))
        refined_global_saliency(i) = min(global_saliency(nei(sp_id==i)));
    else
        refined_global_saliency(i) = global_saliency(i);
    end
end
global_saliency = refined_global_saliency; % Refined saliency

% Incorporation of center prior map
height = size(ls,1); width = size(ls,2);
superpixel_center_prior = zeros([no_of_sp,1]);
center_position = [height/2,width/2];
sigma2 = min(height,width)/2.5;
for i = 1:no_of_sp
        d = norm(superpixel_location(i,:)-center_position);
        superpixel_center_prior(i) = exp(-d^2/(2*sigma2^2));
end
global_saliency = global_saliency.*superpixel_center_prior;
minval = min(global_saliency); maxval = max(global_saliency);
global_saliency(:) = (global_saliency(:)-minval)/(maxval-minval);

% Saliency assignment
for i = 1 : no_of_sp  
    saliency_map(saliency_map==i) = global_saliency(i);
end
end