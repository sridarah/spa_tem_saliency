%% This code is used to cluster image patches using spectral clustering
% Input : 
%       ls - labelled patches in image
%       sp - patch color and location
%       am - adjacency matrix of neighbouring patches
% Output :
%       lc - labelled clusters in image
%       csp - clustered patches
%       sc - labels of patches
%
% Thanks to Asad Ali for their code to implement spectral clustering algorithm in the paper
% Ng, A., Jordan, M., and Weiss, Y. (2002). On spectral clustering: analysis and an algorithm. In T. Dietterich,
% S. Becker, and Z. Ghahramani (Eds.), Advances in Neural Information Processing Systems 14 
% (pp. 849 – 856). MIT Press.
%% Copyright (c) 2010, Asad Ali
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%
% Ng, A., Jordan, M., and Weiss, Y. (2002). On spectral clustering: analysis and an algorithm. In T. Dietterich,
% S. Becker, and Z. Ghahramani (Eds.), Advances in Neural Information Processing Systems 14 
% (pp. 849 – 856). MIT Press.
% Asad Ali
% GIK Institute of Engineering Sciences & Technology, Pakistan
% Email: asad_82@yahoo.com
% CONCEPT: Introduced the normalization process of affinity matrix(D-1/2 A D-1/2), 
% eigenvectors orthonormal conversion and clustering by kmeans 

function [lc, csp, sc] = regionSpectralClustering(ls, sp, am, k_min, k_max)
% Affinity Matrix Calculation
% read color values of superpixels
no_of_sp = length(sp);    % Number of superpixels
super_pixel_color = zeros([no_of_sp,3]); % Creating a vector to represent color each superpixel (l,a,b)
for i = 1 : no_of_sp
    super_pixel_color(i,1) = sp(1,i).L;    super_pixel_color(i,2) = sp(1,i).a;     super_pixel_color(i,3) = sp(1,i).b;
end;
% initialize affinity matrix
affinity = zeros([no_of_sp,no_of_sp]);
% read ajdacency matrix
no_of_adjsp = size(nonzeros(am),1);     % number of nonzero elements in the adjacency matrix
[a,b] = find(am);   % store indices of adjacent superpixels
% calculate affinity matrix
for i = 1 : no_of_adjsp
     affinity(a(i),b(i)) = norm(super_pixel_color(a(i),:) - super_pixel_color(b(i),:)); % color distance
end
% Distance normalization
sigma = 0.4 ; % weight for distance to similairy conversion
max_dis=max(affinity); min_dis=min(affinity);
for i = 1 : no_of_adjsp
     affinity(a(i),b(i)) = (affinity(a(i),b(i))-min_dis)/(max_dis-min_dis); % color distance normalization
     affinity(a(i),b(i)) = exp(-affinity(a(i),b(i))^2/(2*sigma^2)); % conversion of distance to similarity
end
% Clustering superpixels using Jordan et al 2002
% compute the degree matrix
for i=1:size(affinity,1)
    D(i,i) = sum(affinity(i,:));
end

% compute the normalized laplacian / affinity matrix (method 1)
%NL1 = D^(-1/2) .* L .* D^(-1/2);
for i=1:size(affinity,1)
    for j=1:size(affinity,2)
        NL1(i,j) = affinity(i,j) / (sqrt(D(i,i)) * sqrt(D(j,j)));  
    end
end

% compute the normalized laplacian (method 2)  eye command is used to
% obtain the identity matrix of size m x n
% NL2 = eye(size(affinity,1),size(affinity,2)) - (D^(-1/2) .* affinity .* D^(-1/2));

% perform the eigen value decomposition
[eigVectors,eigValues] = eig(NL1);

% choosing value of k
% calculate the spectrum of the Laplacian and sort the eigenvalues in ascending order
v=diag(eigValues);
[vs, is] = sort(v,'ascend');

% Here we implement the eigengap heuristic
eigengaps = zeros(length(vs)-1,1);
for i=1:1:length(eigengaps)
    if ((i<k_min) || (i> k_max))
        eigengaps(i)=-1;
    else
        eigengaps(i)=vs(i+1)-vs(i);
    end
end
[junk k] = max(eigengaps); % choosing the value of k

% select k largest eigen vectors
nEigVec = eigVectors(:,(size(eigVectors,1)-(k-1)): size(eigVectors,1));

% construct the normalized matrix U from the obtained eigen vectors
for i=1:size(nEigVec,1)
    n = sqrt(sum(nEigVec(i,:).^2));    
    U(i,:) = nEigVec(i,:) ./ n; 
end
warning('off');
% perform kmeans clustering on the matrix U
[sc,C] = kmeans(U,k,'start','sample','emptyaction','singleton','Replicates',100);
% assigning labels to the superpixel clusters
lc=ls;
for i = 1 : length(sc)
lc(lc == i) = sc(i);
end
% save the superpixels in each cluster using cell array
csp = {};
for i = 1 : k
    [temp] = find(sc==i);
    temp = transpose(temp); % inversion to make it row matrix
    csp(i) = {temp};
end
end
