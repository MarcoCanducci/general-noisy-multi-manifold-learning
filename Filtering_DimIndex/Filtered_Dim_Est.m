function [B,LNew] = Filtered_Dim_Est(Data,Radius)
% function LNew = Filtered_Dim_Est(Data,Radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% FILTERED_DIM_EST_NEW determines the intrinsic dimensionality of a cloud % 
% of points.                                                              %
%                                                                         %
%  Copyright (C) 2020  Marco Canducci                                     %
%  Email: marco.canducci91@gmail.com                                      %
%                                                                         %
%  This program is free software: you can redistribute it and/or modify   %
%  it under the terms of the GNU Affero General Public License as         %
%  published by the Free Software Foundation, either version 3 of the     %
%  License, or (at your option) any later version.                        %
%                                                                         %
%  This program is distributed in the hope that it will be useful,        %
%  but WITHOUT ANY WARRANTY; without even the implied warranty of         %
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          %
%  GNU Affero General Public License for more details.                    %
%                                                                         %
% You should have received a copy of the GNU Affero General Public License%
% along with this program.  If not, see <https://www.gnu.org/licenses/>.  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Description
%
%   Input parameters:
%       -   Data_Original: Data not yet processed using SAF.
%       -   Data: Matrix Nx3 containing the observed points in a 3D space.
%                 N is the number of observations, 3 is the dimensionality 
%                 of the embedding space.
%       -   Radius: Radius for the neighbouring search.
%
%   Output parameters:
%       -   LNew: Matrix containing eigenvalues associated to each  
%           obs' neighbourhood, renormalized by their sum. 
%       -   NewLabels: Logic array containing "1" for the valid obs and
%           "0" for the rejected ones.
%       -   NewData: Matrix containing the valid obs of which a
%           dimensionality index has been found.
%       -   idx: column vector with dimensionality index for each valid
%           observation.
%       -   C: 3x3 matrix with the new centres after kmedoids.
%       -   DistM: N_Valid_obs x 3 matrix with the distances of each valid
%           obs to the initial clustering centres.
%   
%   LNew = Dim_Est_New(Data,Radius) gives the matrix containing the
%       renormalized eigenvectors.
%   [LNew,NewLabels] = Dim_Est_New(Data,Radius) also returns the labels
%       for the selection of valid obs;
%   [LNew,NewLabels,NewData] = Dim_Est_New(Data,Radius) also returns the
%       valid obs
%   [LNew,NewLabels,NewData, idx] = Dim_Est_New(Data,Radius) returns the
%       dimensionality index associated to each obs in Data.
%   [LNew,NewLabels,NewData, idx, C] = Dim_Est_New(Data,Radius) gives the
%       centres after the clustering through k-medoids is performed.
%   [LNew,NewLabels,NewData, idx, C, DistM] = Dim_Est_New(Data,Radius)
%       returns the matrix with all distances of the obs to the centres as
%       computed with FisherMetricDist.
%
%   For each obs in Data, the covariance matrix of its neighbourhood is
%   decomposed into its eigenvalues. These get then normalized by their sum,
%   so that each triade of eigenvalues lies on a portion of the
%   3d-simplex. A further step is applied for distributing the eigenvalues
%   in the complete simplex ([1 0 0, 0 1 0; 0 0 1]) so that the its mean is
%   equidistant from each vertex and coincides with the point [1/3 1/3
%   1/3]. Here probabilities vary between [0 1] for each one of the
%   dimensions. The distance of each triade of new rescaled eigenvaòues
%   to the vertices is computed and stored in DistM, while the MATLAB
%   implementation of the algorithm k-medoids is used for clustering obs
%   with respect to the centres (used as clusters' prototypes).

Neighbors2 = rangesearch(Data,Data,Radius);
LNew = zeros(size(Data));
disp(size(Neighbors2));
k = 1;
for i=1:length(Data)
    % By evaluating neighborhoods in the original dataset we can discard
    % all points which didn't have neighbors populated enough to survive
    % before appyling SAF. This step removes artifacts of the SAF
    % algorithm.
        [~,~,Lambda2] = pca(Data(Neighbors2{i},:));
        Lambda2 = sort(Lambda2,'descend');

        LambdaTot = zeros(1,size(Data,2));
        LambdaTot(1:size(Lambda2,1)) = Lambda2;
        LambdaTot = LambdaTot./sum(LambdaTot);
        LNew(i,:) = LambdaTot;
        if(any(isnan(LambdaTot))==1)
            disp(i);
            C(k) = i;
            k = k+1;
        end        
end

P = [1 0 0; 0.5 0.5 0; 1/3 1/3 1/3];
T = [1 2 3];
TR = triangulation(T,P);
B = zeros(size(LNew));

% Rescaling of the eigenvalues using the barycentrc coordinates with
% respect to the portion of simplex of vertices [1 0 0; 0.5 0.5 0; 1/3 1/3
% 1/3].

for i=1:size(B,1)
    B(i,:) = cartesianToBarycentric(TR,1,LNew(i,:));
end