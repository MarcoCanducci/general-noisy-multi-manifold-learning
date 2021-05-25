function [M,Noisy_Man,NoisyData,finaldata] = MultiManifold_UnDens(Man_ID,NoiseDSize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% MULTIMANIFOLD_UNDENS builds the toy data set described in the paper.     %
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
%   Man_ID is a row vector containing the ID for the selected manifolds:
%       ID = 1: Parab arm (1D);
%       ID = 2: Spiral arm (1D);
%       ID = 3: S - surface (2D);
%       ID = 4: Bow shock/ Hyperbolic surface (2D);
%       ID = 5: Empty toroid (2D);
%       ID = 6: Mobius strip (2D);
%       ID = 7: Gaussian point cloud (3D);
%
%   finaldata is the dataset composed by the noisy manifolds and the
%   background noise with double the number of points then those of the
%   union of all selcted manifolds; Noise is centered so that the domain of
%   the datasets contains widely all manifolds. The width of the noise
%   can be changed via the parameter "eps", while the percentage of noise via
%   "percentage".
%


if max(Man_ID>7)
%     disp('Only 7 types of manifold are available. The ID in Man_ID cannot be greater than 7.')
    error('Man_ID(:) > 7')
end
N_points = 5000;
X1 = 0.9*rand(N_points,1);
Y1 = (-2*X1.^2 + 2.*X1+1);
Z1 = zeros(N_points,1);
Parab = [X1 Y1 Z1-0.2];
Parab(:,1) = Parab(:,1) - 0.1;
Parab(:,2) = Parab(:,2) + 0.1;

M{1} = Parab;

clearvars Parab X1 Y1 Z1

t1 = 24*rand(N_points,1) + 7;
t1=sort(t1);
Radius = linspace(0.1,1.,N_points)';
X2 = (Radius.*cos(t1)+1)/2 -0.2;
Y2 = (Radius.*sin(t1)+1)/2;
Z2 = (linspace(2,1,N_points)'+2);
Spiral = [X2 Y2 Z2];

Spiral_Shift = zeros(size(Spiral,1),3);

Spiral_Shift(:,1) = 0.1;
Spiral_Shift(:,2) = 0.4;
Spiral_Shift(:,3) = -0.9;

Spiral = Spiral.*0.4 + Spiral_Shift;

M{2} = Spiral;

clearvars Spiral Spiral_Shift X2 Y2 Z2 Radius

R_SShape = [1 0 0; 0 0 -1; 0 1 0];
SShape = Surf_SShape(N_points);

SShape = SShape*R_SShape;

SShape_Shift = zeros(size(SShape,1),3);
SShape_Shift(:,3) = SShape_Shift(:,3) + 0.5;

SShape = SShape + SShape_Shift;

M{3} = SShape;

clearvars SShape SShape_Shift R_SShape

B = Surf_Hyperbole(N_points);

B_Shift = zeros(size(B,1),3);

B_Shift(:,1) = -0.4;
B_Shift(:,2) = 1.3;
B_Shift(:,3) = 0.3;

B = 1.2.*B + B_Shift;

M{4} = B;

clearvars B B_Shift 

T = Surf_Toroid(N_points);

T_Shift = zeros(size(T,1),3);
T_Shift(:,1) = 0.2;
T_Shift(:,2) = 0.68;
T_Shift(:,3) = 0.5;

% T = (T + T_Shift).*1.1;
T = T.*1.2 + T_Shift;

M{5} = T;

clearvars T T_Shift

Mob = (Surf_Mobius(N_points).*0.25)*[1 0 0; 0 0 -1; 0 1 0]*[0 -1 0; 1 0 0; 0 0 1];
% Mob = Surf_Mobius(N_points);
Mob(:,1) = Mob(:,1) + 0.4;
Mob(:,2) = Mob(:,2) + 1.5;
Mob(:,3) = Mob(:,3) - 0.2;

M{6} = Mob;

clearvars Mob

G = randn(N_points,3)./20;

G(:,1) = G(:,1) - 0.2;
G(:,2) = G(:,2) - 0.5;

M{7} = G;

clearvars G

% U_Data(:,1) = -1 + 3.5*rand(1e6,1);
% U_Data(:,2) = -1 + 4.5*rand(1e6,1);
% U_Data(:,3) = -1 + 2*rand(1e6,1);
U_Data(:,1) = -1.1 + 2.5*rand(NoiseDSize,1);
U_Data(:,2) = -0.9 + 2.7*rand(NoiseDSize,1);
U_Data(:,3) = -1.1 + 2.5*rand(NoiseDSize,1);

% Note that now data on many manifold is not distributed uniformly. To make
% a uniform distribution we produce a densly sampled uniform distribution
% and choose all point in this distribution which are near to point on
% manifold in a small radius, so we are completely sure that the
% distribution on the manifold is uniform. We replace new point instead of
% manifold points.
for j = 1:size(Man_ID,2)
    NN = rangesearch(U_Data,M{Man_ID(j)},0.04);
    Noisy_Man{j} = U_Data(NN{1},:);
    for i = 2:size(NN,1)
        Noisy_Man{j} = union(Noisy_Man{j},U_Data(NN{i},:),'rows','stable');
    end
    disp(j)
    clearvars NN
    Noisy_Man{j} = unique(Noisy_Man{j},'rows','stable');
    size(Noisy_Man{j})
end
% figure;
% plot3(Noisy_Man{1}(:,1),Noisy_Man{1}(:,2),Noisy_Man{1}(:,3),'.')
% hold on
% plot3(Noisy_Man{2}(:,1),Noisy_Man{2}(:,2),Noisy_Man{2}(:,3),'.')
% plot3(Noisy_Man{3}(:,1),Noisy_Man{3}(:,2),Noisy_Man{3}(:,3),'.')
% plot3(Noisy_Man{4}(:,1),Noisy_Man{4}(:,2),Noisy_Man{4}(:,3),'.')
% plot3(Noisy_Man{5}(:,1),Noisy_Man{5}(:,2),Noisy_Man{5}(:,3),'.')
% plot3(Noisy_Man{6}(:,1),Noisy_Man{6}(:,2),Noisy_Man{6}(:,3),'.')
% plot3(Noisy_Man{7}(:,1),Noisy_Man{7}(:,2),Noisy_Man{7}(:,3),'.')

if size(Man_ID,2) == 1
    NoisyData = Noisy_Man{1};
elseif size(Man_ID,2) == 2
    NoisyData = union(Noisy_Man{1},Noisy_Man{2},'rows','stable');
else
    NoisyData = union(Noisy_Man{1},Noisy_Man{2},'rows','stable');
    for i =3:size(Man_ID,2)
        NoisyData = union(NoisyData,Noisy_Man{i},'rows','stable');
    end
end
eps = 0.2;
percentage = 200;
Max = max(NoisyData);
Min = min(NoisyData);
Max = Max + eps.*(Max-Min);
Min = Min - eps.*(Max-Min);
% Mid = (Max + Min).*0.5;

backgroundnoise = zeros((percentage/100)*size(NoisyData,1),size(NoisyData,2));
backgroundnoise(:,1) = (1.2*(Max(1)-Min(1)).*rand((percentage/100)*size(NoisyData,1),1)+Min(1))';
backgroundnoise(:,2) = (1.2*(Max(2)-Min(2)).*rand((percentage/100)*size(NoisyData,1),1)+Min(2))';
backgroundnoise(:,3) = (1.2*(Max(3)-Min(3)).*rand((percentage/100)*size(NoisyData,1),1)+Min(3))';

% figure; plot3(NoisyData(:,1),NoisyData(:,2),NoisyData(:,3),'.','markersize',5)
%making noise for the background

% backgroundnoise = 2*rand((percentage/100)*size(NoisyData,1),size(NoisyData,2))-1;


finaldata = [NoisyData;backgroundnoise];
% finaldata = NoisyData;