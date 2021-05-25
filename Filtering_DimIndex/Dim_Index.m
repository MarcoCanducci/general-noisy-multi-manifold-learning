function [Struct,idx] = Dim_Index(SAF_D,r,Simplex,smooth)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Dim_Index, assigns to every point in a data set a dimensionality index. %
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

%Find rescaled eigenvectors using barycentric cohordinates (B) and
%denoised eigenvalues.
r_Idx = r;
NewData = SAF_D;
[~,d] = size(SAF_D);
[B,L] = Filtered_Dim_Est(SAF_D,r_Idx);
% L = Filtered_Dim_Est(SAF_D,r_Idx);
Struct.Eig = L;
Struct.Baryc_Eig = B;
idx = zeros(size(SAF_D,1),4);


switch Simplex    
    case 'Original'
        % Definition of Original simplex (portion) vertices and equidistant point. 
        Eig = abs(L);
        V = [1 0 0; 0.5 0.5 0; 1/3 1/3 1/3];
        for i=1:d
            V(i,1:i) = 1/i;
        end
        m0 = mean(V);
        Ain = [-1 1 0; 0 -1 1; 0 0 -1];
        bin = [0; 0; 0];
        Aeq = [1 1 1; 0 0 0; 0 0 0];        
        beq = [1; 0; 0];
        m = patternsearch(ErrF,m0,Ain,bin,Aeq,beq);
    case 'Barycentric'
        % Definition of whole simplex vertices and equidistant points. 
        
        Eig = abs(B);
        V = [1 0 0; 0 1 0; 0 0 1];

        m = mean(V);
        
end
% First index: kmedoids on eigenvalues using vertices as propototypes.

% Index as smallest distance from vertex
disp('Computing second Index')
scale = pdist2(V,m,@FisherMetricDist);

Dist = pdist2(Eig,V,@FisherMetricDist);

[~,ID2] = min(Dist,[],2);
idx(:,1) = ID2;


Struct.Idx1 = idx(:,1);
Struct.Idx2 = idx(:,2);

% l1 and l2 kernels computation.
disp('Computing third and fourth Index')
K1 = exp(-Dist./(2.*scale(1)));
K2 = exp(-(Dist.^2)./(2.*scale(1)^2));


K1 = K1./repmat(sum(K1,2),1,size(K1,2));
K2 = K2./repmat(sum(K2,2),1,size(K2,2));

% Entropies computation

log3P1 = log(K1)./repmat(log(3),size(K1,1),size(K1,2));
H1 = -sum(K1.*log3P1,2);

log3P2 = log(K2)./repmat(log(3),size(K2,1),size(K2,2));
H2 = -sum(K2.*log3P2,2);

T1 = 1;

DNewH1 = NewData(H1<=T1,:);
EigH1 = Eig(H1<=T1,:);
[~,H1_Idx] = max(K1(H1<=T1,:),[],2);


idx(:,2) = H1_Idx;
Struct.Kernel1_Data = DNewH1;
Struct.Kernel1_Idx = H1_Idx;

T2 = 1;
DNewH2 = NewData(H2<=T2,:);
EigH2 = Eig(H2<=T2,:);
[~,H2_Idx] = max(K2(H2<=T2,:),[],2);


idx(:,3) = H2_Idx;
Struct.Kernel1_Data = DNewH2;
Struct.Kernel1_Idx = H2_Idx;

% Smoothing kernels on data space;
disp('Computing fifth Index')
switch smooth
    case 'l1'
        K = K1;
    case 'l2'
        K = K2;
end

r_smooth = 2*r;
K_scale2 = 2*(r_smooth.^2);
[NN,Dist] = rangesearch(NewData,NewData,r_smooth);
S = zeros(size(K));

for i = 1:size(NN,1)
    W = exp(-(Dist{i}.^2)./K_scale2);
    W2 = repmat(W',1,3).*K(NN{i},:);
    S(i,:) = sum(W2)./sum(W);

    clearvars W W2;
end

log3S = log(S)./repmat(log(3),size(S,1),size(S,2));
HS = -sum(S.*log3S,2);  

TS = 2;
DNewS = NewData(HS<=TS,:);
EigS = Eig(HS<=TS,:);
[~,S_Idx] = max(S(HS<=TS,:),[],2);


idx(:,4) = S_Idx;
Struct.K_Smooth_Data = DNewS;
Struct.K_Smooth_Idx = S_Idx;
