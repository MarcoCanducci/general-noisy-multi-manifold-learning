function net = gtminit(Graph,LatentD,c_width,Data,r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% GTMINIT initializes the AGTM net using the graph, latent dimension, size%
% of RBFs kernels, noisy data set and neighbourhood search radius inputs. %
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
%	Description
%   This functions firstly sets up the AGTM net. It then procedes to
%   compute a the weighted covariance matrix of a small neighbourhood
%   centered at each node of the graph.

% [Z,Mu,S] = zscore(Data);

net = agtm(Graph,LatentD,c_width);
% net.gmmnet.V_ = str2num(char(Graph.Nodes.Name));
V_ = net.gmmnet.V_;

[m,n] = size(V_);

KCov_width = 2*r;

N = rangesearch(Data,V_,KCov_width);
Sigma1 = ones(n,n,m);
Sigma = Sigma1;

for i=1:size(N,1)
    Sigma1(:,:,i) = KCov(V_(i,:),Data(N{i},:),KCov_width);
    Sigma1(:,:,i) = 0.5.*(Sigma1(:,:,i) + Sigma1(:,:,i)') + 1e-5.*eye(n);
    Sigma(:,:,i) = cov(Data(N{i},:));
    Sigma(:,:,i) = 0.5.*(Sigma(:,:,i) + Sigma(:,:,i)');
    if(any(any(isnan(Sigma1(:,:,i)))) || cond(Sigma1(:,:,i))>1e3)
        Sigma1(:,:,i) = eye(3).*1e-2;
    end
end

net.gmmnet.T = Data;
net.gmmnet.Sigma1 = Sigma1;
net.gmmnet.Sigma = Sigma;
net.gmmnet.KCov_width = KCov_width;
net.gmmnet.priors = ones(size(V_,1),1).*(1/size(V_,1));
net.gmmnet.alpha = 1;