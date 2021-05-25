function rbfnet = rbf(G,LatentD,c_width,actfcn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% RBF sets teh RBF field in AGTM.                                         %
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
%
%	Description
%   This functions sets the Radial Basis Function net as a field of AGTM.
%   It first callas GRAPH_EMESH to identify the centers of the RBFs, then
%   it computes the geodesic distances of all nodes on abstract graph G
%   w.r.t. these centers and computes the activation functions.
%

rbfnet.type = 'RBFnet';
rbfnet.actfcn = actfcn;
rbfnet.LatentD = LatentD;
rbfnet.graph = G;

% Build the epsilon-mesh on the graph
RBF_c = Graph_EMesh(G,c_width,LatentD);
rbfnet.c = RBF_c;

% Recover the size of the basis function centers (geodesic distances 
% nbetween nodes of the graph)

cdist = distances(G,RBF_c,RBF_c);
cdist = cdist + 1e20.*eye(size(RBF_c,1));
c_width = min(cdist(:));


rbfnet.c_width = c_width;

V = 1:size(G.Nodes,1);
rbfnet.V = V';

Dist = distances(G,RBF_c,V');

factor = 3;%size(G.Nodes,1)./size(RBF_c,1);

Dist2 = Dist.^2;

% computation of Phi
c_width2 = (factor*c_width).^2;
rbfnet.Phi = exp(-Dist2./(2*c_width2));