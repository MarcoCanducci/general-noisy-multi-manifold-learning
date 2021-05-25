function net = agtm(Graph,LatentD,c_width)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% AGTM builds the struct array net, containing all AGTM fields            %
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

%	Description
%   This functions sets the AGTM net struct array.
%   The net array contains the RBFs net and GMMNET net plus information on
%   the global setup of AGTM:
%   - Initialization Graph;
%   - Latent dimensionality;
%   - Type of activation function;
%

net.type = 'AGTM';
net.LatentD = LatentD;
net.graph = Graph;

%setup rbfnet

actfcn = 'Gaussian';
Adj = adjacency(Graph);
LGraph = graph(Adj);
% LGraph.Edges.Weight = Graph.Edges.Distance;

net.rbfnet = rbf(LGraph,LatentD,c_width,actfcn);

net.gmmnet = gmm(Graph,LatentD);

net = rbffwd(net);