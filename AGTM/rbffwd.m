function net = rbffwd(net)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% RBFFWD propagates the abstract graph nodes onto the embedding via the   %
% RBF network.                                                            %
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
%   This function propagates the abstract graph nodes onto the emebedding
%   using the RBF network. The new embedded nodes are stired in the field
%   GMMNET.
gmm = net.gmmnet;

rbfnet = net.rbfnet;

V_= gmm.V_;
Phi = rbfnet.Phi;

% if net.LatentD == 1
%     lambda = 1e-5;
% else
%     lambda = 1e-2;
%     
% 
% end
lambda =1e-5;
H = Phi*Phi';
%     A1 = inv(H + 1e-2.*eye(size(H,1)));
W = (H + lambda.*eye(size(H,1)))\Phi*V_;
% W = Phi' \ V_;
net.gmmnet.W = W;
net.gmmnet.V_ = (W'*Phi)';




