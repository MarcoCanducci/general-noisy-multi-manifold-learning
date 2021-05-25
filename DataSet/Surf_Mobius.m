function Data = Surf_Mobius(N_Data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% SURF_MOBIUS builds a Mobius strip sampled from crisp manifold.          %
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

u = linspace(0,2*pi,500);
    v = linspace(-0.5,0.5,500);
    [u,v] = meshgrid(u,v);
    x = (1+v.*cos(u/2)).*cos(u);
    y = (1+v.*cos(u/2)).*sin(u);
    z = v.*sin(u/2);

    M = zeros(size(u,2)^2,3);
    L1 = length(x(:,1,1));
    L2 = length(y(1,:,1));
    L3 = length(z(1,1,:));
    T = 1;
    for i = 1:1:L1
        for j = 1:1:L2
            for k = 1:1:L3
            M(T,1) = x(i,j,k) ;
            M(T,2) = y(i,j,k) ;
            M(T,3) = z(i,j,k) ;
            T = T+1;
            end
        end
    end

    Data = datasample(M,N_Data);
    clearvars M T L1 L2 L3 u v x y z i j k
end