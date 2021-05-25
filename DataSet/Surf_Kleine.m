function Data = Surf_Kleine(N_Data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% SURF_KLEINE builds a Klein bottle sampled from crisp manifold.          %
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
    a = 3;
    b = 4;
    c = 2;
    
    u1 = pi.*rand(N_Data,1);
    u2 = pi + pi.*rand(N_Data/2,1);
    v1 = 2*pi.*rand(N_Data,1);
    v2 = 2*pi.*rand(N_Data/2,1);
    
    r1 = c.*(1 - cos(u1)./2);
    r2 = c.*(1 - cos(u2)./2);
    
    x1 = (a.*(1+sin(u1)) + r1.*cos(v1)).*cos(u1);
    y1 = (b + r1.*cos(v1)).*sin(u1);
    z1 = r1.*sin(v1);
    x2 = a.*(1 + sin(u2)).*cos(u2) - r2.*cos(v2);
    y2 = b.*sin(u2);
    z2 = r2.*sin(v2);
    X = [x1; x2];
    Y = [y1; y2];
    Z = [z1; z2];
    
    XN = 0.5.*(X-min(X))/(max(X)-min(X));
    YN = 0.5.*(Y-min(Y))/(max(Y)-min(Y));
    ZN = 0.5.*(Z-min(Z))/(max(Z)-min(Z));
    
    Data = [XN YN ZN];

    
end