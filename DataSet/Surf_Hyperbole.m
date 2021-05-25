function Data = Surf_Hyperbole(N_Data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% SURF_HYPERBOLE builds an hyperbolic line sampled from crisp manifold.   %
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
    N_points = N_Data;
    ZB = linspace(0.4,0.6,N_points)';
    R = linspace(0,0.5,N_points)';
    RadiusB = -R.^2 + 0.25;
    tB = 2*pi.*rand(N_points,1);
    XB = RadiusB.*sin(tB) + 0.5;
    YB = RadiusB.*cos(tB) -0.5 ;

    Data = [XB YB ZB];

    clearvars ZB R RadiusB tB xB YB N_points

end