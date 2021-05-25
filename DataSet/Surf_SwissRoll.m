function Data = Surf_SwissRoll(N_Data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% SURF_SWISSROLL builds a Swiss roll sampled from crisp manifold.         %
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

%SURF_SWISSROLL builds a Swiss roll sampled from crisp manifold.


    N = N_Data;
    %noise = 0.05;
    t = 3*pi/2 * (1 + 2*rand(N,1));
    h = 11 * rand(N,1);
    Data = [t.*cos(t), h, t.*sin(t)]; %+ noise*randn(N,3);
    Data = Data./10;

    clearvars N noise t h
end