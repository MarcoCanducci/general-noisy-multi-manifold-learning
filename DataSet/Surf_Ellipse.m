function Data = Surf_Ellipse(N_Data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%SURF_ELLIPSE builds an elliptical cup surface sampled from crisp manifold%
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

    a = 2;
    b = 1;
    z = -1 + 2.*rand(N_Data,1);
    t = 2*pi.*rand(N_Data,1);
    x = sqrt(1-z.^2).*cos(t).*a^2;
    y = sqrt(1-z.^2).*sin(t).*b^2;
    Data = [x,y,z];
    clearvars x y z t a b c

end