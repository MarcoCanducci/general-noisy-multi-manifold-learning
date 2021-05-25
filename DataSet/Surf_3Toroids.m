function [Data,Labels] = Surf_3Toroids(N_Data,noisy,sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% SURF_TOROID builds a toroidal surface sampled from crisp manifold.      %
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


    N_Data = 2*N_Data;
    R = 0.7;
    r = 0.15;

    U = rand(N_Data,1);
    V = rand(N_Data,1);
    W = rand(N_Data,1);

    theta = 2*pi.*U;
    phi = 2*pi.*V;

    Const = (R + r.*cos(theta))./(R + r);

    j = 1;
    for i=1:size(W,1)
        if(W(i) <= Const(i))
            x{j} = (R + r*cos(theta(i)))*cos(phi(i));
            y{j} = (R + r*cos(theta(i)))*sin(phi(i));
            z{j} = r*sin(theta(i));
            j = j+1;
        end
    end
    
    if noisy
        x1 = cell2mat(x')  + sigma.*rand(size(x,2),1);
        y1 = cell2mat(y')  + sigma.*rand(size(x,2),1);
        z1 = cell2mat(z')  + sigma.*rand(size(x,2),1);
    else
        x1 = cell2mat(x');
        y1 = cell2mat(y');
        z1 = cell2mat(z');
    end
    
    T1 = [x1 y1 z1] - [0.95 0 0];
    
    T2 = [x1 y1 z1]*[1 0 0; 0 0 -1; 0 1 0];
    T2 = T2.*[1.4 1.4 1.4] + [0 0 0.12];
    T3 = T1 + [1.9 0 0];
    
    Data = [T1; T2; T3];
    Labels = [ones(size(T1,1),1); 2.*ones(size(T2,1),1); 3.*ones(size(T3,1),1)];
    
    clearvars r R x y z Const j phi theta x1 y1 z1 U V W i
end