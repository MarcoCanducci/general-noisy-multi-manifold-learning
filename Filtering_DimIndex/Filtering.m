function [Labels, Neighbors1, Neighbors2] = Filtering(Data_Original,Data,Radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Filtering, given a noisy data set and a diffused one (obtained by SAF)  %
%filters out noisy data points.                                           %
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

Neighbors1 = rangesearch(Data_Original,Data_Original,Radius);
Neighbors2 = rangesearch(Data,Data,Radius);
Labels = zeros(size(Data,1),1);
for i=1:size(Data,1)
    % By evaluating neighborhoods in the original dataset we can discard
    % all points which didn't have neighbors populated enough to survive
    % before appyling SAF. This step removes artifacts of the SAF
    % algorithm.
    if(Labels(i) == 0)
        if(size(Neighbors2{i},2)<=5 || size(Neighbors1{i},2)<=5)
%         cond =abs(size(Neighbors2{i},2) - size(Neighbors1{i},2))/size(Neighbors1{i},2);
% %         cond2 = 
%         if(cond<=0.5)    
            continue;
        else        
            Labels(Neighbors2{i}) = 1;
        end
    end     
end