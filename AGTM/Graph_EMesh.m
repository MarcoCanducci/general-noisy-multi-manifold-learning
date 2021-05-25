function [centres,MInd] = Graph_EMesh(G,r,ldim) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Graph_EMesh constructs an epsilon-mesh of size r on graph G             %
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
%     NEdges = ldim*2;
%     d = degree(G);
%     N4 = find(d>NEdges - 1);
%     i = 1;
%     while(isempty(N4)==0)
%         centres(i,:) = N4(1);
%         MInd{i} = nearest(G,N4(1),r);
%         M4NodeInt = intersect(N4,MInd{i});
%         N4 = setdiff(N4,N4(1));
%         N4 = setdiff(N4,M4NodeInt);
%         i = i + 1;
%     end
% G2 = graph(adjacency(G));
if r==0
    centres = 1:size(G.Nodes,1);
    centres = centres';
else
    
    deg = degree(G);
    NId = find(deg>ldim);
    Id1 = setdiff(1:size(G.Nodes,1),NId);
%     G = G.subgraph(NId);
    N = 1:size(G.Nodes);
    i = 1;
%     for j=1:size(Id,2)
%         centres(j,:) = Id1(j);
%         MInd{j} = nearest(G,Id1(j),r);
%         N = setdiff(N,Id1(j));
%     end
%     i = size(Id,2)+1;
    while(isempty(N) == 0)
        if ismember(N(1),NId)
            centres(i,:) = N(1);
            MInd{i} = nearest(G,N(1),r);
            MNodeInt = intersect(N,MInd{i});
            N = setdiff(N,N(1));
            N = setdiff(N,MNodeInt);
            i = i + 1;
        else
            N = setdiff(N,N(1));
        end
    end
end
%disp(size(centres));
