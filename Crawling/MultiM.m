function [FG,NoisyMan,TotR] = MultiM(D,NoisyD,r,ldim,betha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% MULTIM is a recursive procedure for trimming graphs recovered through   % 
% crawling and recovering of the noisy sampled manifold.                  %
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
%   Description:
%   This function, given a point cloud sampling multiple noisy manifolds
%   and a data set containing the same manifolds but collapsed onto theri
%   respective spine, recovers a set of graphs and a set of point clouds
%   describing the geometry of the manifolds and the noisy distributions of
%   points surrounding them. It applies Manifold crawling iteratively until
%   the residual set's size is minimized.
%
%   See also: Crawling_TwoPhases

TotR = D;
it = 1;
rN = r;
while(size(TotR,1)>20)
%     if ldim==2
        [G,TotR,c] = Crawling_TwoPhases(TotR,r,ldim,betha);
%     else        
%         [G,TotR,c] = Crawling_TwoPhases(TotR,r,ldim,betha);
%     end
    if c == 0
%         R2 = TotR;
        Nodes = str2num(char(G.Nodes.Name));
        if(size(Nodes,1)>ldim*2+1)

            NN = rangesearch(TotR,Nodes,rN);
            TotR = TotR(setdiff(1:size(TotR,1),unique(cell2mat(NN'))),:);
%             for i=1:size(NN,1)
%                 R2 = setdiff(R2,TotR(NN{i},:),'rows','stable');
%             end
%             TotR = R2;
%             clearvars R2;
            Dist = G.Edges.Distance;
            CollG{it} = G;
%             CollG{it} = rmedge(CollG{it},find(Dist>r));
            it = it + 1;
        end
    end
        
end

clearvars NN;

k = 1;
for i=1:size(CollG,2)
    bin = conncomp(CollG{i});
    if size(unique(bin),2) > 1        
        for j=1:size(unique(bin),2)
            SubGraph = subgraph(CollG{i},find(bin==j));
            Nodes = str2num(char(SubGraph.Nodes.Name));
            if(size(Nodes,1) > ldim*2 + 1)
                NN = rangesearch(NoisyD,Nodes,rN);
%                 for l=1:size(NN,2)
%                     if(size(NN{l},2)<=5)
%                         break
%                     end
%                 end
                NoisyMan{k} = NoisyD(unique(cell2mat(NN')),:);
%                 for l = 2:size(NN,1)
%                     NoisyMan{k} = union(NoisyMan{k},...
%                         NoisyD(NN{l},:),'rows','stable');
%                 end
                FG{k} = SubGraph;
                k = k + 1;
                clearvars NN
            end
                
        end
    else
        SubGraph = CollG{i};
        Nodes = str2num(char(SubGraph.Nodes.Name));
        if(size(Nodes,1) > ldim*2 + 1)
            NN = rangesearch(NoisyD,Nodes,rN);
            NoisyMan{k} = NoisyD(unique(cell2mat(NN')),:);
%             NoisyMan{k} = NoisyD(NN{1},:);
%             for l = 2:size(NN,1)
%                 NoisyMan{k} = union(NoisyMan{k},...
%                     NoisyD(NN{l},:),'rows','stable');
%             end
            clearvars NN;
            FG{k} = SubGraph;
            k = k + 1;
        end
    end
end