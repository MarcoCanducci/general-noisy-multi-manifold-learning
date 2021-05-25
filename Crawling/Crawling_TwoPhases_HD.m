function [G,R,c] = Crawling_TwoPhases_HD(D,r,ldim,betha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% CRAWLING_TWOPHASES is a recursive procedure for building graphs locally %
% aligned with the local tangent space of a manifold.                     %
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
%   Description:
%   This function, given a manifold sampled in point cloud, recovers an 
%   abstract and an embedded graph lying on the local tangent spaces of the
%   manifold and describing its geometry.
%   It is compose of three phases:
%   - Initialiation: selects randomly a seed and computes the local tangent
%   space of the manifold on that seed. It then identifies principal
%   directions by performing decomposition of the local covariance matrix
%   and estimates new candidate nodes.
%   - Expansion: given the set of nodes of degree 1 (having only 1 incident
%   edge), it expands the crawling on the new locally linear principal
%   directions, looking for new candidate nodes of the graph.
%   - Contraction: Aims at contracting the set of new nodes by collapsing
%   neighbouring nodes into one single node, given an accuracy threshold.
%
%   Output:
%   - Residual set R;
%   - Embedded graph G (from which the abstract graph is extracted);
%
    % Initialization
    NEdges = ldim*2;
    alpha = 0.75;
    %betha = 0.4;
    R = D;
    G = graph;
    c = 0;
    %disp(size(R))
%     M = rangesearch(R,R,r);
%     for i = 1:size(M,1)        
%         if(size(M{i},2)<5)
%             R = setdiff(R,D(M{i},:),'rows','stable');
%         end
%     end
    %disp(size(R))
    
    
    %Random initialization
    t0 = datasample(R,1);
%    t0 = R(1,:);
    N = rangesearch(R,t0,r);
    
    G = addnode(G,num2str(t0));
    Point_ID = 1;
    
    % First iteration (no eigenvector projection)

    DN = R(N{1},:);
    R = setdiff(R,t0,'rows','stable');
    U = pca(DN);
    if(size(U,2)<ldim)
        return
    end
    tNp = repmat(t0,ldim,1) + (alpha*r).*U(:,1:ldim)';
    tNn = repmat(t0,ldim,1) - (alpha*r).*U(:,1:ldim)';

    tN = [tNp; tNn];
    tNN_Idx = knnsearch(R,tN,'k',1);
    tNN = R(tNN_Idx,:);
    PDir = tNN - repmat(t0,NEdges,1);
    
    Mod = vecnorm(PDir');
    PDir = PDir./Mod';
    
    DotPDirMat = abs(PDir*PDir');
    m = size(PDir,1);
    tNN = unique(tNN,'rows','stable');
    for i=1:size(tNN,1)
        Nodes = str2num(char(G.Nodes.Name));
        [ISM,ISM_ID] = ismembertol(tNN(i,:),Nodes,1e-4,'ByRows',true);
        if(ISM == 1)
                x = Nodes(ISM_ID,:);
                G = CreateEdge(G,x,ISM_ID,PDir,m,i,1,t0);
        else
                G = addnode(G,num2str(tNN(i,:)));
                R = setdiff(R,tNN(i,:),'rows','stable');
                Point_ID = Point_ID + 1;
                
                [~,ord] = sort(DotPDirMat(i,:),'ascend');
                SubSp(1,:) = PDir(i,:);
                SubSp(2,:) = PDir(ord(1),:);
                if ldim>2
                    for k=2:size(ord,2)
                        if DotPDirMat(ord(k),ord(1)) < 0.2 && ...
                                DotPDirMat(ord(k),i) < 0.2
                            SubSp(3,:) = PDir(ord(k),:);
                            break
                        else
                            SubSp(3,:) = PDir(ord(2),:);
                        end
                    end
                end                
                
                EdgeTable = table([1,i],{SubSp},Mod(i),...
                     'VariableNames',{'EndNodes','Subspace_Basis','Distance'});
                G = addedge(G,EdgeTable);
        end
        Nodes = str2num(char(G.Nodes.Name));    
    end

    d = degree(G);

    for i=1:NEdges
        NewGen{i} = str2num(char(G.Nodes{d==i,:}));
    end
    Nodes = str2num(char(G.Nodes.Name));
    A = size(Nodes,1); 
    B = A - 1;
    it = 1;

    % Alternating Expansion and contraction phases    
    while(isempty(NewGen{1})==0 && size(R,1)>20 && B < A)
        
        % Expansion Phase
        B = A;
        for i=1:size(NewGen{1},1)
            t0 = NewGen{1}(i,:);

            N = rangesearch(R,t0,r);
            if(size(N{1},2) <= 5)
                % border found
                R = setdiff(R,R(N{1},:),'rows','stable');
                continue 
            end

            DN = R(N{1},:);
            V = pca(DN);

            % Projection of parent subspace onto new eigenvectors
            
            NodeID = findnode(G,num2str(t0));
            GN = neighbors(G,NodeID);
            GE = findedge(G,NodeID, GN);

            S = G.Edges.Subspace_Basis(GE);
            Sub = S{1}';

            P = V(:,1:ldim)*pinv(V(:,1:ldim));
            U = P*Sub(:,1:ldim);
            for j=1:ldim
                U(:,j) = U(:,j)./norm(U(:,j));
            end
            
            tNp = repmat(t0,ldim,1) + (alpha*r).*U(:,1:ldim)';
            tNn = repmat(t0,ldim,1) - (alpha*r).*U(:,1:ldim)';
            
            tN = [tNp; tNn];
            tNN_Idx = knnsearch(R,tN,'k',1);
            tNN = R(tNN_Idx,:);

            PDir = tNN - repmat(t0,NEdges,1);
            Mod = sqrt(sum(PDir.^2,2));

            PDir = PDir./repmat(Mod,1,size(PDir,2));
            m = size(PDir,1);        

            % Contraction phase

            for j = 1:size(tN,1)
                tNNeigh = rangesearch(R,tNN(j,:),r);
                if(size(tNNeigh{1},2)<5)
                    continue
                end
                [ISM,ISM_ID] = ismembertol(tNN(j,:),Nodes,1e-4,'ByRows',true);
                if ISM == 1
                    NNN = neighbors(G,ISM_ID);
                    if(size(NNN,1) < NEdges)
                        x = Nodes(ISM_ID,:);
                        x_ID = ISM_ID;
                        G = CreateEdge(G,x,x_ID,PDir,m,j,NodeID,t0);
                    end
                else
                    RS = rangesearch(Nodes,tN(j,:),betha*r);
                    if(size(RS{1},2)~=0 )
                        NNN = neighbors(G,RS{1}(1));
                        if(size(NNN,1) < NEdges)
                            x = Nodes(RS{1}(1),:);
                            x_ID = RS{1}(1);
                            G = CreateEdge(G,x,x_ID,PDir,m,j,NodeID,t0);
                        end
                    else
                        % Create new node and connect to t0;
                        NNN = neighbors(G,NodeID);
                        if(size(NNN,1) < NEdges)
                            G = addnode(G,num2str(tNN(j,:)));
                            R = setdiff(R,tNN(j,:),'rows','stable');
                            Point_ID = Point_ID + 1;
                            
                            DotPDirMat = abs(PDir*PDir');
                            [~,ord] = sort(DotPDirMat(j,:),'ascend');
                            SubSp(1,:) = PDir(j,:);
                            SubSp(2,:) = PDir(ord(1),:);
                            if ldim>2
                                for k=2:size(ord,2)
                                    if DotPDirMat(ord(k),ord(1)) < 0.2 && ...
                                            DotPDirMat(ord(k),j) < 0.2
                                        SubSp(3,:) = PDir(ord(k),:);
                                        break
                                    end
                                end
                            end                                  
                            
                            EdgeTable = table([NodeID,Point_ID],...
                                {SubSp},Mod(j),...
                                         'VariableNames',...
                                         {'EndNodes','Subspace_Basis','Distance'});
                            G = addedge(G,EdgeTable);
                        end
                    end
                end
                Nodes = str2num(char(G.Nodes.Name));
            end
        end
        G = rmedge(G, 1:numnodes(G), 1:numnodes(G));        
        d = degree(G);
        Nodes = str2num(char(G.Nodes.Name));
        
        % Updating the set of 1-edge nodes.
        for i=1:NEdges
            NewGen{i} = str2num(char(G.Nodes{d==i,:}));
        end
        if(~isempty(NewGen{NEdges}))
            R4 = rangesearch(R,NewGen{NEdges},betha*r);
%             R2 = R;
            R = setdiff(R,R(unique(cell2mat(R4')),:),'rows','stable');
%             for j=1:size(NewGen{NEdges},1)
%                 if(isempty(R4{j})==0)
%                 R2 = setdiff(R2,R(R4{j},:),'rows','stable');
%                 end
%             end
%             R = R2; clearvars R2;
        end
        % Updating size of nodes for conditional break.
        A = size(Nodes,1);
        it = it + 1;
    end
    % Removing self-loops from graph G

end

function G = CreateEdge(G,x,x_ID,PDir,m,j,NodeID,t0)
    
    RecE = findedge(G, NodeID, x_ID);
    if(RecE == 0)
        DotPDirMat = abs(PDir*PDir');
        [~,ord] = sort(DotPDirMat(j,:),'ascend');
        SubSp(1,:) = PDir(j,:);
        SubSp(2,:) = PDir(ord(1),:);
        if size(PDir,1)*0.5 > 2
            for k=2:size(ord,2)
                if DotPDirMat(ord(k),ord(1)) < 0.2 && ...
                        DotPDirMat(ord(k),j) < 0.2
                    SubSp(3,:) = PDir(ord(k),:);
                    break
                else
                    SubSp(3,:) = PDir(ord(2),:);
                end
            end
        end       

        %Create edge 
        Dir = x - t0;
        Mod2 = norm(Dir);
%         Dir = Dir./Mod2;
        EdgeTable = table([NodeID,x_ID],{SubSp},Mod2,...
              'VariableNames',{'EndNodes','Subspace_Basis','Distance'});
        G = addedge(G,EdgeTable);
    end 
end  