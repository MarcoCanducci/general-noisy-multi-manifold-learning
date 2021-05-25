function [G,R,c] = Crawling_Partition(SP,NP,D,r,ldim,betha)
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
    %betha = 0.4;
    R = D;
    G = graph;
%     c = 0;
%     disp(size(R))
%     M = rangesearch(R,R,r);
%     for i = 1:size(M,1)        
%         if(size(M{i},2)<5)
%             R = setdiff(R,D(M{i},:),'rows','stable');
%         end
%     end
    
%     [SP,NP] = DatasetPartition(R,r*0.5);
    
    t0 = SP(1,:);
    DN = NP{1};
    
    G = addnode(G,num2str(t0));
    Point_ID = 1;
    
    % First iteration (no eigenvector projection)
    
    SP(1,:) = [];
    NP(1) = [];
    U = pca(DN);
    if(size(U,2)<ldim)
        return
    end
    
    [kN,kNDist] = knnsearch(SP,t0,'k',2*ldim);
    
    DirV = (SP(kN,:) - t0)./kNDist';
    
    Proj_Coeff = U(:,1:ldim)'*DirV';
    k = Point_ID + 1;

    for j=1:size(kN,2)
        Proj(j,:) = sum(Proj_Coeff(:,j).*U(:,1:ldim)',1);
        angle(j) = subspace(Proj(j,:)',DirV(j,:)')*180/pi;
        if (kNDist(j)<=r && angle(j) < 90)
            
            G = addnode(G,num2str(SP(kN(j),:)));
            EdgeTable = table([Point_ID,k],'VariableNames',{'EndNodes'});
            G = addedge(G,EdgeTable);
            SP(kN(j,:)) = []
            k = k + 1;
        end
    end
    
    d = degree(G);

    for i=1:NEdges
        NewGen{i} = str2num(char(G.Nodes{d==i,:}));
    end
    
    
    P = U(:,1:ldim)*pinv(U(:,1:ldim));
    U = P*Sub(:,1:ldim);
    
    tNp = repmat(t0,ldim,1) + r.*U(:,1:ldim)';
    tNn = repmat(t0,ldim,1) - r.*U(:,1:ldim)';

    tN = [tNp; tNn];
    tNN_Idx = knnsearch(SP,tN,'k',1);
    tNN = R(tNN_Idx,:);
    PDir = tNN - repmat(t0,NEdges,1);
    Mod = sqrt(sum(PDir.^2,2));

    PDir = PDir./repmat(Mod,1,size(PDir,2));
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
                EdgeTable = table([1,i],{[PDir(i,:); PDir(mod(i,m) + 1,:)]},Mod(i),...
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

                            EdgeTable = table([NodeID,Point_ID],...
                                {[PDir(j,:); PDir(mod(j,m) + 1,:)]},Mod(j),...
                                         'VariableNames',...
                                         {'EndNodes','Subspace_Basis','Distance'});
                            G = addedge(G,EdgeTable);
                        end
                    end
                end
                Nodes = str2num(char(G.Nodes.Name));
            end
        end
        d = degree(G);
        Nodes = str2num(char(G.Nodes.Name));
        
        % Updating the set of 1-edge nodes.
        for i=1:4
            NewGen{i} = str2num(char(G.Nodes{d==i,:}));
        end
        if(isempty(NewGen{NEdges})==0)
            R4 = rangesearch(R,NewGen{NEdges},betha*r);
            R2 = R;
            for j=1:size(NewGen{NEdges},1)
                if(isempty(R4{j})==0)
                R2 = setdiff(R2,R(R4{j},:),'rows','stable');
                end
            end
            R = R2; clearvars R2;
        end
        % Updating size of nodes for conditional break.
        A = size(Nodes,1);
        it = it + 1;
    end
    % Removing self-loops from graph G
    G = rmedge(G, 1:numnodes(G), 1:numnodes(G));
end

function G = CreateEdge(G,x,x_ID,PDir,m,j,NodeID,t0)
    
    RecE = findedge(G, NodeID, x_ID);
    if(RecE == 0)
        %Create edge 
        Dir = x - t0;
        Mod2 = norm(Dir);
        Dir = Dir./Mod2;
        EdgeTable = table([NodeID,x_ID],{[Dir; PDir(mod(j,m) + 1,:)]},Mod2,...
              'VariableNames',{'EndNodes','Subspace_Basis','Distance'});
        G = addedge(G,EdgeTable);
    end 
end  