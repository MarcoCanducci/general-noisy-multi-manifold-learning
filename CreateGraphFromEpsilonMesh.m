function [SubGraph,SNodes] = CreateGraphFromEpsilonMesh(G1,r)

    % Define radius r for epsilon mesh
    Nodes = str2num(char(G1.Nodes.Name));

    DG = G1.Edges.Distance;
    N1 = str2num(char(G1.Edges.EndNodes(:,1)));
    N2 = str2num(char(G1.Edges.EndNodes(:,2)));
    s = zeros(size(N1,1),1);
    t = s;
    for i=1:size(N1,1)
        s(i) = findnode(G1,num2str(N1(i,:)));
        t(i) = findnode(G1,num2str(N2(i,:)));
    end
    size(s)
    G = graph(s,t,DG);

    [centres,MInd] = Graph_EMesh(G,r,2);

    Dist = distances(G,centres,centres);
    AdjS = Dist;
    AdjS(AdjS>2*r) = 0;
    if issymmetric(AdjS) ==0
        AdjS = 0.5.*(AdjS + AdjS');
    end
%     figure; imagesc(AdjS)
%     figure; gplot3(AdjS,Nodes(centres,:)); axis equal
    SNodes = Nodes(centres,:);
    SNodesC = num2cell(num2str(Nodes(centres,:)),2);
    SubGraph = graph(AdjS,SNodesC);