function [G,SC] = NeighbourhoodGraph(D,r,ldim)

% r = 0.05;
% D = D1;
% ldim = 1;
% eps = 0.1*r
r = 0.5*r;
NN = rangesearch(D,D,r);
CN = {};
k =1;
for i=1:size(NN,1)
    if(size(NN{i},2)>=5)
        CN{k} = NN{i};
        k = k+1;
    end
end

F_D = D(unique(cell2mat(CN)),:);
Dim = size(F_D,2);
epsilon = 1e-2;
[SP,NP] = DatasetPartition(F_D,r);

j = 1;
k = 1;
Id_Remove = [];
UC = {};
for i=1:size(NP,2)
    U = pca(NP{i});
    if(rank(U)<ldim)
        Id_Remove(k) = i;
        k = k + 1;
    else
        UC{j} = U;
        j =j+1;
    end
end

Keep_Id = setdiff(1:size(SP,1),Id_Remove);
SC = SP(Keep_Id,:);

Adj_r = pdist2(SC,SC);
Adj_r(Adj_r > 2*r + eps) = 0;
Adj_angle = zeros(size(Adj_r));
Ind = find(Adj_r);
[Row,Col] = ind2sub(size(Adj_r),Ind);

for i=1:size(Ind,1)
    Adj_angle(Row(i),Col(i)) = subspace(UC{Col(i)}(:,1:ldim),UC{Row(i)}(:,1:ldim));
    Adj_angle(Col(i),Row(i)) = Adj_angle(Row(i),Col(i));
end
% for i=1:size(Adj_r,1)
%     if(nnz(Adj_angle(i,:))>2*ldim)
%         [angle_sorted, sort_order] = sort(Adj_angle(i,:),'descend');
%         Adj_angle(i,sort_order(2*ldim+1:end)) = 0;
%         Adj_angle(sort_order(2*ldim+1:end),i) = 0;
% %         Adj_angle(Row(i),sort_order(2*ldim+1:end)) = 0;
%     end
% %     
% end
% Adj_angle = Adj_angle*180/pi;
T = 0.8;
Adj_angle(Adj_angle > T*pi/2) = 0;
% TAdj_angle = Adj_angle';
% Adj_angle(Adj_angle~=TAdj_angle) = 0;
Adj = Adj_r;
Adj(Adj_angle==0) = 0;
Adj(Adj_angle~=0) = r;

% Adj(Adj_angle*180/pi > 70) = 0;

G = graph(Adj,'upper');
G.Nodes.Name = num2cell(num2str(SC),2);
G.Nodes.TSpace = UC';

CC = G.conncomp;
unique(CC)