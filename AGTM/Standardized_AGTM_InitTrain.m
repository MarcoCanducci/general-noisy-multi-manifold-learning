function [net,logL,GMDist,NoisyMan] = Standardized_AGTM_InitTrain(FG,F_NoisyData,IntDim,r,epsilon,mem)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Standardized_AGTM_InitTrain sets up AGTM and performs training calling  %
% AGTM_EM after standardization of datasets and graphs                    %
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

%figure; hold on
F_NoisyData = double(F_NoisyData);
k =1;
for i=1:size(FG,2)  
    %tic;
    %disp(i)
    disp(['net n. ', int2str(i)]);
    Adj = adjacency(FG{i});    
    Nodes = str2num(char(FG{i}.Nodes.Name));
%     NN = rangesearch(F_NoisyData{i},Nodes,r);
    NN = rangesearch(F_NoisyData,Nodes,r);
%     size(NN)
    T = F_NoisyData(unique(cell2mat(NN')),:);
    %size(T,1);
%     T = F_NoisyData{i};%(NN{1},:);
%     T = [];
%     for j = 1:size(NN,1)
%         T = [T; F_NoisyData(NN{j},:)];
%     end
%     T = unique(T,'rows','stable');
%     if(size(T,1) < 2*size(Nodes,1))
%         disp('Not enough training points')
%         continue
%     end
%     plot3(T(:,1),T(:,2),T(:,3),'b.','markersize',3)

%     [Zscore,Mu,S] = zscore(T);
%     Nodes_Z = (Nodes - Mu)./S;
%     r_AGTM = max(min(pdist2(Nodes_Z,Nodes_Z) + 1e10.*eye(size(Nodes_Z,1))));
    G = graph(Adj);
%     G.Nodes.Name = num2cell(num2str(Nodes_Z),2);
    G.Nodes.Name = num2cell(num2str(Nodes),2);
    
%     if (size(Nodes_Z,1)<=2*epsilon)
%         epsilon = 1;
%     end
%     net{k} = gtminit(G,IntDim,epsilon,Zscore,r_AGTM);
%     [net{k},logL{k}] = agtm_em(net{k},30,Zscore,mem);
    r_AGTM = mean(min(pdist2(Nodes,Nodes)+1e10.*eye(size(Nodes,1))));
    net{k} = gtminit(G,IntDim,epsilon,T,r);
    [net{k},logL{k}] = agtm_em(net{k},20,T,mem);    
%     SigmaNew = zeros(size(net{k}.gmmnet.Sigma,1),size(net{k}.gmmnet.Sigma,2),size(net{k}.gmmnet.Sigma,3));
%     for j=1:size(net{k}.gmmnet.Sigma,3)
%         SigmaNew(:,:,j) = diag(S)*net{k}.gmmnet.Sigma(:,:,j)*diag(S);
%         if(issymmetric(SigmaNew(:,:,j)) == 0)
%             SigmaNew(:,:,j) = (SigmaNew(:,:,j) + SigmaNew(:,:,j)').*0.5;
%         end
%     end
%     Nodes = net{k}.gmmnet.V_.*S + Mu;
%     net{k}.gmmnet.V_ = Nodes;
%     net{k}.gmmnet.Sigma = SigmaNew;
% %     SigmaNew = net{k}.gmmnet.Sigma;
% %     Nodes = net{k}.gmmnet.V_;
% %     gplot3(Adj,Nodes);
% %     axis equal
% %     drawnow
    priors = net{k}.gmmnet.GMdist.ComponentProportion;
%     GMDist{k} = gmdistribution(Nodes,SigmaNew,priors);
    GMDist{k} = gmdistribution(net{k}.gmmnet.V_,net{k}.gmmnet.Sigma,priors);
    NoisyMan{k} = T;
    k = k+1;
    
   % fprintf('Time for agtm model n. %d \n',k-1)
    %toc; 

end