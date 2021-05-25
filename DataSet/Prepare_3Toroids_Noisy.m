function [D,DN,DNN] = Prepare_3Toroids_Noisy(N_Data)

[D,L] = Surf_3Toroids(N_Data,false);
D = [D L];
figure; plot3(D(:,1),D(:,2),D(:,3),'.','markersize',3)
axis equal
ax = gca;
xlim = ax.XLim;
ylim = ax.YLim;
zlim = ax.ZLim;

close

Lim = [xlim(1) - 0.2 xlim(2) + 0.2; ylim(1) - 0.2 ylim(2) + 0.2; zlim(1) - 0.2 zlim(2) + 0.2];
Len = Lim(:,2) - Lim(:,1);

% UNoise = Lim(:,1)' + Len'.*rand(100000,3);

[DN,Labels] = Surf_3Toroids(N_Data,true,0.1);

DN = [DN Labels];
UNoise2 = Lim(:,1)' + Len'.*rand(N_Data,3);
IdTol = ismembertol(UNoise2(:,1:3),DN(:,1:3),1e-2,'ByRows',true);
NN = UNoise2(~IdTol,:);

DNN = [DN; [NN zeros(size(NN,1),1)]];
% LNN = [Labels; zeros(size(NN,1),1)]; 
IdPerm = randperm(size(DNN,1));
DNN = DNN(IdPerm,:);

