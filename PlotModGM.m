function [X,Y,Z,P] = PlotModGM(DataTot,TotMod,factor)

Int = max(DataTot) - min(DataTot);
bin_x = max(Int)/70;
x = min(DataTot(:,1)) - 5*bin_x:bin_x:max(DataTot(:,1)) + 5*bin_x;
y = min(DataTot(:,2)) - 5*bin_x:bin_x:max(DataTot(:,2)) + 5*bin_x;
z = min(DataTot(:,3)) - 5*bin_x:bin_x:max(DataTot(:,3)) + 5*bin_x;
[X,Y,Z] = meshgrid(x,y,z);
P = zeros(size(X));
for i=1:size(X,1)
% tic;
disp(i)
for j=1:size(X,2)
%         for k=1:size(X,3)
Data_Post = reshape([X(i,j,:) Y(i,j,:) Z(i,j,:)],3,size(X,3))';
P(i,j,:) = pdf(TotMod,Data_Post);
%         end
end
% toc;
end
% figure;
isovalue = max(max(max(P)))*factor;
IsoSurf = isosurface(X,Y,Z,P,isovalue);
p2 = patch(isosurface(X,Y,Z,P,isovalue));
isonormals(X,Y,Z,P,p2);
set(p2,'FaceColor',[1 0 0],'EdgeColor',[0.1 0.1 0.1]);
p2.FaceAlpha = 0.1 ;
p2.EdgeAlpha = 0.1;
drawnow