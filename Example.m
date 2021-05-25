%This is an example on how to apply filtering and dimensionality index on a
%pair of noisy and diffused (obtained e.g. through SAF algorithm) data sets.

% DTot: noisy data set;
% SAF_D: diffused data set (otained through SAF)
% r: radius for neighbourhood search and dimensionality index. I've found
% good results when the same radius is used for both SAF and Dimensionality
% index, however it can be freely set. Watch out for RAM overload, if it 
% gets to use the swap it will take a long time to run, although sometimes
% it is unavoidable.

% First apply filtering using both the diffused and noisy data sets (SAF_D 
% and DTot respectively). The output is a boolean the size of the data set.
% A label 1 means that point has survived the filtering.

Labels = Filtering(DTot,SAF_D,r); %Set the radius r

% Construct the filtered diffused and noisy data sets
F_Data = SAF_D(Labels==1,:);
DN = DTot(Labels==1,:);

% For each point evaluate its dimensionality index. The output is a struct
% array containing different versions of the index (Read function's header).
% Option 'Barycentric' computes the index with the barycentric
% cohordinates of eigenvalues triades in the total simplex (check the paper)
% Option 'l2' sets a Gaussian Kernel on the vertices of the total simplex
% for the soft index.
Struct = Dim_Index(F_Data,r,'Barycentric','l2');

% The field "K_Smooth_Idx" contains the soft dime index spatially smoothed
% for considering manifolds' edges. You can play around with the other
% fields as well, as well as with the parameters.
Idx = Struct.K_Smooth_Idx;

% Create partitions of both the diffused and noisy filtered data sets.
D1 = F_Data(Idx==1,:);
D2 = F_Data(Idx==2,:);
D3 = F_Data(Idx==3,:);

DN1 = DN(Idx==1,:);
DN2 = DN(Idx==2,:);
DN3 = DN(Idx==3,:);