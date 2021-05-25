function [S,N] = DatasetPartition(D,r)


Res = D;
T = [];
x = datasample(Res,1);
S(1,:) = x;
Nr = rangesearch(Res,x,r);
N{1} = Res(Nr{1},:);
T = [T; N{1}];
i = 2;

while(size(Res,1) > 5)
    x = datasample(Res,1);
%    S(i,:) = x;
    Nr = rangesearch(Res,x,r);
    N{i} = Res(Nr{1},:);
    S(i,:) = mean(N{i});
    T = [T; N{i}];
    Res = setdiff(D,T,'rows','stable');
    i = i +1;

end
    
    