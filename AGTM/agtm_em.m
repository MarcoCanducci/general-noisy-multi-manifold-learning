function [net,logL] = agtm_em(net,NIter,T,mem)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% AGTM_EM performs E-M algorithm on the struct array net using the data   %
% in T as training data set.                                              %
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
%	Description
%   This function computes the parameters updates for AGTM at each
%   iteration. Parameters W and alpha are updated by computing the
%   responsibilities of net.gmmnet w.r.t. to data set T. This is performed
%   taking advantage of the GMDISTRIBUTION function in MATLAB.
%   At each iteration the new centers and covariance matrices are
%   substitued to the provious ones. The outputs are the updated AGTM and
%   the logL function. logL can be used to inspect convergence.

priors = net.gmmnet.priors;
V_ = net.gmmnet.V_;
Sigma = net.gmmnet.Sigma1;
[N,D] = size(T);
SigmaNew = Sigma;
Norm = 1/D;

Phi = net.rbfnet.Phi;

[DimM,K] = size(Phi);

GM = gmdistribution(V_,Sigma,priors);

R = posterior(GM,T);

logL(1) = Compute_logL(GM,T);


W = net.gmmnet.W;
alpha = repmat(net.gmmnet.alpha,1,size(R,2));
invSigma = multinv(Sigma);
% G = sum(R,1);
T2 = repmat(T',1,1,K);
Phi2 = reshape(Phi,1,DimM,K);


for i=1:NIter
    % E-step;
    GM = gmdistribution(V_,SigmaNew,priors);
    R = posterior(GM,T);
    
    
    prodT = T'*R;
    
    % M-step
    disp(i)
    G = sum(R,1);
    invSigma = multiprod(invSigma,1./alpha,[1 2], 1);
    Kron2 = multiprod(invSigma,G, [1 2], 1);
    KronS = zeros(D*DimM,D*DimM);
    
    % W update (three different ways for computing the product in the 
    % equation. They should be equivalent). Suggested mem==1, that is quite
    % memory intensive but faster (epsecially on low-dimensional problems).
    
    if mem==0

        for j=1:DimM
            Psi_i = reshape(Phi(j,:).*Phi,DimM,1,1,K);
            KronS(D*(j-1)+1:D*(j-1)+D,:) = reshape(permute(sum(multiprod(Psi_i,...
                reshape(Kron2,1,D,D,K),[2 3],[2 3]),4),[2 3 1]),D,DimM*D);
        end
    elseif mem==1
        
        Psi_i = reshape(Phi,DimM,1,K);
        Const_Kron = permute(multiprod(Psi_i,reshape(Kron2,D*D,K),2,1),[3 2 1]);

        KronS = reshape(permute(reshape(multiprod(Phi,Const_Kron,[1 2],[1 2]),...
            DimM,D,D,DimM),[2 1 3 4]),DimM*D,DimM*D);
%         toc;
    else
        for k=1:K
             Psi = Phi(:,k)*Phi(:,k)';% + 1e-7.*eye(size(Phi,1));
             KronS = KronS + kron(Psi',Kron2(:,:,k));
        end
    end
    
    % Compute Right hand side of the W uodate equation
    R2 = reshape(R,N,1,K);
    rhs = sum(multiprod(multiprod(multiprod(invSigma,T2,[1 2],[1 2]),...
        R2,[1 2],[1 2]),Phi2,[1 2],[1 2]),3);


    % Invert the Kron matrix by svd, removing 
    [U,S,V] = svd(KronS);
    [cholDcmp, ~] = chol(KronS);
    s = diag(S);
    tolerance = 1e-4;
    p = sum(s>tolerance);

    
    Up = U(:,1:p);
    Vp = V(:,1:p);
    Inv = diag( 1./s(1:p));
    AInv = Vp * Inv * Up';
    x = AInv * rhs(:);

    W = reshape(x,D,DimM)';

    % Recompute the new centers of gaussian mixture.
    
    V_ = Phi'*W;
    
    % alpha update 
    
    d = mahal_d(T,V_,SigmaNew);
    
    M = R.*d;
    alpha = Norm.*sum(M,1)./G;%+ 1e-2.*ones(1,size(Sigma,3));
%     alpha = max([alpha',1e-2.*ones(size(alpha,2),1)],[],2);
%     alpha = ((Norm./size(G,1)).*sum(M,1))./G;
    alpha = min([max([alpha',1e-1.*ones(size(alpha,2),1)],[],2)...,
        ones(size(alpha,2),1).*2],[],2)';
%     alpha = alpha';
    SigmaNew = multiprod(Sigma,alpha,[1 2],1);


    % Symmetrize the covariance matrix.
    for j=1:size(SigmaNew,3)
        EigDec = eig(SigmaNew(:,:,j));
        if(~issymmetric(SigmaNew(:,:,j)) || any(EigDec<0))
            SigmaNew(:,:,j) = 0.5.*(SigmaNew(:,:,j) + SigmaNew(:,:,j)');
        end
    end 
    
    % Update the gaussian mixture model
    GMNew = gmdistribution(V_,SigmaNew,priors);

    [~,nlogL] = posterior(GMNew,T);
    logL(i) = -nlogL;
%     if i>2
%         if(abs(logL(i) - logL(i-1))< 1e-3)
%             break
%         end
%     end
end
net.gmmnet.alpha = alpha; 
net.gmmnet.V_ = V_;
net.gmmnet.W = W;
net.gmmnet.Sigma = SigmaNew;

net.gmmnet.GMdist = GMNew;