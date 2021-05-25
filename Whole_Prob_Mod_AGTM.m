% function [TotMod,Data] = Whole_Prob_Mod(net,NoisyDataMan_Fin)
function [TotMod,Data] = Whole_Prob_Mod_AGTM(net,MixingCoeff)

MixingCoeff = MixingCoeff./sum(MixingCoeff);
s = 0;
s_D = 0;
for i=1:size(net,2)
    Nc = size(net{i}.gmmnet.V_,1);
    C = net{i}.gmmnet.Sigma;
    mu = net{i}.gmmnet.V_;    
    C_Tot{i} = C;
    mu_Tot{i} = mu;
    
    s = s + Nc;
end

C_Fin = zeros(size(C,1), size(C,1), s);
mu_Fin = zeros(s,size(C,1));
prior = zeros(s,1);
Data = zeros(s_D,size(C,1));

In_S = 1;
% In_D = 1;
for i=1:size(mu_Tot,2)
    C_Fin(:,:,In_S:In_S + size(C_Tot{i},3) - 1) = C_Tot{i};
    mu_Fin(In_S:In_S + size(mu_Tot{i},1) - 1,:) = mu_Tot{i};
    prior(In_S:In_S + size(mu_Tot{i},1) - 1) = ones(size(mu_Tot{i},1),1).*MixingCoeff(i);
    
    In_S = In_S + size(C_Tot{i},3);    
%     Data(In_D:In_D + size(NoisyDataMan_Fin{i},1) - 1,:) = NoisyDataMan_Fin{i};
%     In_D = In_D + size(NoisyDataMan_Fin{i},1);
%     disp(i)
end
    
    

% prior = ones(s,1).*1/s;

TotMod = gmdistribution(mu_Fin,C_Fin,prior);




