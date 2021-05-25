function [samples_output] = SAF_Mod(samples,original,iter_number,radius,mu)

% clc
% iter_n = 1;

% samples = datasample(original,fix(size(original,1)/2));

for iter = 1:iter_number
    disp(iter);
    
    %compute neighborhood

    
    [samples_original_neighbor_idx, samples_original_neighbor_dist] = rangesearch(original, samples, radius*2);
    [samples_samples_neighbor_idx, samples_samples_neighbor_dist] = rangesearch(samples, samples, radius*2);
    
    
 
        
    %Compute Average Term
    radius2 = radius^2;
    iradius16 = (-1)/(2*radius2);
    average_term = samples;
    repulsion_term = zeros(size(samples));
    
    for i = 1:size(samples)
        
        v_neighbors_idx = samples_original_neighbor_idx{i};
        v_neighbors_dist = samples_original_neighbor_dist{i};
        
        neighbor_size = length(v_neighbors_idx);
        if neighbor_size < 1
            continue;
        end
        
        if iter < 2 && i < 5
            disp(neighbor_size);
        end
        
        t = original(v_neighbors_idx,:);
        dist2 = v_neighbors_dist.^2;
        average_weight = exp(dist2.*iradius16);
        average_weight_sum = sum(average_weight);
        prov = average_weight * t;
        average = sum(prov,1);
%         Idx = find(average_weight_sum > 1e-20);
        if(average_weight_sum >1e-20)
            average_term(i,:) = average./average_weight_sum;
        end
        
%     end
    
    %Compute Repulsion Term
    
    %repulsion_term(:,:) = 0.0;
%     for i = 1:size(samples)
        v = samples(i, :);
        v_neighbors_idx = samples_samples_neighbor_idx{i};
        v_neighbors_dist = samples_samples_neighbor_dist{i};
        
        neighbor_size = length(v_neighbors_idx);
        if neighbor_size < 2
            continue;
        end
        
        
        t = samples(v_neighbors_idx(2:end), :);
        dist = v_neighbors_dist(2:end);
        dist2 = dist.^2;
        dist(dist<1e-8) = 1e-8;
        repulsion_weight = exp(dist2 .* iradius16) .* (1.0./dist).^1.0;
        diff = repmat(v,size(t,1),1) - t;
        repulsion_weight_sum = sum(repulsion_weight);
        repulsion = sum(diff.*repmat(repulsion_weight',1,3),1);
        if(repulsion_weight_sum >1e-20)
            repulsion_term(i,:) = repulsion./repulsion_weight_sum;
        end
        
%         repulsion_weight = 0.;
%         repulsion_weight_sum = 0.;
%         repulsion = zeros(size(v));
%         for j = 2:neighbor_size
%             t = samples(v_neighbors_idx(j), :);
%             dist = v_neighbors_dist(j);
%             dist2 = dist^2;
%             
%             if dist < 1e-8
%                 dist = 1e-8;
%                 continue;
%             end
%             
%             repulsion_weight = exp(dist2 * iradius16) * (1.0/dist)^1.0;
%             diff = v - t;
%             
%             repulsion_weight_sum = repulsion_weight_sum + repulsion_weight;
%             repulsion = repulsion + diff * repulsion_weight;
%         end
%         
%         if (repulsion_weight_sum > 1e-20)
%             repulsion_term(i, :) = repulsion / repulsion_weight_sum;
%         end
%         
    end
    
    samples = average_term + mu * repulsion_term;
    
    
           % figure result
    samples_output = samples;
end