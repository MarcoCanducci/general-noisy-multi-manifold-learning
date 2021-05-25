function [samples_output] = SAF_TrueOriginal(samples,original,iter_number,radius,mu)
copy_original = original;

% clc
tic;%tic1
t1=clock;

iter_n = 1;


for iter = 1:iter_number
    display(iter);
    
    %compute neighborhood
    samples_original_neighbor_idx = [];
    samples_original_neighbor_dist = [];
    samples_samples_neighbor_idx = [];
    samples_samples_neighbor_dist = [];
    
    [samples_original_neighbor_idx, samples_original_neighbor_dist] = rangesearch(original, samples, radius*3);
    [samples_samples_neighbor_idx, samples_samples_neighbor_dist] = rangesearch(samples, samples, radius*3);
    
    
 
        
    %Compute Average Term
    radius2 = radius^2;
    iradius16 = (-1)/(2*radius^2);
    average_term = samples;
    
    for i = 1:size(samples)
        
        v = samples(i, :);
        v_neighbors_idx = samples_original_neighbor_idx{i};
        v_neighbors_dist = samples_original_neighbor_dist{i};
        
        neighbor_size = length(v_neighbors_idx);
        if neighbor_size < 1
            continue;
        end
        
        if iter < 2 && i < 5
            display(neighbor_size);
        end
        
        average_weight = 0.;
        average_weight_sum = 0.;
        average = zeros(size(v));
        for j = 1:neighbor_size
            t = original(v_neighbors_idx(j), :);
            dist2 = v_neighbors_dist(j)^2;
            average_weight = exp(dist2 * iradius16);
            
            average_weight_sum = average_weight_sum + average_weight;
            average = average + t * average_weight;
        end
        
        if (average_weight_sum > 1e-20)
            average_term(i, :) = average / average_weight_sum;
        end
        
    end
    
    %Compute Repulsion Term
    repulsion_term = zeros(size(samples));
    %repulsion_term(:,:) = 0.0;
    for i = 1:size(samples)
        v = samples(i, :);
        v_neighbors_idx = samples_samples_neighbor_idx{i};
        v_neighbors_dist = samples_samples_neighbor_dist{i};
        
        neighbor_size = length(v_neighbors_idx);
        if neighbor_size < 2
            continue;
        end
        
        repulsion_weight = 0.;
        repulsion_weight_sum = 0.;
        repulsion = zeros(size(v));
        for j = 2:neighbor_size
            t = samples(v_neighbors_idx(j), :);
            dist = v_neighbors_dist(j);
            dist2 = dist^2;
            
            if dist < 1e-8
                dist = 1e-8;
                continue;
            end
            
            repulsion_weight = exp(dist2 * iradius16) * (1.0/dist)^1.0;
            diff = v - t;
            
            repulsion_weight_sum = repulsion_weight_sum + repulsion_weight;
            repulsion = repulsion + diff * repulsion_weight;
        end
        
        if (repulsion_weight_sum > 1e-20)
            repulsion_term(i, :) = repulsion / repulsion_weight_sum;
        end
        
    end
    
    samples = average_term + mu * repulsion_term;
    
    
           % figure result
%     samples_output = samples;
%     if mod(int16(iter), 5) == 1
% 
%         figure ('Position', [50 200 1020 420])
%         colors = [1,0,0; 0,0,0; 0,1,0;0,0,1;1,1,0;1,0,1;0,1,1;0,0,0;0.5,0,0;0,0.5,0;0,0,0.5;0.5,0.5,0;0,0.5,0.5;0.5,0,0.5;0.1,0.3,0.9;1,0,0;0,1,0;0,0,1;1,1,0;1,0,1;0,1,1;0,0,0;0.5,0,0;0,0.5,0;0,0,0.5;0.5,0.5,0;0,0.5,0.5;0.5,0,0.5;0.1,0.3,0.9];
% 
%         samples_dot_size = 12;
%         original_dot_size = 2;
%         fig_lim = 2.5;
% 
% 
% 
%         h=subplot(1, 2, 1);
%         hold on;
%         plot(copy_original(:,1),copy_original(:,2), '.','Color',colors(2,:),'MarkerSize',original_dot_size);
%         xlim([-fig_lim fig_lim]);
%         ylim([-fig_lim fig_lim]);  
%         hold off;
%         drawnow;
% 
% 
%         h=subplot(1, 2, 1);
%         hold on;
%         plot(samples_output(:,1),samples_output(:,2),'.','Color',colors(1,:),'MarkerSize',samples_dot_size);
%         xlim([-fig_lim fig_lim]);
%         ylim([-fig_lim fig_lim]);
%         hold off;
%         drawnow;
% 
% 
%         h=subplot(1, 2, 2);
%         hold on;
%         plot(samples_output(:,1),samples_output(:,2),'.','Color',colors(1,:),'MarkerSize',samples_dot_size);
%         xlim([-fig_lim fig_lim]);
%         ylim([-fig_lim fig_lim]);
%         hold off;
%         drawnow;
%         hold off
% 
%         %print(iter);
%         number = num2str(iter);
%         file_name = strcat('../results/iter_', number,'.jpeg');
%         print(gcf,'-r100','-djpeg', file_name);

%    end
end

samples_output = samples;

%close all;

end